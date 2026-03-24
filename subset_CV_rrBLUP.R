
### Subset-population cross-validation (rrBLUP)####

library(rrBLUP)
library(writexl)

set.seed(123)

format_geno <- function(file) {
  geno <- read.delim(file, header = TRUE)
  geno[geno == "0"] <- -1
  geno[geno == "1"] <- 0
  geno[geno == "2"] <- 1
  rownames(geno) <- geno$taxa
  geno <- geno[, -1, drop = FALSE]
  geno[] <- lapply(geno, function(x) as.numeric(as.character(x)))
  return(geno)
}

format_pheno <- function(file) {
  pheno <- read.delim(file, header = TRUE)
  rownames(pheno) <- pheno$Genotype
  pheno <- pheno[, "BLUP", drop = FALSE]
  return(pheno)
}

impute_mean <- function(df) {
  df[] <- lapply(df, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x
  })
  return(df)
}

# Data files ####

geno_2R  <- format_geno("2row_geno")
geno_6R  <- format_geno("6row_geno")
geno_germ <- format_geno("genebank_geno")

pheno_2R  <- format_pheno("2row_pheno")
pheno_6R  <- format_pheno("6row_pheno")
pheno_germ <- format_pheno("genebank_pheno")


# Find common makrers####

common_markers <- Reduce(intersect, list(
  colnames(geno_2R),
  colnames(geno_6R),
  colnames(geno_germ)
))

geno_2R  <- geno_2R[, common_markers]
geno_6R  <- geno_6R[, common_markers]
geno_germ <- geno_germ[, common_markers]


#CV scenario:
# Training = all populations
# Validation target = genebank

geno_train_pop <- rbind(geno_2R, geno_6R, geno_germ)
pheno_train_pop <- rbind(pheno_2R, pheno_6R, pheno_germ)

geno_val_pop <- geno_germ
pheno_val_pop <- pheno_germ

val_prop <- 0.2
n_runs <- 20


# Run CV ####

all_GEBV <- data.frame()
PA_values <- numeric(n_runs)

for (i in 1:n_runs) {
  
  cat("Run:", i, "\n")
  
  val_ids <- sample(
    rownames(geno_val_pop),
    size = round(val_prop * nrow(geno_val_pop))
  )
  
  train_ids <- setdiff(rownames(geno_train_pop), val_ids)
  
  geno_train <- geno_train_pop[train_ids, , drop = FALSE]
  geno_valid <- geno_val_pop[val_ids, , drop = FALSE]
  
  pheno_train <- pheno_train_pop[train_ids, , drop = FALSE]
  pheno_valid <- pheno_val_pop[val_ids, , drop = FALSE]
  
  # Impute
  geno_train <- impute_mean(geno_train)
  geno_valid <- impute_mean(geno_valid)
  
  # Model
  model <- mixed.solve(y = pheno_train$BLUP, Z = as.matrix(geno_train))
  
  X_valid <- cbind(1, as.matrix(geno_valid))
  u <- c(model$beta, model$u)
  pred <- as.vector(X_valid %*% u)
  
  obs <- pheno_valid$BLUP
  
  PA <- if (sd(obs) > 0 && sd(pred) > 0) cor(obs, pred) else NA
  PA_values[i] <- PA
  
  temp <- data.frame(
    Run = i,
    Genotype = rownames(geno_valid),
    Observed = obs,
    Predicted = pred
  )
  
  all_GEBV <- rbind(all_GEBV, temp)
}

# Save results ####

PA_summary <- data.frame(
  Run = 1:n_runs,
  Predictive_Ability = PA_values
)

write_xlsx(
  list(PA = PA_summary, GEBV = all_GEBV),
  "subset_CV_results.xlsx"
)

cat("Mean PA:", mean(PA_values, na.rm = TRUE), "\n")