
# Within-population cross-validation (rrBLUP) ####

library(openxlsx)
library(rrBLUP)

set.seed(123)

# Data files ####

genotype_file <- "2row_geno"
phenotype_file <- "2row_pheno"

# Other options:
# "6row_geno" / "6row_pheno"
# "genebank_geno" / "genebank_pheno"
# "2row_6row_geno" / "2row_6row_pheno"
# "all_geno" / "all_pheno"

# Data files_1 ####

genotype_data <- read.table(genotype_file, header = TRUE, check.names = FALSE)
phenotype_data <- read.table(phenotype_file, header = TRUE, check.names = FALSE)

# Genotypic data ####

genotype_data[genotype_data == "0"] <- -1
genotype_data[genotype_data == "1"] <- 0
genotype_data[genotype_data == "2"] <- 1

rownames(genotype_data) <- genotype_data$taxa
genotype_data <- genotype_data[, setdiff(names(genotype_data), "taxa"), drop = FALSE]

genotype_data[] <- lapply(genotype_data, function(x) as.numeric(as.character(x)))

# Match genotype and phenotype ####

rownames(phenotype_data) <- phenotype_data$Genotype

common_ids <- intersect(rownames(phenotype_data), rownames(genotype_data))
phenotype_data <- phenotype_data[common_ids, , drop = FALSE]
genotype_data  <- genotype_data [common_ids, , drop = FALSE]

if (!"BLUP" %in% names(phenotype_data)) stop("Missing BLUP column")
if (any(!is.finite(phenotype_data$BLUP))) stop("Non-finite BLUP values")


impute_col_mean <- function(v) {
  mu <- mean(v, na.rm = TRUE)
  if (!is.finite(mu)) mu <- 0
  v[!is.finite(v)] <- mu
  v
}

genotype_data[] <- lapply(genotype_data, impute_col_mean)


# Remove zero-varianec markers ####

marker_sd <- vapply(genotype_data, sd, numeric(1))
genotype_data <- genotype_data[, marker_sd > 0, drop = FALSE]


pheno_geno <- cbind(phenotype_data, genotype_data)
pheno_geno$Genotype <- factor(pheno_geno$Genotype)

markernames <- colnames(genotype_data)

# Run CV ####

runs <- 500
train_prop <- 0.8

res <- NULL
validation_results <- vector("list", runs)


all_ids <- unique(as.character(pheno_geno$Genotype))
train_size <- round(length(all_ids) * train_prop)

for (r in seq_len(runs)) {
  
  cat("Run:", r, "\n")
  
  train_ids <- sample(all_ids, train_size)
  test_ids  <- setdiff(all_ids, train_ids)
  
  train_data <- pheno_geno[pheno_geno$Genotype %in% train_ids, ]
  test_data  <- pheno_geno[pheno_geno$Genotype %in% test_ids, ]
  
  y_train <- as.numeric(train_data$BLUP)
  
  Z_train <- as.matrix(train_data[, markernames])
  Z_test  <- as.matrix(test_data[, markernames])
  
  model <- mixed.solve(y = y_train, Z = Z_train)
  
  X <- cbind(1, Z_test)
  u <- c(model$beta, model$u)
  pred <- as.vector(X %*% u)
  
  obs <- as.numeric(test_data$BLUP)
  
  pa <- if (sd(obs) > 0 && sd(pred) > 0) cor(obs, pred) else NA
  
  res <- rbind(res, data.frame(Run = r, Predictive_Ability = pa))
  
  validation_results[[r]] <- data.frame(
    Run = r,
    Genotype = test_data$Genotype,
    Observed = obs,
    Predicted = pred
  )
}


# Save results ####

all_validation_results <- do.call(rbind, validation_results)

wb <- createWorkbook()
addWorksheet(wb, "Predictive_Abilities")
addWorksheet(wb, "Validation_GEBVs")

writeData(wb, "Predictive_Abilities", res)
writeData(wb, "Validation_GEBVs", all_validation_results)

saveWorkbook(wb, "within_CV_results.xlsx", overwrite = TRUE)

cat("Mean PA:", mean(res$Predictive_Ability, na.rm = TRUE), "\n")