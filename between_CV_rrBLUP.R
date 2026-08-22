# BETWEEN-POPULATION CROSS-VALIDATION

library(rrBLUP)
library(writexl)

# 1. LOAD DATA

geno_2R <- read.delim("2row_geno.txt", header = TRUE)
geno_6R <- read.delim("6row_geno.txt", header = TRUE)

pheno_2R <- read.delim("2row_pheno.txt",header = TRUE,sep = "\t",na.strings = "NA")
pheno_6R <- read.delim("6row_pheno.txt",header = TRUE,sep = "\t",na.strings = "NA")

# 2. FORMAT GENOTYPE DATA

format_geno <- function(geno) {
  geno[geno == "0"] <- -1
  geno[geno == "1"] <- 0
  geno[geno == "2"] <- 1
  
  rownames(geno) <- geno$taxa
  geno <- geno[, -1, drop = FALSE]
  
  # Ensure marker data are numeric
  geno[] <- lapply(geno, function(x) as.numeric(as.character(x)))
  
  return(geno)
}

geno_2R <- format_geno(geno_2R)
geno_6R <- format_geno(geno_6R)

# 3. FORMAT PHENOTYPE DATA

format_pheno <- function(pheno) {
  rownames(pheno) <- pheno$Genotype
  pheno <- pheno[, "BLUP", drop = FALSE]
  return(pheno)
}

pheno_2R <- format_pheno(pheno_2R)
pheno_6R <- format_pheno(pheno_6R)

# 4. KEEP COMMON MARKERS

common_markers <- Reduce(
  intersect,
  list(
    colnames(geno_2R),
    colnames(geno_6R)
  ))

geno_2R <- geno_2R[, common_markers, drop = FALSE]
geno_6R <- geno_6R[, common_markers, drop = FALSE]

cat("Starting common markers:", length(common_markers), "\n")

# 5. CHOOSE TRAINING AND VALIDATION POPULATIONS

geno_train <- geno_6R
pheno_train <- pheno_6R

geno_valid <- geno_2R
pheno_valid <- pheno_2R

scenario_name <- "6R_to_2R"

# 6. MATCH GENOTYPE AND PHENOTYPE IDs
                   
train_ids <- intersect(
  rownames(geno_train),
  rownames(pheno_train))

valid_ids <- intersect(
  rownames(geno_valid),
  rownames(pheno_valid))

geno_train <- geno_train[
  train_ids,
  ,drop = FALSE]

pheno_train <- pheno_train[
  train_ids,
  ,drop = FALSE]

geno_valid <- geno_valid[
  valid_ids,
  ,drop = FALSE]

pheno_valid <- pheno_valid[
  valid_ids,
  , drop = FALSE]

# 7. CHECK DATA

cat("\nScenario:", scenario_name, "\n")
cat("Training genotypes:", nrow(geno_train), "\n")
cat("Validation genotypes:", nrow(geno_valid), "\n")
cat("Starting markers:", ncol(geno_train), "\n")

if (
  anyNA(geno_train) ||
  anyNA(geno_valid)
) {
  stop("Missing values")
}

# 8. REMOVE ZERO-VARIANCE MARKERS
#    BASED ON TRAINING POPULATION ONLY

marker_sd <- apply(geno_train,2, sd,na.rm = TRUE)
keep_markers <- is.finite(marker_sd) & marker_sd > 0

cat(
  "Zero-variance markers removed:",
  sum(!keep_markers),
  "\n")

cat(
  "Markers retained:",
  sum(keep_markers),
  "\n")

# Use exactly the same retained markers in TS and VS

geno_train <- geno_train[,keep_markers,drop = FALSE]
geno_valid <- geno_valid[, keep_markers,drop = FALSE]

# 9. FIT RR-BLUP MODEL

pmodel <- mixed.solve(
  y = pheno_train$BLUP,
  Z = as.matrix(geno_train))

# 10. PREDICT VALIDATION POPULATION

X_valid <- cbind(1,as.matrix(geno_valid))
u_effects <- c(pmodel$beta,pmodel$u)
y_pred <- X_valid %*% u_effects

# 11. PREDICTIVE ABILITY

PA <- cor(pheno_valid$BLUP,y_pred,use = "pairwise.complete.obs")

cat("\nPredictive Ability:", round(PA, 3), "\n")

# 12. SAVE RESULTS
                   
results <- data.frame(
  Genotype = rownames(geno_valid),
  Observed = pheno_valid$BLUP,
  Predicted = as.vector(y_pred))

summary_results <- data.frame(
  Scenario = scenario_name,
  Training_N = nrow(geno_train),
  Validation_N = nrow(geno_valid),
  Starting_Markers = length(common_markers),
  Zero_Variance_Removed = sum(!keep_markers),
  Markers_Used = ncol(geno_train),
  Predictive_Ability = PA)

output_name <- paste0(scenario_name,"_between_population_CV.xlsx")
write_xlsx(list(Summary = summary_results,Predictions = results),output_name)
