
# WITHIN-POPULATION CROSS-VALIDATION USING RR-BLUP #
# Random 80% training / 20% validation, repeated 500 times #

library(rrBLUP)
library(openxlsx)

set.seed(123)

# 1. LOAD DATA #

genotype_data <- read.delim("genotypic_data.txt", header = TRUE, check.names = FALSE)
phenotype_data <- read.delim("whole_phenotype.txt", header = TRUE, sep = "\t", na.strings = "NA")

# 2. FORMAT GENOTYPE DATA #

genotype_data[genotype_data == "0"] <- -1
genotype_data[genotype_data == "1"] <- 0
genotype_data[genotype_data == "2"] <- 1

rownames(genotype_data) <- genotype_data$taxa
genotype_data <- genotype_data[, -1, drop = FALSE]

genotype_data[] <- lapply(
  genotype_data,
  function(x) as.numeric(as.character(x)))

# 3. FORMAT PHENOTYPE DATA #

rownames(phenotype_data) <- phenotype_data$Genotype
phenotype_data <- phenotype_data[, "BLUP", drop = FALSE]

# Match genotype and phenotype IDs

common_ids <- intersect(rownames(genotype_data), rownames(phenotype_data))
genotype_data <- genotype_data[common_ids,,drop = FALSE]
phenotype_data <- phenotype_data[common_ids,,drop = FALSE]

if (any(!is.finite(phenotype_data$BLUP))) {
  stop("Non-finite BLUP values detected.")
}

if (anyNA(genotype_data)) {
  stop("Missing values detected in genotype matrix.")
}

# 4. CROSS-VALIDATION SETTINGS #
  
runs <- 500
train_prop <- 0.80

all_ids <- rownames(genotype_data)
train_size <- round(length(all_ids) * train_prop)

res <- data.frame()
marker_results <- data.frame()
validation_results <- vector("list", runs)

# 5. RUN CROSS-VALIDATION #

for (r in seq_len(runs)) {

  # Random 80% training / 20% validation split

  train_ids <- sample(all_ids, train_size)
  test_ids <- setdiff(all_ids, train_ids)

  geno_train <- genotype_data[rain_ids, ,drop = FALSE]
  geno_test <- genotype_data[test_ids,,drop = FALSE]

  pheno_train <- phenotype_data[train_ids,,drop = FALSE]
  pheno_test <- phenotype_data[test_ids, ,drop = FALSE]

  # Remove zero-variance markers based on training set

  marker_sd <- apply(geno_train,2,sd,na.rm = TRUE)
  keep_markers <- is.finite(marker_sd) & marker_sd > 0

  geno_train <- geno_train[, keep_markers,drop = FALSE]
  geno_test <- geno_test[,keep_markers,drop = FALSE]

  # Fit RR-BLUP model

  model <- mixed.solve(y = pheno_train$BLUP,Z = as.matrix(geno_train))

  # Predict validation set

  X_test <- cbind(1,as.matrix(geno_test))
  u_effects <- c(model$beta,model$u)
  pred <- as.vector(X_test %*% u_effects)
  obs <- pheno_test$BLUP

  # Calculate predictive ability

  pa <- if (sd(obs) > 0 && sd(pred) > 0) {cor(obs, pred) } else {NA}

  # Save results

  res <- rbind(res, data.frame(Run = r, Predictive_Ability = pa))

  marker_results <- rbind(marker_results, data.frame(Run = r,Starting_Markers = ncol(genotype_data), 
                  Zero_Variance_Removed = sum(!keep_markers), Markers_Used = sum(keep_markers)))

  validation_results[[r]] <- data.frame(
    Run = r,Genotype = rownames(geno_test),
    Observed = obs, Predicted = pred)}

# 6. SUMMARIZE RESULTS #

all_validation_results <- do.call(
  rbind, validation_results)

overall_summary <- data.frame(
  Total_Genotypes = nrow(genotype_data),
  Starting_Markers = ncol(genotype_data),
  Repetitions = runs,
  Training_Proportion = train_prop,
  Mean_Markers_Used = mean(marker_results$Markers_Used),
  Minimum_Markers_Used = min(marker_results$Markers_Used),
  Maximum_Markers_Used = max(marker_results$Markers_Used),
  Mean_PA = mean(res$Predictive_Ability, na.rm = TRUE),
  SD_PA = sd(res$Predictive_Ability, na.rm = TRUE)
)

# 7. SAVE RESULTS #

wb <- createWorkbook()

addWorksheet(wb, "Overall_Summary")
addWorksheet(wb, "Predictive_Abilities")
addWorksheet(wb, "Marker_Counts")
addWorksheet(wb, "Validation_GEBVs")

writeData(wb, "Overall_Summary", overall_summary)
writeData(wb, "Predictive_Abilities", res)
writeData(wb, "Marker_Counts", marker_results)
writeData(wb, "Validation_GEBVs", all_validation_results)

saveWorkbook( wb, "whole_population_within_CV_500runs.xlsx",
  overwrite = TRUE)
