# Between-population cross-validation (rrBLUP) ####

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

geno_TS <- format_geno("2row_geno")
pheno_TS <- format_pheno("2row_pheno")

geno_BS <- format_geno("6row_geno")
pheno_BS <- format_pheno("6row_pheno")

# Other options:
# "genebank_geno" / "genebank_pheno"
# "2row_6row_geno" / "2row_6row_pheno"
# "all_geno" / "all_pheno"
# "2row_genebank_geno" / "2row_genebank_pheno"
# "6row_genebank_geno" / "6row_genebank_pheno"

# Match genotype and phenotype ####

common_TS_ids <- intersect(rownames(geno_TS), rownames(pheno_TS))
geno_TS <- geno_TS[common_TS_ids, ]
pheno_TS <- pheno_TS[common_TS_ids, , drop = FALSE]

common_BS_ids <- intersect(rownames(geno_BS), rownames(pheno_BS))
geno_BS <- geno_BS[common_BS_ids, ]
pheno_BS <- pheno_BS[common_BS_ids, , drop = FALSE]

# Find common makrers####

common_markers <- intersect(colnames(geno_TS), colnames(geno_BS))
geno_TS <- geno_TS[, common_markers]
geno_BS <- geno_BS[, common_markers]


geno_TS <- impute_mean(geno_TS)
geno_BS <- impute_mean(geno_BS)


geno_combined <- rbind(geno_TS, geno_BS)
geno_mat <- as.matrix(geno_combined)

maf <- apply(geno_mat, 2, function(x) {
  x <- x[!is.na(x)]
  p <- (mean(x) + 1) / 2
  min(p, 1 - p)
})

maf_threshold <- 0.05
valid_markers <- names(maf[maf > 0 & maf >= maf_threshold])

geno_TS <- geno_TS[, valid_markers]
geno_BS <- geno_BS[, valid_markers]

cat("Markers retained:", length(valid_markers), "\n")

# Run CV ####

model <- mixed.solve(y = pheno_TS$BLUP, Z = as.matrix(geno_TS))


X <- cbind(1, as.matrix(geno_BS))
u <- c(model$beta, model$u)
pred <- as.vector(X %*% u)

# Save results ####

results <- data.frame(
  Genotype = rownames(geno_BS),
  Observed = pheno_BS$BLUP,
  Predicted = pred
)

write_xlsx(results, "between_CV_results.xlsx")