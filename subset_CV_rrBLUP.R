# SUBSET-POPULATION CROSS-VALIDATION

library(rrBLUP)
library(writexl)

set.seed(123)

# 1. LOAD ALL DATA

geno_2R <- read.delim("2row_geno.txt",header = TRUE,check.names = FALSE)
geno_6R <- read.delim("6row_geno.txt",header = TRUE,check.names = FALSE)

pheno_2R <- read.delim("2row_pheno.txt",header = TRUE,sep = "\t",na.strings = "NA")

pheno_6R <- read.delim("6row_pheno.txt",header = TRUE,sep = "\t", na.strings = "NA")

# 2. FORMAT GENOTYPE DATA

format_geno <- function(geno) {
  
  geno[geno == "0"] <- -1
  geno[geno == "1"] <- 0
  geno[geno == "2"] <- 1
  
  rownames(geno) <- geno$taxa
  geno <- geno[, -1, drop = FALSE]
  
  # Ensure marker data are numeric
  geno[] <- lapply(
    geno,
    function(x) as.numeric(as.character(x))
  ) return(geno)
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
  )
)

geno_2R <- geno_2R[
  ,common_markers,drop = FALSE]

geno_6R <- geno_6R[
  ,common_markers,drop = FALSE]


cat(
  "Starting common markers:",
  length(common_markers),
  "\n")

# Check for missing genotype values

if (
  anyNA(geno_2R) ||
  anyNA(geno_6R) ||
  anyNA(geno_germ)
) {
  stop("Missing values ")
}

# Check phenotype values

if (
  any(!is.finite(pheno_2R$BLUP)) ||
  any(!is.finite(pheno_6R$BLUP))
) {
  stop("Non-finite BLUP values detected.")
}

# 5. CHOOSE TRAINING PANEL AND TARGET POPULATION

# Breeding panel (2R + 6R) -> 2R breeding
#
# In each run:
# Training = all 6R + 80% of 2R
# Validation = remaining 20% of 2R

geno_train_pop <- rbind(
  geno_2R,
  geno_6R)

pheno_train_pop <- rbind(
  pheno_2R,
  pheno_6R)

geno_val_pop <- geno_2R
pheno_val_pop <- pheno_2R

scenario_name <- "Breeding Panel_to_2R"

# 6. MATCH GENOTYPE AND PHENOTYPE IDs

train_common_ids <- intersect(
  rownames(geno_train_pop),
  rownames(pheno_train_pop)
)

geno_train_pop <- geno_train_pop[
  train_common_ids,
  ,
  drop = FALSE
]

pheno_train_pop <- pheno_train_pop[
  train_common_ids,
  ,
  drop = FALSE
]

val_common_ids <- intersect(
  rownames(geno_val_pop),
  rownames(pheno_val_pop)
)

geno_val_pop <- geno_val_pop[
  val_common_ids,
  ,
  drop = FALSE
]

pheno_val_pop <- pheno_val_pop[
  val_common_ids,
  ,
  drop = FALSE
]

if (
  !all(
    rownames(geno_val_pop) %in%
    rownames(geno_train_pop)
  )
) {
  stop(
    "Target population is not in the training panel."
  )
}

# 7. CROSS-VALIDATION SETTINGS

val_prop <- 0.20
n_runs <- 500

# 8. PREPARE RESULT OBJECTS

PA_values <- numeric(n_runs)
all_GEBV <- vector("list",n_runs)

marker_results <- data.frame()
training_size_results <- data.frame()

# 9. RUN SUBSET-POPULATION CROSS-VALIDATION

for (i in seq_len(n_runs)) {
  
  # Randomly select 20% of TARGET population for validation
  
  val_ids <- sample(
  rownames(geno_val_pop),
    size = round(
    val_prop * nrow(geno_val_pop)
    ))
  
  train_ids <- setdiff(
  rownames(geno_train_pop),
   val_ids)
    
  # Construct training data
 
  geno_train <- geno_train_pop[
  train_ids,,drop = FALSE]
  
  pheno_train <- pheno_train_pop[
  train_ids,,drop = FALSE]
  
  # Construct validation data
  
  geno_valid <- geno_val_pop[
    val_ids,,drop = FALSE]
  
  pheno_valid <- pheno_val_pop[
    val_ids,,drop = FALSE]
  
  # --------------------------------------------------------
  # Remove zero-variance markers
  # BASED ON TRAINING SET ONLY
  # --------------------------------------------------------
  
  marker_sd <- apply(
    geno_train,2,sd,na.rm = TRUE)
  
  keep_markers <- (is.finite(marker_sd) & marker_sd > 0)
  
  markers_removed <- sum(!keep_markers)
  
  markers_used <- sum(keep_markers)
  
  # Use exactly the same retained markers
  # in training and validation
  
  geno_train <- geno_train[,
    keep_markers,drop = FALSE]
  
  geno_valid <- geno_valid[
    , keep_markers,drop = FALSE]
  
 
  #  check
  
  if (
    ncol(geno_train) !=
    ncol(geno_valid)
  ) {
    stop(
      paste(
        "Marker mismatch in run",
        i ))
  }
  
  # FIT RR-BLUP
  
  pmodel <- mixed.solve(
    y = pheno_train$BLUP,
    Z = as.matrix(geno_train))
  
  # PREDICT VALIDATION SET
  
  X_valid <- cbind(
    1,as.matrix(geno_valid))
  
  u_effects <- c(
    pmodel$beta,
    pmodel$u)
  
  y_pred <- as.vector(
    X_valid %*% u_effects)
 
  # PREDICTIVE ABILITY
  
  obs <- pheno_valid$BLUP
  
  PA <- if (
    sd(obs) > 0 &&
    sd(y_pred) > 0
  ) {
    cor(
      obs,
      y_pred
    )
  } else {
    NA
  }
  
  PA_values[i] <- PA
  
  # STORE MARKER INFORMATION
  
  marker_results <- rbind(
    marker_results,
    data.frame(
      Run = i,
      Starting_Markers = length(common_markers),
      Zero_Variance_Removed = markers_removed,
      Markers_Used = markers_used
    )
  )
 
  # STORE TRAINING / VALIDATION SIZE
  
  training_size_results <- rbind(
    training_size_results,
    data.frame(
      Run = i,
      Training_N = nrow(geno_train),
      Validation_N = nrow(geno_valid)
    )
  )
  
  # --------------------------------------------------------
  # STORE PREDICTIONS
  # --------------------------------------------------------
  
  all_GEBV[[i]] <- data.frame(
    Run = i,
    Genotype = rownames(geno_valid),
    Observed = obs,
    Predicted = y_pred
  )
  
# 10. COMBINE STORED PREDICTIONS

all_GEBV <- do.call(
  rbind,
  all_GEBV
)

# 11. SUMMARY STATISTICS

mean_PA <- mean(
  PA_values,
  na.rm = TRUE
)

sd_PA <- sd(
  PA_values,
  na.rm = TRUE
)

se_PA <- sd_PA / sqrt(
  sum(
    !is.na(PA_values)
  )
)

mean_markers <- mean(
  marker_results$Markers_Used
)

min_markers <- min(
  marker_results$Markers_Used
)

max_markers <- max(
  marker_results$Markers_Used
)

mean_zero_removed <- mean(
  marker_results$Zero_Variance_Removed
)

PA_summary <- data.frame(
  Run = seq_len(n_runs),
  Predictive_Ability = PA_values
)

Overall_summary <- data.frame(
  Scenario = scenario_name,
  Training_Panel_N = nrow(geno_train_pop),
  Target_Population_N = nrow(geno_val_pop),
  Validation_Proportion = val_prop,
  Repetitions = n_runs,
  Starting_Markers = length(common_markers),
  Mean_Zero_Variance_Removed = mean_zero_removed,
  Mean_Markers_Used = mean_markers,
  Minimum_Markers_Used = min_markers,
  Maximum_Markers_Used = max_markers,
  Mean_PA = mean_PA,
  SD_PA = sd_PA,
  SE_PA = se_PA
)

# 12. SAVE RESULTS

output_name <- paste0(
  scenario_name,
  "_subset_CV_500runs.xlsx"
)

write_xlsx(
  list(Overall_Summary = Overall_summary,
    PA_per_Run = PA_summary,
    Marker_Counts = marker_results,
    Sample_Sizes = training_size_results,
    All_GEBV = all_GEBV
  ),
  output_name)
