# ==== Load Required Libraries ====
library(data.table)
library(emdi)
library(SAEforest)   # (if needed for other functions; MERFranger is assumed to be available)
library(lme4)
library(reshape2)

# ==== Modified Generic predict Function for MERFranger ====
# This function now requires a 'newdata' argument so that non-quantile predictions
# are provided with a data argument.
predict.MERFranger <- function(object, newdata, ...) {
  # Predict using the fixed effect part (random forest) and the random effects model
  retval <- predict(object$Forest, data = newdata, ...)$predictions +
    predict(object$EffectModel, newdata = newdata, allow.new.levels = TRUE, ...)
  return(retval)
}

# ==== Step 1: Load Data ====
# Use an interactive file chooser to load your population data.
pop <- fread(file.choose())  
pop[, Domain := as.integer(as.character(Domain))]
colnames(pop)
# ==== Step 2: Set Parameters ====
B <- 500                       # Number of simulations
y_col <- "Y_Inter_Par"         # Study variable
predictors <- c("X1_IPAR", "X2_IPAR", paste0("V_", 1:10))

pop$Y_Normal_Par

# Precompute true domain means and population sizes
true_vals <- pop[, .(True_Mean = mean(get(y_col))), by = Domain]
pop_size  <- pop[, .(N_pop = .N), by = Domain]

# ==== Prepare Lists to Store Results ====
direct_estimates_list <- list()
lmm_estimates_list    <- list()
merf_estimates_list   <- list()
sample_counts_list    <- list()

# Define MERF configurations (for MERF1, MERF2, and MERF3)
merf_cfgs <- list(
  list(Model = "MERF1", num.trees = 500,  MaxIterations = 100, ErrorTolerance = 0.001, 
       mtry = 8, min.node.size = 10, sample.fraction = 0.7, B_adj = 50),
  list(Model = "MERF2", num.trees = 750,  MaxIterations = 150, ErrorTolerance = 0.0008,
       mtry = 10, min.node.size = 5,  sample.fraction = 0.8, B_adj = 100),
  list(Model = "MERF3", num.trees = 1000, MaxIterations = 200, ErrorTolerance = 0.0005, 
       mtry = 12, min.node.size = 3,  sample.fraction = 0.9, B_adj = 150)
)

# ==== Step 3: Run Simulation ====
for (b in 1:B) {
  set.seed(100 + b)
  sampled_data <- pop[sample(.N, 1200)]
  
  # Record sample counts by Domain
  sc <- sampled_data[, .(Sample_Count = .N), by = Domain]
  sc[, Simulation := b]
  sample_counts_list[[b]] <- sc
  
  # Ensure Domain column is integer in the sample
  sampled_data[, Domain := as.integer(as.character(Domain))]
  
  ## --- 3.1 Direct Estimation ---
  direct_obj <- direct(y = y_col, 
                       smp_data = sampled_data, 
                       smp_domains = "Domain", 
                       weights = "weight", 
                       var = TRUE)
  direct_df <- as.data.table(estimators(direct_obj, indicator = "Mean", MSE = TRUE, CV = TRUE))
  setnames(direct_df, old = c("Mean", "Mean_MSE", "Mean_CV"),
           new = c("Mean", "MSE", "CV"))
  direct_df[, Domain := as.integer(as.character(Domain))]
  direct_df <- merge(direct_df, true_vals, by = "Domain", all.x = TRUE)
  direct_df[, `:=`(Simulation = b, Model = "Direct")]
  direct_estimates_list[[b]] <- direct_df[, .(Domain, Mean, MSE, CV, True_Mean, Model, Simulation)]
  
  ## --- 3.2 LMM Estimation --- GREG Generalized Regression....Align with Battesse harter and fuler BHF model.
  pop_lmm <- copy(pop)
  sampled_data_lmm <- copy(sampled_data)
  fml <- as.formula(paste(y_col, "~", paste(predictors, collapse = "+"), "+ (1|Domain)"))
  lmm_fit <- lmer(fml, data = sampled_data_lmm)
  
  pop_lmm[, pred_lmm := predict(lmm_fit, newdata = pop_lmm, allow.new.levels = TRUE)]
  sampled_data_lmm[, pred_lmm := predict(lmm_fit, newdata = sampled_data_lmm)]
  
  pop_sum_lmm <- pop_lmm[, .(PopSum = sum(pred_lmm)), by = Domain]
  adj_lmm <- sampled_data_lmm[, .(Sample_Adjustment = sum(weight * (get(y_col) - pred_lmm))), by = Domain]
  
  lmm_df <- merge(pop_sum_lmm, adj_lmm, by = "Domain")
  lmm_df <- merge(lmm_df, pop_size, by = "Domain")
  lmm_df <- merge(lmm_df, true_vals, by = "Domain")
  lmm_df[, MA_Est := (PopSum + Sample_Adjustment) / N_pop]
  lmm_df[, `:=`(Model = "LMM", Simulation = b)]
  lmm_estimates_list[[b]] <- lmm_df[, .(Domain, MA_Est, Model, Simulation, True_Mean)]
  
  ## --- 3.3 MERF Estimation using MERFranger ---
  for (cfg in merf_cfgs) {
    merf_result <- tryCatch({
      pop_merf <- copy(pop)
      sampled_data_merf <- copy(sampled_data)
      
      # Fit the MERFranger model
      mrf <- MERFranger(
        Y = sampled_data_merf[[y_col]],
        X = sampled_data_merf[, ..predictors],
        random = "(1|Domain)",
        data = sampled_data_merf,
        num.trees = cfg$num.trees,
        MaxIterations = cfg$MaxIterations,
        ErrorTolerance = cfg$ErrorTolerance,
        B_adj = cfg$B_adj,
        mtry = cfg$mtry,
        min.node.size = cfg$min.node.size,
        sample.fraction = cfg$sample.fraction,
        importance = "impurity",
        na.rm = TRUE
      )
      
      # Obtain predictions on the population and sample data using the new predict function:
      pop_merf[, pred_merf := predict(mrf, newdata = pop_merf)]
      sampled_data_merf[, pred_merf := predict(mrf, newdata = sampled_data_merf)]
      
      # Aggregate predictions at domain level
      pop_sum_merf <- pop_merf[, .(PopSum = sum(pred_merf)), by = Domain]
      adj_merf <- sampled_data_merf[, .(Sample_Adjustment = sum(weight * (get(y_col) - pred_merf))), by = Domain]
      
      est_df <- merge(pop_sum_merf, adj_merf, by = "Domain")
      est_df <- merge(est_df, pop_size, by = "Domain")
      est_df <- merge(est_df, true_vals, by = "Domain")
      est_df[, MA_Est := (PopSum + Sample_Adjustment) / N_pop]
      est_df[, `:=`(Model = cfg$Model, Simulation = b)]
      
      est_df[, .(Domain, MA_Est, Model, Simulation, True_Mean)]
      
    }, error = function(e) {
      message("MERF model failed for ", cfg$Model, " in simulation ", b, ": ", e$message)
      data.table(Domain = unique(pop$Domain),
                 MA_Est = NA_real_,
                 Model = cfg$Model,
                 Simulation = b,
                 True_Mean = NA_real_)
    })
    merf_key <- paste0(cfg$Model, "_", b)
    merf_estimates_list[[merf_key]] <- merf_result
  }
}

# ---- Remove any NULL elements from lists ----
direct_estimates_list <- Filter(Negate(is.null), direct_estimates_list)
lmm_estimates_list    <- Filter(Negate(is.null), lmm_estimates_list)
merf_estimates_list   <- Filter(Negate(is.null), merf_estimates_list)
sample_counts_list    <- Filter(Negate(is.null), sample_counts_list)

# ==== Step 4: Summarize Results ====

# 4.1 Direct Summary
all_direct <- rbindlist(direct_estimates_list)
summary_direct <- all_direct[, .(
  direct_Mean  = mean(Mean, na.rm = TRUE),
  direct_SE    = sd(Mean, na.rm = TRUE),
  direct_MSE   = mean((Mean - True_Mean)^2, na.rm = TRUE),
  direct_RMSE  = sqrt(mean((Mean - True_Mean)^2, na.rm = TRUE)),
  direct_RB    = mean((Mean - True_Mean) / True_Mean, na.rm = TRUE),
  direct_RRMSE = sqrt(mean((Mean - True_Mean)^2, na.rm = TRUE)) / mean(True_Mean, na.rm = TRUE),
  direct_CV    = 100 * sd(Mean, na.rm = TRUE) / mean(Mean, na.rm = TRUE)
), by = Domain]

# 4.2 LMM Summary
all_lmm <- rbindlist(lmm_estimates_list)
summary_lmm <- all_lmm[, .(
  LMM_Mean  = mean(MA_Est, na.rm = TRUE),
  LMM_SE    = sd(MA_Est, na.rm = TRUE),
  LMM_MSE   = mean((MA_Est - True_Mean)^2, na.rm = TRUE),
  LMM_RMSE  = sqrt(mean((MA_Est - True_Mean)^2, na.rm = TRUE)),
  LMM_RB    = mean((MA_Est - True_Mean) / True_Mean, na.rm = TRUE),
  LMM_RRMSE = sqrt(mean((MA_Est - True_Mean)^2, na.rm = TRUE)) / mean(True_Mean, na.rm = TRUE),
  LMM_CV    = 100 * sd(MA_Est, na.rm = TRUE) / mean(MA_Est, na.rm = TRUE)
), by = Domain]

# 4.3 MERF Summary (aggregating across MERF1, MERF2, MERF3)
all_merf <- rbindlist(merf_estimates_list)
summary_merf <- all_merf[, .(
  MERF_Mean  = mean(MA_Est, na.rm = TRUE),
  MERF_SE    = sd(MA_Est, na.rm = TRUE),
  MERF_MSE   = mean((MA_Est - True_Mean)^2, na.rm = TRUE),
  MERF_RMSE  = sqrt(mean((MA_Est - True_Mean)^2, na.rm = TRUE)),
  MERF_RB    = mean((MA_Est - True_Mean) / True_Mean, na.rm = TRUE),
  MERF_RRMSE = sqrt(mean((MA_Est - True_Mean)^2, na.rm = TRUE)) / mean(True_Mean, na.rm = TRUE),
  MERF_CV    = 100 * sd(MA_Est, na.rm = TRUE) / mean(MA_Est, na.rm = TRUE)
), by = .(Domain, Model)]

# ==== Step 5: (Optional) Reshape MERF Summary Metrics for easier comparison ====
reshape_metric <- function(metric, suffix) {
  models_present <- unique(summary_merf$Model)
  df <- dcast(summary_merf, Domain ~ Model, value.var = metric, drop = FALSE)
  new_names <- paste0(models_present, "_", suffix)
  setnames(df, old = models_present, new = new_names)
  return(df)
}

merf_mean_wide   <- reshape_metric("MERF_Mean", "Mean")
merf_se_wide     <- reshape_metric("MERF_SE", "SE")
merf_mse_wide    <- reshape_metric("MERF_MSE", "MSE")
merf_rmse_wide   <- reshape_metric("MERF_RMSE", "RMSE")
merf_rb_wide     <- reshape_metric("MERF_RB", "RB")
merf_rrmse_wide  <- reshape_metric("MERF_RRMSE", "RRMSE")
merf_cv_wide     <- reshape_metric("MERF_CV", "CV")

merged_merf_metrics <- Reduce(function(x, y) merge(x, y, by = "Domain"), 
                              list(merf_mean_wide, merf_se_wide, merf_mse_wide, 
                                   merf_rmse_wide, merf_rb_wide, merf_rrmse_wide, merf_cv_wide))

# ==== Step 6: Merge All Summaries into Final Table ====
final_table <- merge(merged_merf_metrics, summary_direct, by = "Domain", all.x = TRUE)
final_table <- merge(final_table, summary_lmm, by = "Domain", all.x = TRUE)

# ==== Step 7: Average Sample Counts per Domain ====
sample_counts_all <- rbindlist(sample_counts_list)
sample_counts_avg <- sample_counts_all[, .(Avg_Sample_Count = mean(Sample_Count)), by = Domain]
final_table <- merge(final_table, sample_counts_avg, by = "Domain", all.x = TRUE)

# (Optional) View the final table in RStudio
View(final_table)

# ==== Step 8: Optionally Save to CSV ==== 
# fwrite(final_table, "C:/path/to/Y_Inter_Par_summary.csv")


############ Summaries (for RMSE)
MERF1_summary = summary(final_table$MERF1_RMSE)
MERF2_summary = summary(final_table$MERF2_RMSE)
MERF3_summary = summary(final_table$MERF3_RMSE)
LMM_summary = summary(final_table$LMM_RMSE)
direct_summary = summary(final_table$direct_RMSE)

print("MERF1 RMSE Summary:")
print(MERF1_summary)
print("MERF2 RMSE Summary:")
print(MERF2_summary)
print("MERF3 RMSE Summary:")
print(MERF3_summary)
print("LMM RMSE Summary:")
print(LMM_summary)
print("Direct RMSE Summary:")
print(direct_summary)

# ==== Plot: Boxplot for RMSE Comparison across all models ====
boxplot(final_table$MERF1_RMSE,
        final_table$MERF2_RMSE,
        final_table$MERF3_RMSE,
        final_table$LMM_RMSE,
        final_table$direct_RMSE,
        names = c("MERF1", "MERF2", "MERF3", "LMM", "Direct"),
        col = c("tomato", "goldenrod", "seagreen", "skyblue", "orchid"),
        main = "RMSE Comparison of All Models",
        ylab = "RMSE",
        outline = TRUE)

fwrite(final_table, "Accurate_Interaction_PAR_results_Low.csv")
