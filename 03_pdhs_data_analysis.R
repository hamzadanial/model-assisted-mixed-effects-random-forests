# Load necessary libraries
library(emdi)
library(SAEforest)
library(data.table)
library(dplyr)
library(LongituRF)
library(lme4)

# ==== Modified Generic predict Function for MERFranger ====
# This function now requires a 'newdata' argument so that non-quantile predictions
# are provided with a data argument.
predict.MERFranger <- function(object, newdata, ...) {
  # Predict using the fixed effect part (random forest) and the random effects model
  retval <- predict(object$Forest, data = newdata, ...)$predictions +
    predict(object$EffectModel, newdata = newdata, allow.new.levels = TRUE, ...)
  return(retval)
}


# Load dataset (Manually Select File)
data = fread(file.choose()) 

# Remove unwanted columns and handle missing values
df = data[, !(names(data) %in% c("V105", "M49Z")), with = FALSE]
ddf = na.omit(df)
setDT(ddf)

ddf[, weight := V005 / 1000000]  
ddf <- ddf[!is.na(weight) & weight > 0]  # filter afterwards

ddf$weight

# Define covariates
ddf[, `:=` (
  AgeChild = B19,
  SexChild = ifelse(B4 == 1, 1, 0),
  Size = M18,
  Breastfeed = ifelse(M34 %in% c(0, 100, 101, 199, 201), 1, 
                      ifelse(M34 %in% c(102:298, 199, 201:299), 2, 0)),
  Diarrhea = ifelse(H11 == 0, 0, 1),
  EduCat = ifelse(V106 == 0, 0, 1),
  Work = ifelse(V714 == 1, 1, 0),
  BMI = V445 / 100,
  Antenatal = M14,
  Resid = ifelse(V102 == 1, 1, 0),
  Wealth = ifelse(V190 <= 3, 0, 1),
  Water = ifelse(V113 %in% c(11, 12, 13, 14, 21, 31, 41, 42, 43), 1, 0),
  Toilet = ifelse(V116 %in% c(11:19, 21:24, 31:36, 41:49), 1, 0),
  child_5 = V137,
  Delivery = ifelse(M15 %in% 11:19, 0, 1),
  TV = ifelse(V120 == 1, 1, 0),
  Radio = ifelse(V121 == 1, 1, 0),
  Internet = ifelse(V171A %in% c(1, 2), 1, 0),
  Media = ifelse(V120 == 1 | V121 == 1 | V171A %in% c(1, 2), 1, 0),
  Dist = as.integer(SDIST),
  Stunting = ifelse(HW70 < -200, 1, 0)
)]

ddf$HW70

# Create final dataset
my_data = ddf[, .(Stunting, AgeChild, SexChild, Size, Breastfeed, Diarrhea,
                  EduCat, Work, BMI, Antenatal, Resid, Wealth, Water, Toilet,
                  child_5, Delivery, Media, Dist, weight)]
cor(my_data)
my_data[, Dist := as.integer(Dist)]

# Define response and predictors
y_col = "Stunting"
x_cols = c("AgeChild", "SexChild", "Size", "Breastfeed", "Diarrhea", 
           "EduCat", "Work", "BMI", "Resid", "Wealth", "Water", "Toilet", 
           "Delivery", "Media")

# Number of simulations
B = 5
s_data_valid$Stunting
dim(my_data)

# True means from full data
true_means = my_data[, .(True_Mean = mean(Stunting, na.rm = TRUE), Domain_Size = .N), by = Dist]
pop_size  <- my_data[, .(N_pop = .N), by = Dist]

# Define MERF configurations (for MERF1, MERF2, and MERF3)
merf_cfgs <- list(
  list(Model = "MERF1", num.trees = 500,  MaxIterations = 100, ErrorTolerance = 0.001, 
       mtry = 10, min.node.size = 10, sample.fraction = 0.7, B_adj = 50),
  list(Model = "MERF2", num.trees = 750,  MaxIterations = 150, ErrorTolerance = 0.0008,
       mtry = 12, min.node.size = 5,  sample.fraction = 0.8, B_adj = 100),
  list(Model = "MERF3", num.trees = 1000, MaxIterations = 200, ErrorTolerance = 0.0005, 
       mtry = 14, min.node.size = 3,  sample.fraction = 0.9, B_adj = 150)
)


# ==== Prepare Lists to Store Results ====
direct_estimates_list <- list()
lmm_estimates_list    <- list()
merf_estimates_list   <- list()
sample_counts_list    <- list()

# Simulation loop
for (b in 1:B) {
  set.seed(123 + b)
  sample_index = sample(1:nrow(my_data), size = 0.55 * nrow(my_data), replace = TRUE)
  s_data = my_data[sample_index]
  s_data_valid = s_data[!is.na(s_data$weight) & weight > 0]


  sc = s_data_valid[, .(Sample_Count = .N), by = Dist]
  sc[, Simulation := b]
  sample_counts_list[[b]] = sc

  if (nrow(s_data_valid) > 0) {
    direct = direct(y = y_col, smp_data = s_data_valid, smp_domains = "Dist", weights = "weight", var = TRUE)
    direct_means = estimators(direct, indicator = "Mean", MSE = TRUE, CV = TRUE)
    direct_means_df = as.data.table(direct_means)
    setnames(direct_means_df, "Domain", "Dist")
    direct_means_df[, Dist := as.integer(as.character(Dist))]
    direct_means_df[, Simulation := b]
    direct_estimates_list[[b]] = direct_means_df[, .(Dist, Mean, MSE, CV, Simulation)]
  }

  # GLMM
  ## --- 3.2 LMM Estimation ---
  pop_lmm <- copy(my_data)
  sampled_data_lmm <- copy(s_data_valid)
  fml <- as.formula(paste(y_col, "~", paste(x_cols, collapse = "+"), "+ (1|Dist)"))
  lmm_fit <- lmer(fml, data = sampled_data_lmm)
  
  pop_lmm[, pred_lmm := predict(lmm_fit, newdata = pop_lmm, allow.new.levels = TRUE)]
  sampled_data_lmm[, pred_lmm := predict(lmm_fit, newdata = sampled_data_lmm)]
  
  pop_sum_lmm <- pop_lmm[, .(PopSum = sum(pred_lmm)), by = Dist]
  adj_lmm <- sampled_data_lmm[, .(Sample_Adjustment = sum(weight * (get(y_col) - pred_lmm))), by = Dist]
  
  lmm_df <- merge(pop_sum_lmm, adj_lmm, by = "Dist")
  lmm_df <- merge(lmm_df, pop_size, by = "Dist")
  lmm_df <- merge(lmm_df, true_means, by = "Dist")
  lmm_df[, MA_Est := (PopSum + Sample_Adjustment) / N_pop]
  lmm_df[, `:=`(Model = "LMM", Simulation = b)]
  lmm_estimates_list[[b]] <- lmm_df[, .(Dist, MA_Est, Model, Simulation, True_Mean)]
  
########## MERF Estimation using MERFranger
for (cfg in merf_cfgs){
	merf_result = tryCatch({
		pop_merf = copy(my_data)
		sampled_data_merf = copy(s_data_valid)

		# Fit the MERFranger model
		mrf = MERFranger(
			Y = sampled_data_merf[[y_col]],
			X = sampled_data_merf[, ..x_cols],
			random = "(1|Dist)",
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

		# Obtain prediction on pop and sample level
		pop_merf[, pred_merf := predict(mrf, newdata = pop_merf)]
		sampled_data_merf[, pred_merf := predict(mrf, newdata = sampled_data_merf)]

		# Aggregate predictions at district level
		pop_sum_merf = pop_merf[, .(PopSum = sum(pred_merf)), by = Dist]
		adj_merf = sampled_data_merf[, .(Sample_Adjustment = sum(weight * (get(y_col) - pred_merf))), by = Dist]

		# Combine results into a single data frame
		est_df = merge(pop_sum_merf, adj_merf, by = "Dist")
		est_df = merge(est_df, pop_size, by = "Dist")
		est_df = merge(est_df, true_means, by = "Dist")
		est_df[, MA_Est := (PopSum + Sample_Adjustment) / N_pop]
		est_df[, `:=`(Model = cfg$Model, Simulation = b)]  # Make sure the Model column is set

		est_df[, .(Dist, MA_Est, Model,  True_Mean)]  # Ensure Model is included in output

	}, error = function(e){
		message("MERF model failed for ", cfg$Model, " in simulation ", b, ": ", e$message)
		data.table(Domain = unique(pop$Dist),
				   MA_Est = NA_real_,
				   Model = cfg$Model,
				   Simulation = b,
				   True_mean = NA_real_)
	})

	# Append the results
	merf_key = paste0(cfg$Model, "_", b)
	merf_estimates_list[[merf_key]] = merf_result
}


}

# 4.1 Direct Summary
all_direct <- rbindlist(direct_estimates_list)
# Merge the true_means into all_direct based on Dist (ensure True_Mean is numeric)
all_direct <- merge(all_direct, true_means, by = "Dist", all.x = TRUE)

# Ensure True_Mean is numeric (this column should now exist after merging)
all_direct[, True_Mean := as.numeric(True_Mean)]

# Now you can calculate the summary statistics
summary_direct <- all_direct[, .(
  direct_Mean  = mean(Mean, na.rm = TRUE),
  direct_SE    = sd(Mean, na.rm = TRUE),
  direct_MSE   = mean((Mean - True_Mean)^2, na.rm = TRUE),
  direct_RMSE  = sqrt(mean((Mean - True_Mean)^2, na.rm = TRUE)),
  direct_RB    = mean((Mean - True_Mean) / True_Mean, na.rm = TRUE),
  direct_RRMSE = sqrt(mean((Mean - True_Mean)^2, na.rm = TRUE)) / mean(True_Mean, na.rm = TRUE),
  direct_CV    = 100 * sd(Mean, na.rm = TRUE) / mean(Mean, na.rm = TRUE)
), by = Dist]



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
), by = Dist]


# Combine all MERF results
all_merf <- rbindlist(merf_estimates_list, fill = TRUE)


# Create summary for MERF
summary_merf <- all_merf[, .(
  MERF_Mean  = mean(MA_Est, na.rm = TRUE),
  MERF_SE    = sd(MA_Est, na.rm = TRUE),
  MERF_MSE   = mean((MA_Est - True_Mean)^2, na.rm = TRUE),
  MERF_RMSE  = sqrt(mean((MA_Est - True_Mean)^2, na.rm = TRUE)),
  MERF_RB    = mean((MA_Est - True_Mean) / True_Mean, na.rm = TRUE),
  MERF_RRMSE = sqrt(mean((MA_Est - True_Mean)^2, na.rm = TRUE)) / mean(True_Mean, na.rm = TRUE),
  MERF_CV    = 100 * sd(MA_Est, na.rm = TRUE) / mean(MA_Est, na.rm = TRUE)
), by = .(Dist, Model)]



# ==== Step 5: (Optional) Reshape MERF Summary Metrics for easier comparison ====
reshape_metric <- function(metric, suffix) {
  models_present <- unique(summary_merf$Model)
  df <- dcast(summary_merf, Dist ~ Model, value.var = metric, drop = FALSE)
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

merged_merf_metrics <- Reduce(function(x, y) merge(x, y, by = "Dist"), 
                              list(merf_mean_wide, merf_se_wide, merf_mse_wide, 
                                   merf_rmse_wide, merf_rb_wide, merf_rrmse_wide, merf_cv_wide))

# ==== Step 6: Merge All Summaries into Final Table ====
final_table <- merge(merged_merf_metrics, summary_direct, by = "Dist", all.x = TRUE)
final_table <- merge(final_table, summary_lmm, by = "Dist", all.x = TRUE)

# ==== Step 7: Average Sample Counts per Domain ====
sample_counts_all <- rbindlist(sample_counts_list)
sample_counts_avg <- sample_counts_all[, .(Avg_Sample_Count = mean(Sample_Count)), by = Dist]
final_table <- merge(final_table, sample_counts_avg, by = "Dist", all.x = TRUE)

# (Optional) View the final table in RStudio
View(final_table)

# ==== Step 8: Optionally Save to CSV ==== 
# fwrite(final_table, "C:/path/to/Y_Inter_Par_summary.csv")


############ Summaries (for RMSE)
MERF1_summary = summary(final_table$MERF1_Mean)
MERF2_summary = summary(final_table$MERF2_Mean)
MERF3_summary = summary(final_table$MERF3_Mean)
LMM_summary = summary(final_table$LMM_Mean)
direct_summary = summary(final_table$direct_Mean)

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
summary(true_means$True_Mean)
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

fwrite(final_table, "final_PAKISTAN_ZINDABAD_results.csv")








