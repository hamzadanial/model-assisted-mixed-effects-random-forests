## Bootstrap MSE
## MERF3 (unweighted) vs MA-MERF3 (weighted)
## ONLY 2 models: MERF and MA-MERF3
## both have same hyperparameters
## 18 CPU cores for simulations (1 thread per worker)

# -------------------------------
# THREAD CONTROL (important on server)
# -------------------------------
Sys.setenv(
  OMP_NUM_THREADS="1",
  MKL_NUM_THREADS="1",
  OPENBLAS_NUM_THREADS="1",
  VECLIB_MAXIMUM_THREADS="1",
  RCPP_PARALLEL_NUM_THREADS="1"
)

suppressPackageStartupMessages({
  library(data.table)
  library(SAEforest)
  library(lme4)
  library(ranger)
  library(parallel)
  library(pbapply)
})

# keep data.table single-threaded (per worker)
data.table::setDTthreads(1)

# -------------------------------
# SAFE SEED (avoids "seed not valid integer")
# -------------------------------
safe_seed <- function(x) {
  sd <- suppressWarnings(as.integer(x))
  if (is.na(sd)) sd <- 123L
  sd <- as.integer(sd %% 2147483646L)
  if (sd <= 0L) sd <- 1L
  sd
}

# -------------------------------
# FIX predict() for MERFranger (force 1 thread)
# -------------------------------
predict.MERFranger <- function(object, newdata, ...) {
  pred_fix <- predict(object$Forest, data = newdata, num.threads = 1, ...)$predictions
  pred_re  <- predict(object$EffectModel, newdata = newdata, allow.new.levels = TRUE)
  as.numeric(pred_fix + pred_re)
}

# -------------------------------
# Fit MERF + compute domain mean estimator:
# Mean = (PopSum + Adj) / N_pop
# Adj  = sum(y - pred) OR sum(w*(y - pred))  (domain-wise)
# -------------------------------
fit_merf_mean <- function(pop, smp, y_col, Xvars, cfg, use_weights) {
  
  fit <- suppressWarnings(MERFranger(
    Y = smp[[y_col]],
    X = smp[, ..Xvars],
    random = "(1|Domain)",
    data = smp,
    num.trees       = cfg$num.trees,
    MaxIterations   = cfg$MaxIterations,
    ErrorTolerance  = cfg$ErrorTolerance,
    B_adj           = cfg$B_adj,
    mtry            = cfg$mtry,
    min.node.size   = cfg$min.node.size,
    sample.fraction = cfg$sample.fraction,
    num.threads     = 1,
    na.rm = TRUE
  ))
  
  # predictions
  pop_pred <- predict(fit, newdata = pop)
  smp_pred <- predict(fit, newdata = smp)
  
  domains <- sort(unique(pop$Domain))
  dom_chr <- as.character(domains)
  
  # PopSum and N_pop
  PopSum <- tapply(pop_pred, pop$Domain, sum)
  N_pop  <- as.numeric(table(pop$Domain))
  
  # domain adjustment
  if (use_weights) {
    Adj <- tapply(seq_len(nrow(smp)), smp$Domain, function(ii) {
      sum(smp$weight[ii] * (smp[[y_col]][ii] - smp_pred[ii]), na.rm = TRUE)
    })
  } else {
    Adj <- tapply(seq_len(nrow(smp)), smp$Domain, function(ii) {
      sum((smp[[y_col]][ii] - smp_pred[ii]), na.rm = TRUE)
    })
  }
  
  # align vectors
  PopSum_a <- as.numeric(PopSum[dom_chr]); PopSum_a[is.na(PopSum_a)] <- 0
  Adj_a    <- as.numeric(Adj[dom_chr]);    Adj_a[is.na(Adj_a)] <- 0
  
  Mean_Est <- (PopSum_a + Adj_a) / N_pop
  
  est <- data.table(
    Domain   = domains,
    Mean_Est = as.numeric(Mean_Est)
  )
  
  list(fit = fit, est = est)
}

# -------------------------------
# Bootstrap MSE with REFIT (same methodology)
# BootMSE_i = mean_b (mu_star_i - mu_hat_i)^2
# -------------------------------
boot_mse_refit <- function(pop, smp, y_col, Xvars, base_fit, cfg,
                           use_weights, Bboot, seed) {
  
  domains <- sort(unique(pop$Domain))
  dom_chr <- as.character(domains)
  D <- length(domains)
  
  # mapping each unit -> domain index (for u_b assignment)
  dom_idx <- match(pop$Domain, domains)
  
  # RF-only predictions from base forest (fixed part)
  fhat_pop_fix <- as.numeric(predict(base_fit$Forest, data = pop[, ..Xvars], num.threads = 1)$predictions)
  fhat_smp_fix <- as.numeric(predict(base_fit$Forest, data = smp[, ..Xvars], num.threads = 1)$predictions)
  
  # residuals on sample
  e_hat <- smp[[y_col]] - fhat_smp_fix
  
  # domain residual means r_hat
  if (use_weights) {
    r_hat <- tapply(seq_len(nrow(smp)), smp$Domain, function(ii) {
      sum(smp$weight[ii] * e_hat[ii]) / sum(smp$weight[ii])
    })
  } else {
    r_hat <- tapply(e_hat, smp$Domain, mean)
  }
  
  # within-domain residuals
  r1_hat <- e_hat - r_hat[match(smp$Domain, names(r_hat))]
  
  # center r1
  r1_mean <- if (use_weights) sum(smp$weight * r1_hat) / sum(smp$weight) else mean(r1_hat)
  r1_c <- r1_hat - r1_mean
  
  # center r2
  if (use_weights) {
    w_dom <- tapply(smp$weight, smp$Domain, sum)
    r2_mean <- sum(w_dom * r_hat) / sum(w_dom)
  } else {
    w_dom <- tapply(rep(1, nrow(smp)), smp$Domain, sum)
    r2_mean <- mean(r_hat)
  }
  r2_c <- r_hat - r2_mean
  
  # variance estimates
  sigma2_eps_hat <- if (use_weights) sum(smp$weight * (r1_c^2)) / sum(smp$weight) else mean(r1_c^2)
  sigma2_u_hat   <- if (use_weights) sum(w_dom * (r2_c^2)) / sum(w_dom)           else mean(r2_c^2)
  
  # OOB bias correction (kept)
  sigma2_eps_bc <- sigma2_eps_hat
  if (!is.null(base_fit$OOBresiduals)) {
    fhat_oob_fix <- as.numeric(base_fit$Forest$predictions)
    k_hat <- if (use_weights) {
      sum(smp$weight * (fhat_smp_fix - fhat_oob_fix)^2) / sum(smp$weight)
    } else {
      mean((fhat_smp_fix - fhat_oob_fix)^2)
    }
    temp <- sigma2_eps_hat - k_hat
    if (is.finite(temp) && temp > 0) sigma2_eps_bc <- temp
  }
  
  # scale r1 and r2
  r1_scaled <- if (is.finite(sigma2_eps_hat) && sigma2_eps_hat > 0) {
    r1_c / sqrt(sigma2_eps_hat) * sqrt(sigma2_eps_bc)
  } else rep(0, length(r1_c))
  
  var_r2c <- var(r2_c)
  r2_scaled <- if (is.finite(var_r2c) && var_r2c > 0 && is.finite(sigma2_u_hat) && sigma2_u_hat > 0) {
    r2_c / sqrt(var_r2c) * sqrt(sigma2_u_hat)
  } else rep(0, length(r2_c))
  
  # domain sizes
  N_i <- as.numeric(table(pop$Domain))  # aligned with domains
  n_i_tab <- table(smp$Domain)
  n_i <- as.integer(n_i_tab[dom_chr]); n_i[is.na(n_i)] <- 0L
  
  # indices by domain (for bootstrap sampling)
  pop_by_dom <- split(seq_len(nrow(pop)), pop$Domain)
  
  # incremental accumulator
  mse_acc <- rep(0, D)
  
  # bootstrap loop
  base_sd <- safe_seed(seed)
  
  for (b in seq_len(Bboot)) {
    
    set.seed(safe_seed(base_sd + b))
    
    # draw one u* per domain and assign to each unit
    r2_b <- sample(r2_scaled, D, replace = TRUE)
    u_b  <- r2_b[dom_idx]
    
    # draw eps* for each population unit from within-domain residual pool
    r1_b <- sample(r1_scaled, nrow(pop), replace = TRUE)
    
    # bootstrap population Y*
    y_star <- fhat_pop_fix + u_b + r1_b
    
    # mu_star (true bootstrap domain means)
    mu_star <- tapply(y_star, pop$Domain, mean)
    mu_star_a <- as.numeric(mu_star[dom_chr]); mu_star_a[is.na(mu_star_a)] <- 0
    
    # bootstrap sample with same n_i per domain
    samp_ids <- unlist(lapply(seq_len(D), function(k) {
      if (n_i[k] <= 0) return(integer(0))
      idx <- pop_by_dom[[dom_chr[k]]]
      sample(idx, n_i[k])
    }), use.names = FALSE)
    
    # minimal bootstrap sample (only needed cols)
    smp_b <- pop[samp_ids, c("Domain", "weight", Xvars), with = FALSE]
    smp_b[, Yb := as.numeric(y_star[samp_ids])]
    
    # REFIT MERF (unchanged methodology)
    merf_b <- suppressWarnings(MERFranger(
      Y = smp_b$Yb,
      X = smp_b[, ..Xvars],
      random = "(1|Domain)",
      data = smp_b,
      num.trees       = cfg$num.trees,
      MaxIterations   = cfg$MaxIterations,
      ErrorTolerance  = cfg$ErrorTolerance,
      B_adj           = cfg$B_adj,
      mtry            = cfg$mtry,
      min.node.size   = cfg$min.node.size,
      sample.fraction = cfg$sample.fraction,
      num.threads     = 1,
      na.rm = TRUE
    ))
    
    pred_pop_b <- predict(merf_b, newdata = pop)
    pred_smp_b <- predict(merf_b, newdata = smp_b)
    
    pop_sum_b <- tapply(pred_pop_b, pop$Domain, sum)
    
    if (use_weights) {
      adj_b <- tapply(seq_len(nrow(smp_b)), smp_b$Domain, function(ii) {
        sum(smp_b$weight[ii] * (smp_b$Yb[ii] - pred_smp_b[ii]), na.rm = TRUE)
      })
    } else {
      adj_b <- tapply(seq_len(nrow(smp_b)), smp_b$Domain, function(ii) {
        sum(smp_b$Yb[ii] - pred_smp_b[ii], na.rm = TRUE)
      })
    }
    
    ps <- as.numeric(pop_sum_b[dom_chr]); ps[is.na(ps)] <- 0
    ad <- as.numeric(adj_b[dom_chr]);     ad[is.na(ad)] <- 0
    
    mu_hat <- (ps + ad) / N_i
    
    mse_acc <- mse_acc + (mu_star_a - mu_hat)^2
  }
  
  BootMSE <- mse_acc / Bboot
  data.table(Domain = domains, BootMSE = as.numeric(BootMSE))
}

# -------------------------------
#simulation
# -------------------------------
run_one_sim <- function(b, pop, n_samp, y_col, Xvars, cfg, Bboot, z) {
  
  cat(sprintf("[%s] START sim %d (pid=%d)\n",
              format(Sys.time(), "%H:%M:%S"), b, Sys.getpid()))
  
  set.seed(safe_seed(100000L + b))
  
  smp <- pop[sample(.N, n_samp)]
  if (!("weight" %in% names(smp))) stop("Missing 'weight' column (required for MA-MERF3).")
  
  out_list_est <- list()
  out_list_mse <- list()
  
  # ---- MERF3 (unweighted)
  fit1 <- fit_merf_mean(pop, smp, y_col, Xvars, cfg, use_weights = FALSE)
  est1 <- copy(fit1$est)[, `:=`(Simulation = b, Model = "MERF3")]
  
  mse1 <- boot_mse_refit(
    pop, smp, y_col, Xvars,
    base_fit = fit1$fit, cfg = cfg,
    use_weights = FALSE,
    Bboot = Bboot,
    seed = safe_seed(500000L + 1000L*b + 3L)
  )
  mse1[, `:=`(Simulation = b, Model = "MERF3")]
  
  cat(sprintf("[%s] DONE  sim %d: MERF3 finished\n", format(Sys.time(), "%H:%M:%S"), b))
  
  # ---- MA-MERF3 (weighted)
  fit2 <- fit_merf_mean(pop, smp, y_col, Xvars, cfg, use_weights = TRUE)
  est2 <- copy(fit2$est)[, `:=`(Simulation = b, Model = "MA-MERF3")]
  
  mse2 <- boot_mse_refit(
    pop, smp, y_col, Xvars,
    base_fit = fit2$fit, cfg = cfg,
    use_weights = TRUE,
    Bboot = Bboot,
    seed = safe_seed(900000L + 1000L*b + 7L)
  )
  mse2[, `:=`(Simulation = b, Model = "MA-MERF3")]
  
  cat(sprintf("[%s] DONE  sim %d: MA-MERF3 finished\n", format(Sys.time(), "%H:%M:%S"), b))
  cat(sprintf("[%s] END   sim %d\n", format(Sys.time(), "%H:%M:%S"), b))
  
  list(
    est = rbindlist(list(est1, est2), fill = TRUE),
    mse = rbindlist(list(mse1, mse2), fill = TRUE)
  )
}

# -------------------------------
# LOAD population data
# -------------------------------
pop <- fread(file.choose())
pop[, id := .I]

# Domain as factor for lme4 grouping
pop[, Domain := factor(Domain)]

# checks
stopifnot("Domain" %in% names(pop))
stopifnot("weight" %in% names(pop))

# -------------------------------
# SETTINGS (EDIT THESE)
# -------------------------------
y_col <- "Y_Normal"
Xvars <- c("X1_Normal", "X2_Normal", paste0("V_", 1:10))

B      <- 500
n_samp <- 1200
Bboot  <- 200

alpha <- 0.05
z <- qnorm(1 - alpha/2)

# SAME hyperparameters for MERF3 and MA-MERF3
cfg_shared <- list(
  num.trees       = 1000,
  MaxIterations   = 200,
  ErrorTolerance  = 0.0005,
  mtry            = 12,
  min.node.size   = 3,
  sample.fraction = 0.9,
  B_adj           = 150
)

# True means (population truth)
true_dt <- pop[, .(True_Mean = mean(get(y_col), na.rm = TRUE)), by = Domain]

# -------------------------------
# PARALLEL SETUP (18 CORES)
# -------------------------------
max_avail <- parallel::detectCores()
n_cores_sim <- min(18, max_avail)

cat("\nDetected cores:", max_avail, "| Using cores:", n_cores_sim, "\n")
cat("R PID:", Sys.getpid(), "\n\n")

# reproducible parallel RNG
RNGkind("L'Ecuyer-CMRG")
set.seed(20260213)

# outfile="" shows worker cat() output in console
cl <- makeCluster(n_cores_sim, outfile = "")

clusterEvalQ(cl, {
  library(data.table)
  library(SAEforest)
  library(lme4)
  library(ranger)
  data.table::setDTthreads(1)
  Sys.setenv(
    OMP_NUM_THREADS="1",
    MKL_NUM_THREADS="1",
    OPENBLAS_NUM_THREADS="1",
    VECLIB_MAXIMUM_THREADS="1",
    RCPP_PARALLEL_NUM_THREADS="1"
  )
})

clusterSetRNGStream(cl, iseed = 20260213)

clusterExport(
  cl,
  varlist = c(
    "pop","n_samp","y_col","Xvars","cfg_shared","Bboot","z","true_dt",
    "safe_seed","predict.MERFranger","fit_merf_mean","boot_mse_refit","run_one_sim"
  ),
  envir = environment()
)

# progress bar (moves when a simulation finishes)
pbapply::pboptions(type = "txt", char = "=")
cat("\nRunning simulations with progress bar...\n")

sim_out <- pbapply::pblapply(
  1:B,
  function(b) run_one_sim(b, pop, n_samp, y_col, Xvars, cfg_shared, Bboot, z),
  cl = cl
)

stopCluster(cl)

# -------------------------------
# COMBINE OUTPUTS
# -------------------------------
est_all <- rbindlist(lapply(sim_out, `[[`, "est"), fill = TRUE)
mse_all <- rbindlist(lapply(sim_out, `[[`, "mse"), fill = TRUE)

# attach truth
est_all <- merge(est_all, true_dt, by = "Domain", all.x = TRUE)

# -------------------------------
# EMPIRICAL RESULTS per domain/model
# -------------------------------
emp_domain <- est_all[, .(
  Mean_Est   = mean(Mean_Est, na.rm = TRUE),
  SE_Est     = sd(Mean_Est, na.rm = TRUE),
  MSE_emp    = mean((Mean_Est - True_Mean)^2, na.rm = TRUE),
  RMSE_emp   = sqrt(mean((Mean_Est - True_Mean)^2, na.rm = TRUE)),
  PARB_pct   = 100 * mean(abs((Mean_Est - True_Mean) / True_Mean), na.rm = TRUE),
  PRRMSE_pct = 100 * (sqrt(mean((Mean_Est - True_Mean)^2, na.rm = TRUE)) / abs(unique(True_Mean)))
), by = .(Domain, Model)]

# mean Bootstrap MSE per domain/model
mse_domain_mean <- mse_all[, .(MSE_boot = mean(BootMSE, na.rm = TRUE)), by = .(Domain, Model)]
mse_domain_mean[, RMSE_boot := sqrt(MSE_boot)]

# -------------------------------
# COVERAGE using BootMSE
# -------------------------------
tmp_cov <- merge(
  est_all[, .(Domain, Model, Simulation, Mean_Est, True_Mean)],
  mse_all[, .(Domain, Model, Simulation, BootMSE)],
  by = c("Domain", "Model", "Simulation"),
  all.x = TRUE
)

tmp_cov[, RB := (Mean_Est - True_Mean) / True_Mean]
tmp_cov[, RB_hw := z * (sqrt(BootMSE) / abs(True_Mean))]
tmp_cov[, Cover_RB := as.integer(RB >= -RB_hw & RB <= RB_hw)]
cov_domain_rb <- tmp_cov[, .(Coverage_RB = mean(Cover_RB, na.rm = TRUE)), by = .(Domain, Model)]

tmp_cov[, lower := Mean_Est - z * sqrt(BootMSE)]
tmp_cov[, upper := Mean_Est + z * sqrt(BootMSE)]
tmp_cov[, Cover := as.integer(True_Mean >= lower & True_Mean <= upper)]
cov_domain_ci <- tmp_cov[, .(Coverage_prob = mean(Cover, na.rm = TRUE)), by = .(Domain, Model)]

# RMSE boot diagnostics
rmse_emp_dt <- emp_domain[, .(Domain, Model, RMSE_emp)]
rmse_bias_all <- rbindlist(lapply(unique(mse_all$Model), function(md) {
  
  mse_md  <- mse_all[Model == md, .(Domain, Simulation, BootMSE)]
  rmse_md <- rmse_emp_dt[Model == md, .(Domain, RMSE_emp)]
  
  tmp <- merge(mse_md, rmse_md, by = "Domain", all.x = TRUE)
  tmp[, s := sqrt(BootMSE)]
  
  out <- tmp[, .(
    RB_RMSE_boot    = (mean(s, na.rm = TRUE) - unique(RMSE_emp)) / unique(RMSE_emp),
    RRMSE_RMSE_boot = sqrt(mean((s - unique(RMSE_emp))^2, na.rm = TRUE)) / unique(RMSE_emp)
  ), by = Domain]
  
  out[, Model := md]
  out
}), fill = TRUE)

# final results table
final_domain <- merge(emp_domain, cov_domain_rb,  by = c("Domain", "Model"), all.x = TRUE)
final_domain <- merge(final_domain, cov_domain_ci, by = c("Domain", "Model"), all.x = TRUE)
final_domain <- merge(final_domain, mse_domain_mean, by = c("Domain", "Model"), all.x = TRUE)
final_domain <- merge(final_domain, rmse_bias_all,   by = c("Domain", "Model"), all.x = TRUE)

final_domain[, Coverage_main := Coverage_RB]

setcolorder(final_domain, c(
  "Domain","Model",
  "Mean_Est","SE_Est",
  "MSE_emp","RMSE_emp",
  "PARB_pct","PRRMSE_pct",
  "MSE_boot","RMSE_boot",
  "RB_RMSE_boot","RRMSE_RMSE_boot",
  "Coverage_main","Coverage_RB","Coverage_prob"
))

print(final_domain)

# summary table (Min, Q1, Median, Mean, Q3, Max)
get_stats <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) == 0) return(c(NA, NA, NA, NA, NA, NA))
  q <- quantile(v, probs = c(0, .25, .5, .75, 1), na.rm = TRUE, names = FALSE)
  c(q[1], q[2], q[3], mean(v, na.rm = TRUE), q[4], q[5])
}

models <- c("MERF3","MA-MERF3")
summary_table <- rbindlist(lapply(models, function(m) {
  
  x <- final_domain[Model == m]
  
  a1 <- get_stats(x$Mean_Est)
  a2 <- get_stats(x$PARB_pct)
  a3 <- get_stats(x$PRRMSE_pct)
  a4 <- get_stats(x$Coverage_main)
  a5 <- get_stats(x$MSE_emp)
  a6 <- get_stats(x$MSE_boot)
  a7 <- get_stats(x$RMSE_boot)
  a8 <- get_stats(x$RB_RMSE_boot)
  a9 <- get_stats(x$RRMSE_RMSE_boot)
  
  rbindlist(list(
    data.table(Model=m, Measure="Mean_Est",        Min=a1[1], Q1=a1[2], Median=a1[3], Mean=a1[4], Q3=a1[5], Max=a1[6]),
    data.table(Model=m, Measure="PARB_pct",        Min=a2[1], Q1=a2[2], Median=a2[3], Mean=a2[4], Q3=a2[5], Max=a2[6]),
    data.table(Model=m, Measure="PRRMSE_pct",      Min=a3[1], Q1=a3[2], Median=a3[3], Mean=a3[4], Q3=a3[5], Max=a3[6]),
    data.table(Model=m, Measure="Coverage_prob",   Min=a4[1], Q1=a4[2], Median=a4[3], Mean=a4[4], Q3=a4[5], Max=a4[6]),
    data.table(Model=m, Measure="MSE_emp",         Min=a5[1], Q1=a5[2], Median=a5[3], Mean=a5[4], Q3=a5[5], Max=a5[6]),
    data.table(Model=m, Measure="MSE_boot",        Min=a6[1], Q1=a6[2], Median=a6[3], Mean=a6[4], Q3=a6[5], Max=a6[6]),
    data.table(Model=m, Measure="RMSE_boot",       Min=a7[1], Q1=a7[2], Median=a7[3], Mean=a7[4], Q3=a7[5], Max=a7[6]),
    data.table(Model=m, Measure="RB_RMSE_boot",    Min=a8[1], Q1=a8[2], Median=a8[3], Mean=a8[4], Q3=a8[5], Max=a8[6]),
    data.table(Model=m, Measure="RRMSE_RMSE_boot", Min=a9[1], Q1=a9[2], Median=a9[3], Mean=a9[4], Q3=a9[5], Max=a9[6])
  ), fill = TRUE)
}), fill = TRUE)

print(summary_table)

# save output files
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
file_final   <- paste0("FinalResults_",  y_col, "_", stamp, ".csv")
file_summary <- paste0("SummaryTable_", y_col, "_", stamp, ".csv")

fwrite(final_domain,  file_final)
fwrite(summary_table, file_summary)

cat("\nSaved files:\n")
cat("1) ", file_final, "\n")
cat("2) ", file_summary, "\n")

