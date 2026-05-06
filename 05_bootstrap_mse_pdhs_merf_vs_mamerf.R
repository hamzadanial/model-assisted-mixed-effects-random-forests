
# PDHS: MERF3 and MA-MERF3 with Bootstrap MSE
# Only MERF3 and MA-MERF3


# -------------------------------
# THREAD CONTROL
# -------------------------------
Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  RCPP_PARALLEL_NUM_THREADS = "1"
)

suppressPackageStartupMessages({
  library(data.table)
  library(SAEforest)
  library(lme4)
  library(ranger)
  library(parallel)
  library(pbapply)
  library(survey)
})

data.table::setDTthreads(1)

# -------------------------------
# SAFE SEED
# -------------------------------
safe_seed <- function(x) {
  sd <- suppressWarnings(as.integer(x))
  if (is.na(sd)) sd <- 123L
  sd <- as.integer(sd %% 2147483646L)
  if (sd <= 0L) sd <- 1L
  sd
}

# -------------------------------
# FIX predict() for MERFranger
# -------------------------------
predict.MERFranger <- function(object, newdata, ...) {
  pred_fix <- predict(object$Forest, data = newdata, num.threads = 1, ...)$predictions
  pred_re  <- predict(object$EffectModel, newdata = newdata, allow.new.levels = TRUE)
  as.numeric(pred_fix + pred_re)
}

# -------------------------------
# HELPER: WEIGHTED DOMAIN MEAN
# -------------------------------
weighted_domain_mean_vec <- function(values, weights, domain, domains = levels(domain)) {
  out <- tapply(seq_along(values), domain, function(ii) {
    weighted.mean(values[ii], weights[ii], na.rm = TRUE)
  })
  out_a <- as.numeric(out[as.character(domains)])
  out_a[is.na(out_a)] <- 0
  out_a
}

# -------------------------------
# LOAD AND PREPARE PDHS DATA
# -------------------------------
data_raw <- fread(file.choose())

drop_cols <- intersect(c("V105", "M49Z"), names(data_raw))
if (length(drop_cols) > 0) {
  data_raw <- data_raw[, !(names(data_raw) %in% drop_cols), with = FALSE]
}

ddf <- na.omit(data_raw)
setDT(ddf)

ddf[, weight := V005 / 1000000]
ddf <- ddf[!is.na(weight) & weight > 0]

ddf[, `:=`(
  AgeChild   = B19,
  SexChild   = ifelse(B4 == 1, 1, 0),
  Size       = M18,
  Breastfeed = ifelse(M34 %in% c(0, 100, 101, 199, 201), 1,
                      ifelse(M34 %in% c(102:298, 201:299), 2, 0)),
  Diarrhea   = ifelse(H11 == 0, 0, 1),
  EduCat     = ifelse(V106 == 0, 0, 1),
  Work       = ifelse(V714 == 1, 1, 0),
  BMI        = V445 / 100,
  Antenatal  = M14,
  Resid      = ifelse(V102 == 1, 1, 0),
  Wealth     = ifelse(V190 <= 3, 0, 1),
  Water      = ifelse(V113 %in% c(11, 12, 13, 14, 21, 31, 41, 42, 43), 1, 0),
  Toilet     = ifelse(V116 %in% c(11:19, 21:24, 31:36, 41:49), 1, 0),
  child_5    = V137,
  Delivery   = ifelse(M15 %in% 11:19, 0, 1),
  TV         = ifelse(V120 == 1, 1, 0),
  Radio      = ifelse(V121 == 1, 1, 0),
  Internet   = ifelse(V171A %in% c(1, 2), 1, 0),
  Media      = ifelse(V120 == 1 | V121 == 1 | V171A %in% c(1, 2), 1, 0),
  Domain     = as.integer(SDIST),
  Stunting   = ifelse(HW70 < -200, 1, 0)
)]

pop <- ddf[, .(
  Stunting, AgeChild, SexChild, Size, Breastfeed, Diarrhea,
  EduCat, Work, BMI, Antenatal, Resid, Wealth, Water, Toilet,
  child_5, Delivery, Media, Domain, weight
)]

y_col <- "Stunting"
Xvars <- c("AgeChild", "SexChild", "Size", "Breastfeed", "Diarrhea",
           "EduCat", "Work", "BMI", "Resid", "Wealth", "Water",
           "Toilet", "Delivery", "Media")

need_vars <- c(y_col, Xvars, "Domain", "weight")
pop <- pop[complete.cases(pop[, ..need_vars])]

pop[, Domain := factor(Domain)]
pop[, id := .I]

cat("Total cleaned pseudo-population size =", nrow(pop), "\n")
cat("Number of domains =", length(unique(pop$Domain)), "\n")

# -------------------------------
# TRUE DOMAIN MEANS FROM FULL PDHS
# -------------------------------
des_pop <- svydesign(
  ids = ~1,
  weights = ~weight,
  data = pop
)

frm_y <- as.formula(paste0("~", y_col))

true_dt <- as.data.table(
  svyby(
    formula = frm_y,
    by = ~Domain,
    design = des_pop,
    FUN = svymean,
    na.rm = TRUE,
    vartype = NULL,
    keep.names = TRUE
  )
)

setnames(true_dt, y_col, "True_Mean")
true_dt[, N_pop := as.numeric(table(pop$Domain)[as.character(Domain)])]

# -------------------------------
# SHARED HYPERPARAMETERS
# -------------------------------
cfg_shared <- list(
  num.trees       = 1000,
  MaxIterations   = 200,
  ErrorTolerance  = 0.0005,
  mtry            = 14,
  min.node.size   = 3,
  sample.fraction = 0.9,
  B_adj           = 150
)

# -------------------------------
# FIT MERF DOMAIN MEAN
# -------------------------------
fit_merf_mean <- function(pop, smp, y_col, Xvars, cfg, use_weights = FALSE) {
  
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
  
  pop_pred <- predict(fit, newdata = pop)
  smp_pred <- predict(fit, newdata = smp)
  
  domains <- levels(pop$Domain)
  dom_chr <- as.character(domains)
  
  pop_sum <- tapply(pop_pred, pop$Domain, sum)
  N_pop   <- as.numeric(table(pop$Domain))
  
  if (use_weights) {
    adj <- tapply(seq_len(nrow(smp)), smp$Domain, function(ii) {
      sum(smp$weight[ii] * (smp[[y_col]][ii] - smp_pred[ii]), na.rm = TRUE)
    })
  } else {
    adj <- tapply(seq_len(nrow(smp)), smp$Domain, function(ii) {
      sum(smp[[y_col]][ii] - smp_pred[ii], na.rm = TRUE)
    })
  }
  
  pop_sum_a <- as.numeric(pop_sum[dom_chr]); pop_sum_a[is.na(pop_sum_a)] <- 0
  adj_a     <- as.numeric(adj[dom_chr]);     adj_a[is.na(adj_a)] <- 0
  
  mean_est <- (pop_sum_a + adj_a) / N_pop
  
  est <- data.table(
    Domain   = factor(domains, levels = domains),
    Mean_Est = as.numeric(mean_est)
  )
  
  list(fit = fit, est = est)
}

# -------------------------------
# REB-STYLE BOOTSTRAP MSE WITH REFIT
# -------------------------------
boot_mse_refit <- function(pop, smp, y_col, Xvars, base_fit, cfg,
                           use_weights = FALSE, Bboot = 20, seed = 123,
                           sim_id = NA, model_name = "") {
  
  domains <- levels(pop$Domain)
  dom_chr <- as.character(domains)
  D <- length(domains)
  
  dom_idx <- match(pop$Domain, domains)
  
  fhat_pop_fix <- as.numeric(
    predict(base_fit$Forest, data = pop[, ..Xvars], num.threads = 1)$predictions
  )
  fhat_smp_fix <- as.numeric(
    predict(base_fit$Forest, data = smp[, ..Xvars], num.threads = 1)$predictions
  )
  
  e_hat <- smp[[y_col]] - fhat_smp_fix
  
  if (use_weights) {
    r_hat <- tapply(seq_len(nrow(smp)), smp$Domain, function(ii) {
      sum(smp$weight[ii] * e_hat[ii], na.rm = TRUE) / sum(smp$weight[ii], na.rm = TRUE)
    })
  } else {
    r_hat <- tapply(e_hat, smp$Domain, mean)
  }
  
  r1_hat <- e_hat - r_hat[match(smp$Domain, names(r_hat))]
  
  r1_mean <- if (use_weights) {
    sum(smp$weight * r1_hat, na.rm = TRUE) / sum(smp$weight, na.rm = TRUE)
  } else {
    mean(r1_hat, na.rm = TRUE)
  }
  r1_c <- r1_hat - r1_mean
  
  if (use_weights) {
    w_dom <- tapply(smp$weight, smp$Domain, sum)
    r2_mean <- sum(w_dom * r_hat, na.rm = TRUE) / sum(w_dom, na.rm = TRUE)
  } else {
    w_dom <- tapply(rep(1, nrow(smp)), smp$Domain, sum)
    r2_mean <- mean(r_hat, na.rm = TRUE)
  }
  r2_c <- r_hat - r2_mean
  
  sigma2_eps_hat <- if (use_weights) {
    sum(smp$weight * (r1_c^2), na.rm = TRUE) / sum(smp$weight, na.rm = TRUE)
  } else {
    mean(r1_c^2, na.rm = TRUE)
  }
  
  sigma2_u_hat <- if (use_weights) {
    sum(w_dom * (r2_c^2), na.rm = TRUE) / sum(w_dom, na.rm = TRUE)
  } else {
    mean(r2_c^2, na.rm = TRUE)
  }
  
  sigma2_eps_bc <- sigma2_eps_hat
  if (!is.null(base_fit$OOBresiduals)) {
    fhat_oob_fix <- as.numeric(base_fit$Forest$predictions)
    
    k_hat <- if (use_weights) {
      sum(smp$weight * (fhat_smp_fix - fhat_oob_fix)^2, na.rm = TRUE) / sum(smp$weight, na.rm = TRUE)
    } else {
      mean((fhat_smp_fix - fhat_oob_fix)^2, na.rm = TRUE)
    }
    
    temp <- sigma2_eps_hat - k_hat
    if (is.finite(temp) && temp > 0) sigma2_eps_bc <- temp
  }
  
  r1_scaled <- if (is.finite(sigma2_eps_hat) && sigma2_eps_hat > 0) {
    r1_c / sqrt(sigma2_eps_hat) * sqrt(sigma2_eps_bc)
  } else {
    rep(0, length(r1_c))
  }
  
  var_r2c <- var(r2_c, na.rm = TRUE)
  r2_scaled <- if (is.finite(var_r2c) && var_r2c > 0 &&
                   is.finite(sigma2_u_hat) && sigma2_u_hat > 0) {
    r2_c / sqrt(var_r2c) * sqrt(sigma2_u_hat)
  } else {
    rep(0, length(r2_c))
  }
  
  N_i <- as.numeric(table(pop$Domain))
  
  n_i_tab <- table(smp$Domain)
  n_i <- as.integer(n_i_tab[dom_chr])
  n_i[is.na(n_i)] <- 0L
  
  pop_by_dom <- split(seq_len(nrow(pop)), pop$Domain)
  
  mse_acc <- rep(0, D)
  base_sd <- safe_seed(seed)
  
  for (bb in seq_len(Bboot)) {
    
    if (bb %% 5 == 0 || bb == 1 || bb == Bboot) {
      cat(sprintf("[%s] simulation %s | %s | bootstrap %d/%d\n",
                  format(Sys.time(), "%H:%M:%S"), sim_id, model_name, bb, Bboot))
    }
    
    set.seed(safe_seed(base_sd + bb))
    
    r2_b <- sample(r2_scaled, D, replace = TRUE)
    u_b  <- r2_b[dom_idx]
    
    r1_b <- sample(r1_scaled, nrow(pop), replace = TRUE)
    
    y_star <- fhat_pop_fix + u_b + r1_b
    y_star <- pmin(pmax(y_star, 0), 1)
    
    # WEIGHTED TRUE MEAN IN BOOTSTRAP
    mu_star_a <- weighted_domain_mean_vec(
      values  = y_star,
      weights = pop$weight,
      domain  = pop$Domain,
      domains = domains
    )
    
    samp_ids <- unlist(lapply(seq_len(D), function(k) {
      if (n_i[k] <= 0) return(integer(0))
      idx <- pop_by_dom[[dom_chr[k]]]
      sample(idx, n_i[k], replace = FALSE)
    }), use.names = FALSE)
    
    smp_b <- copy(pop[samp_ids, c("Domain", "weight", Xvars), with = FALSE])
    smp_b[, Yb := as.numeric(y_star[samp_ids])]
    
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
    
    mu_hat_b <- (ps + ad) / N_i
    mse_acc <- mse_acc + (mu_star_a - mu_hat_b)^2
  }
  
  BootMSE <- mse_acc / Bboot
  
  data.table(
    Domain  = factor(domains, levels = domains),
    BootMSE = as.numeric(BootMSE)
  )
}

# -------------------------------
#SIMULATION
# -------------------------------
run_one_sim <- function(b, pop, n_samp, y_col, Xvars, cfg, Bboot) {
  
  cat(sprintf("[%s] START simulation %d (pid=%d)\n",
              format(Sys.time(), "%H:%M:%S"), b, Sys.getpid()))
  
  set.seed(safe_seed(10000L + b))
  
  smp_ids <- sample(seq_len(nrow(pop)), size = n_samp, replace = FALSE)
  smp <- copy(pop[smp_ids])
  
  # DOMAIN-WISE SAMPLE COUNTS
  samp_counts <- data.table(
    Domain = factor(levels(pop$Domain), levels = levels(pop$Domain))
  )
  
  tmp_n <- as.data.table(table(smp$Domain))
  setnames(tmp_n, c("Domain", "Sample_Size"))
  tmp_n[, Sample_Size := as.integer(Sample_Size)]
  
  samp_counts <- merge(samp_counts, tmp_n, by = "Domain", all.x = TRUE)
  samp_counts[is.na(Sample_Size), Sample_Size := 0L]
  samp_counts[, Simulation := b]
  
  fit_merf3 <- fit_merf_mean(
    pop = pop, smp = smp, y_col = y_col, Xvars = Xvars,
    cfg = cfg, use_weights = FALSE
  )
  
  est_merf3 <- copy(fit_merf3$est)
  est_merf3[, `:=`(Simulation = b, Model = "MERF3")]
  
  mse_merf3 <- boot_mse_refit(
    pop = pop, smp = smp, y_col = y_col, Xvars = Xvars,
    base_fit = fit_merf3$fit, cfg = cfg,
    use_weights = FALSE, Bboot = Bboot,
    seed = safe_seed(500000L + 1000L * b + 3L),
    sim_id = b, model_name = "MERF3"
  )
  mse_merf3[, `:=`(Simulation = b, Model = "MERF3")]
  
  cat(sprintf("[%s] DONE simulation %d: MERF3 finished\n",
              format(Sys.time(), "%H:%M:%S"), b))
  
  fit_mamerf3 <- fit_merf_mean(
    pop = pop, smp = smp, y_col = y_col, Xvars = Xvars,
    cfg = cfg, use_weights = TRUE
  )
  
  est_mamerf3 <- copy(fit_mamerf3$est)
  est_mamerf3[, `:=`(Simulation = b, Model = "MA-MERF3")]
  
  mse_mamerf3 <- boot_mse_refit(
    pop = pop, smp = smp, y_col = y_col, Xvars = Xvars,
    base_fit = fit_mamerf3$fit, cfg = cfg,
    use_weights = TRUE, Bboot = Bboot,
    seed = safe_seed(900000L + 1000L * b + 7L),
    sim_id = b, model_name = "MA-MERF3"
  )
  mse_mamerf3[, `:=`(Simulation = b, Model = "MA-MERF3")]
  
  cat(sprintf("[%s] DONE simulation %d: MA-MERF3 finished\n",
              format(Sys.time(), "%H:%M:%S"), b))
  cat(sprintf("[%s] END simulation %d\n",
              format(Sys.time(), "%H:%M:%S"), b))
  
  list(
    est    = rbindlist(list(est_merf3, est_mamerf3), fill = TRUE),
    mse    = rbindlist(list(mse_merf3, mse_mamerf3), fill = TRUE),
    counts = samp_counts
  )
}

# -------------------------------
# SETTINGS
# SAMPLE SIZE 15% or 35%
# -------------------------------
B <- 500
n_samp <- round(0.35 * nrow(pop))
Bboot <- 200

alpha <- 0.05
z <- qnorm(1 - alpha / 2)

cat("Selected sample size per simulation =", n_samp, "\n")
cat("Number of simulations B =", B, "\n")
cat("Bootstrap replicates Bboot =", Bboot, "\n")

# -------------------------------
# PARALLEL SETUP
# -------------------------------
max_avail <- parallel::detectCores()
n_cores_sim <- min(18, B, max_avail)

cat("\nDetected cores:", max_avail, "| Using cores:", n_cores_sim, "\n")
cat("R PID:", Sys.getpid(), "\n\n")

RNGkind("L'Ecuyer-CMRG")
set.seed(20260316)

cl <- makeCluster(n_cores_sim, outfile = "")

clusterEvalQ(cl, {
  library(data.table)
  library(SAEforest)
  library(lme4)
  library(ranger)
  data.table::setDTthreads(1)
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    RCPP_PARALLEL_NUM_THREADS = "1"
  )
})

clusterSetRNGStream(cl, iseed = 20260316)

clusterExport(
  cl,
  varlist = c(
    "pop", "n_samp", "y_col", "Xvars", "cfg_shared", "Bboot", "z", "true_dt",
    "safe_seed", "predict.MERFranger", "weighted_domain_mean_vec",
    "fit_merf_mean", "boot_mse_refit", "run_one_sim"
  ),
  envir = environment()
)

pbapply::pboptions(type = "txt", char = "=")
cat("\nRunning PDHS simulations with progress bar...\n")

sim_out <- pbapply::pblapply(
  1:B,
  function(b) run_one_sim(b, pop, n_samp, y_col, Xvars, cfg_shared, Bboot),
  cl = cl
)

stopCluster(cl)

# -------------------------------
# COMBINE OUTPUTS
# -------------------------------
est_all <- rbindlist(lapply(sim_out, `[[`, "est"), fill = TRUE)
mse_all <- rbindlist(lapply(sim_out, `[[`, "mse"), fill = TRUE)
samp_counts_all <- rbindlist(lapply(sim_out, `[[`, "counts"), fill = TRUE)

est_all <- merge(est_all, true_dt, by = "Domain", all.x = TRUE)

# -------------------------------
# SCATTER MSE DATA BY SAMPLE SIZE
# -------------------------------
scatter_mse_raw <- merge(
  est_all[, .(Domain, Simulation, Model, Mean_Est, True_Mean)],
  samp_counts_all,
  by = c("Domain", "Simulation"),
  all.x = TRUE
)

scatter_mse_raw[, `:=`(
  Sq_Error  = (Mean_Est - True_Mean)^2,
  Abs_Error = abs(Mean_Est - True_Mean)
)]

scatter_mse_by_n <- scatter_mse_raw[, .(
  MSE   = mean(Sq_Error, na.rm = TRUE),
  RMSE  = sqrt(mean(Sq_Error, na.rm = TRUE)),
  Mean_Abs_Error = mean(Abs_Error, na.rm = TRUE),
  N_points = .N
), by = .(Model, Sample_Size)]

scatter_mse_by_domain <- scatter_mse_raw[, .(
  Avg_Sample_Size    = mean(Sample_Size, na.rm = TRUE),
  Median_Sample_Size = median(Sample_Size, na.rm = TRUE),
  MSE_emp            = mean(Sq_Error, na.rm = TRUE),
  RMSE_emp           = sqrt(mean(Sq_Error, na.rm = TRUE))
), by = .(Domain, Model)]

# -------------------------------
# EMPIRICAL RESULTS
# -------------------------------
emp_domain <- est_all[, .(
  Mean_Est = mean(Mean_Est, na.rm = TRUE),
  SE_Est   = sd(Mean_Est, na.rm = TRUE),
  MSE_emp  = mean((Mean_Est - True_Mean)^2, na.rm = TRUE),
  RMSE_emp = sqrt(mean((Mean_Est - True_Mean)^2, na.rm = TRUE)),
  RB       = mean((Mean_Est - True_Mean) / True_Mean, na.rm = TRUE),
  RRMSE    = sqrt(mean((Mean_Est - True_Mean)^2, na.rm = TRUE)) /
    mean(True_Mean, na.rm = TRUE),
  CV       = 100 * sd(Mean_Est, na.rm = TRUE) / mean(Mean_Est, na.rm = TRUE)
), by = .(Domain, Model)]

mse_domain_mean <- mse_all[, .(
  MSE_boot = mean(BootMSE, na.rm = TRUE)
), by = .(Domain, Model)]

mse_domain_mean[, RMSE_boot := sqrt(MSE_boot)]

# -------------------------------
# COVERAGE
# -------------------------------
tmp_cov <- merge(
  est_all[, .(Domain, Model, Simulation, Mean_Est, True_Mean)],
  mse_all[, .(Domain, Model, Simulation, BootMSE)],
  by = c("Domain", "Model", "Simulation"),
  all.x = TRUE
)

tmp_cov[, lower := Mean_Est - z * sqrt(BootMSE)]
tmp_cov[, upper := Mean_Est + z * sqrt(BootMSE)]
tmp_cov[, Cover := as.integer(True_Mean >= lower & True_Mean <= upper)]

cov_domain <- tmp_cov[, .(
  Coverage = mean(Cover, na.rm = TRUE)
), by = .(Domain, Model)]

# -------------------------------
# RMSE BOOT DIAGNOSTICS
# -------------------------------
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

# -------------------------------
# FINAL TABLE
# -------------------------------
final_domain <- merge(emp_domain, cov_domain,        by = c("Domain", "Model"), all.x = TRUE)
final_domain <- merge(final_domain, mse_domain_mean, by = c("Domain", "Model"), all.x = TRUE)
final_domain <- merge(final_domain, rmse_bias_all,   by = c("Domain", "Model"), all.x = TRUE)

setcolorder(final_domain, c(
  "Domain", "Model",
  "Mean_Est", "SE_Est",
  "MSE_emp", "RMSE_emp",
  "RB", "RRMSE", "CV",
  "MSE_boot", "RMSE_boot",
  "RB_RMSE_boot", "RRMSE_RMSE_boot",
  "Coverage"
))

print(final_domain)

# -------------------------------
# SUMMARY TABLE
# -------------------------------
get_stats <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) == 0) return(c(NA, NA, NA, NA, NA, NA))
  q <- quantile(v, probs = c(0, .25, .5, .75, 1), na.rm = TRUE, names = FALSE)
  c(q[1], q[2], q[3], mean(v, na.rm = TRUE), q[4], q[5])
}

models <- c("MERF3", "MA-MERF3")

summary_table <- rbindlist(lapply(models, function(m) {
  
  x <- final_domain[Model == m]
  
  a1  <- get_stats(x$Mean_Est)
  aSE <- get_stats(x$SE_Est)
  a2  <- get_stats(x$MSE_emp)
  a3  <- get_stats(x$RMSE_emp)
  a4  <- get_stats(x$RB)
  a5  <- get_stats(x$RRMSE)
  a6  <- get_stats(x$CV)
  a7  <- get_stats(x$MSE_boot)
  a8  <- get_stats(x$RMSE_boot)
  a9  <- get_stats(x$Coverage)
  a10 <- get_stats(x$RB_RMSE_boot)
  a11 <- get_stats(x$RRMSE_RMSE_boot)
  
  rbindlist(list(
    data.table(Model = m, Measure = "Mean_Est",        Min = a1[1],  Q1 = a1[2],  Median = a1[3],  Mean = a1[4],  Q3 = a1[5],  Max = a1[6]),
    data.table(Model = m, Measure = "SE_Est",          Min = aSE[1], Q1 = aSE[2], Median = aSE[3], Mean = aSE[4], Q3 = aSE[5], Max = aSE[6]),
    data.table(Model = m, Measure = "MSE_emp",         Min = a2[1],  Q1 = a2[2],  Median = a2[3],  Mean = a2[4],  Q3 = a2[5],  Max = a2[6]),
    data.table(Model = m, Measure = "RMSE_emp",        Min = a3[1],  Q1 = a3[2],  Median = a3[3],  Mean = a3[4],  Q3 = a3[5],  Max = a3[6]),
    data.table(Model = m, Measure = "RB",              Min = a4[1],  Q1 = a4[2],  Median = a4[3],  Mean = a4[4],  Q3 = a4[5],  Max = a4[6]),
    data.table(Model = m, Measure = "RRMSE",           Min = a5[1],  Q1 = a5[2],  Median = a5[3],  Mean = a5[4],  Q3 = a5[5],  Max = a5[6]),
    data.table(Model = m, Measure = "CV",              Min = a6[1],  Q1 = a6[2],  Median = a6[3],  Mean = a6[4],  Q3 = a6[5],  Max = a6[6]),
    data.table(Model = m, Measure = "MSE_boot",        Min = a7[1],  Q1 = a7[2],  Median = a7[3],  Mean = a7[4],  Q3 = a7[5],  Max = a7[6]),
    data.table(Model = m, Measure = "RMSE_boot",       Min = a8[1],  Q1 = a8[2],  Median = a8[3],  Mean = a8[4],  Q3 = a8[5],  Max = a8[6]),
    data.table(Model = m, Measure = "Coverage",        Min = a9[1],  Q1 = a9[2],  Median = a9[3],  Mean = a9[4],  Q3 = a9[5],  Max = a9[6]),
    data.table(Model = m, Measure = "RB_RMSE_boot",    Min = a10[1], Q1 = a10[2], Median = a10[3], Mean = a10[4], Q3 = a10[5], Max = a10[6]),
    data.table(Model = m, Measure = "RRMSE_RMSE_boot", Min = a11[1], Q1 = a11[2], Median = a11[3], Mean = a11[4], Q3 = a11[5], Max = a11[6])
  ), fill = TRUE)
}), fill = TRUE)

print(summary_table)

# -------------------------------
# SAVE RESULTS
# -------------------------------
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

file_final         <- paste0("PDHS_MERF3_MAMERF3_WeightedTruth_35pct_Final_", stamp, ".csv")
file_summary       <- paste0("PDHS_MERF3_MAMERF3_WeightedTruth_35pct_Summary_", stamp, ".csv")
file_scatter_raw   <- paste0("PDHS_MERF3_MAMERF3_WeightedTruth_35pct_ScatterRaw_", stamp, ".csv")
file_scatter_by_n  <- paste0("PDHS_MERF3_MAMERF3_WeightedTruth_35pct_ScatterBySampleSize_", stamp, ".csv")
file_scatter_dom   <- paste0("PDHS_MERF3_MAMERF3_WeightedTruth_35pct_ScatterByDomain_", stamp, ".csv")

fwrite(final_domain, file_final)
fwrite(summary_table, file_summary)
fwrite(scatter_mse_raw, file_scatter_raw)
fwrite(scatter_mse_by_n, file_scatter_by_n)
fwrite(scatter_mse_by_domain, file_scatter_dom)

cat("\nSaved files:\n")
cat("1) ", file_final, "\n")
cat("2) ", file_summary, "\n")
cat("3) ", file_scatter_raw, "\n")
cat("4) ", file_scatter_by_n, "\n")
cat("5) ", file_scatter_dom, "\n")
