library(data.table)
library(actuar)  # For rpareto

# ==== Parameters ====
set.seed(123)
D = 50
N_i = 1000
N = D * N_i
domain = rep(1:D, each = N_i)
mu_d = runif(D, -1, 1)
mu_vec = rep(mu_d, each = N_i)

# ==== Generate Predictors for Each Scenario ====
gen_X_pair = function(mu, sd1, sd2) {
  list(X1 = rnorm(N, mean = mu, sd = sd1), X2 = rnorm(N, mean = mu, sd = sd2))
}

X_Normal      = gen_X_pair(mu_vec, 3, 3)
X_Interaction = gen_X_pair(mu_vec, 4, 2)
X_NPar        = gen_X_pair(mu_vec, 3, 3)
X_IPAR        = gen_X_pair(mu_vec, 2, 2)

# ==== Random Effects & Errors ====
v_list = list(
  normal = rep(rnorm(D, 0, 500), each = N_i),
  inter  = rep(rnorm(D, 0, 500), each = N_i),
  npar   = rep(rnorm(D, 0, 500), each = N_i),
  ipar   = rep(rnorm(D, 0, 1000), each = N_i)
)
eps_list = list(
  normal = rnorm(N, 0, 1000),
  inter  = rnorm(N, 0, 1000),
  npar   = rpareto(N, 3, 800),
  ipar   = rpareto(N, 3, 800)
)

hist(rpareto(50, 3, 800))

# ==== Study Variables ====
Y = data.table(
  Y_Normal      = 5000 - 500 * X_Normal$X1 - 500 * X_Normal$X2 + v_list$normal + eps_list$normal,
  Y_Interaction = 15000 - 500 * X_Interaction$X1 * X_Interaction$X2 - 250 * X_Interaction$X2^2 + v_list$inter + eps_list$inter,
  Y_Normal_Par  = 5000 - 500 * X_NPar$X1 - 500 * X_NPar$X2 + v_list$npar + eps_list$npar,
  Y_Inter_Par   = 20000 - 500 * X_IPAR$X1 * X_IPAR$X2 - 250 * X_IPAR$X2^2 + v_list$ipar + eps_list$ipar
)

# ==== Auxiliary Variables ====
V = as.data.table(replicate(100, runif(N, -1, 1)))
setnames(V, paste0("V_", 1:100))

# ==== Final Population Dataset ====
pop = data.table(Domain = domain, weight = runif(N, 0.5, 1.5))
pop = cbind(pop,
             X1_Normal = X_Normal$X1, X2_Normal = X_Normal$X2,
             X1_Inter  = X_Interaction$X1, X2_Inter = X_Interaction$X2,
             X1_NPar   = X_NPar$X1, X2_NPar = X_NPar$X2,
             X1_IPAR   = X_IPAR$X1, X2_IPAR = X_IPAR$X2,
             Y, V)
colnames(pop)
# ==== Save to CSV ====
##fwrite(pop, "C:/Users/user/Desktop/Lab/Hamza Paper MERF/RESULTS/Hamza MERF Stunting April 25/codes/Simulated data/simulated_population_data.csv")
##cat("✅ Data saved at your specified folder.\n")
