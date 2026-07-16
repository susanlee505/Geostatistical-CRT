source("sensitivity_libs_and_functions.R")

library(foreach)
library(doParallel)
library(MASS)
library(Matrix)

# ==== Parse scenario index from command line ====
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("No scenario index provided")
row <- as.integer(args[1])

# ==== Load scenario ====
current_scenario <- scenario_grid[row, ]
scenario_id <- current_scenario$scenario_id
kernel <- current_scenario$kernel
prior_strength <- current_scenario$prior_strength
icc <- current_scenario$ICC
sigma2B <- current_scenario$sigma2B_true
phi <- current_scenario$phi_true
tau2 <- current_scenario$tau2_true
beta <- current_scenario$beta_true
m <- current_scenario$m
grid_size <- current_scenario$grid_size

n_sim <- 10000

# ==== Create template for NA fallback ====
set.seed(1234)
tmp <- SimulateData(m=40, grid_size, sigma2B, phi, tau2, beta)
template <- FitModels(tmp, Delta = 0, ICC_true = icc, beta_true = beta, kernel = "exp", prior_strength = "weak")

# ==== Set up parallel backend ====
num_cores <- min(10, parallel::detectCores() - 1)
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# ==== Run simulations in parallel ====
results_list <- foreach(
  i = 1:n_sim,
  .combine = rbind,
  .packages = c("INLA", "MASS", "Matrix", "dplyr")
) %dopar% {
  data <- SimulateData(m, grid_size, sigma2B, phi, tau2, beta)
  
  fit_result <- tryCatch({
    FitModels(data, Delta = 0, icc, beta, kernel, prior_strength)
  }, error = function(e) NULL)
  
  if (!is.null(fit_result)) {
    fit_result$scenario_id <- scenario_id
    fit_result$kernel <- kernel
    fit_result$prior_strength <- prior_strength
    fit_result$ICC <- icc
    fit_result$phi <- phi
    fit_result$sigma2B <- sigma2B
    fit_result$tau2 <- tau2
    fit_result$beta <- beta
    return(fit_result)
  } else {
    na_row <- template[1,]
    na_row[] <- NA
    na_row$scenario_id <- NA
    na_row$kernel <- NA
    na_row$prior_strength <- NA
    na_row$ICC <- NA
    na_row$phi <- NA
    na_row$sigma2B <- NA
    na_row$tau2 <- NA
    na_row$beta <- NA
    return(na_row)
  } 
}

# ==== Save output ====
write.csv(results_list, file = sprintf("output_sensitivity/scenario_id%d.csv", scenario_id), row.names = FALSE)

# ==== Cleanup ====
stopCluster(cl)