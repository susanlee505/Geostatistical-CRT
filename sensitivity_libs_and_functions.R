library(INLA)
inla.setOption(inla.mode="experimental")
library(fields)  
library(MASS)    
library(sp)
library(ggplot2)
library(dplyr)
library(tidyr)

# Define global parameters
gamma_true <- 0.1  # Fixed covariate effect
delta_true <- 0.1  # Interaction effect

phi_true <- c(3.5) # Decay parameter for exponential covariance. should be consider the grid size.

ICC <- c(0.05)
#ICC_cluster = 0.02

# For sample size calculation
m <- 40 # cluster size
desired_power <- 0.85
desired_theta <- 0.6 # 0.4*sqrt(sigma2W_true)
sigma2W_true <- (desired_theta/0.4)^2; sigma2W_true
z_1minus_beta <- qnorm(desired_power)
beta_true <- c(0, 0.3, 0.6)

#'@prop sigma2B/(sigma2B + tau2)
compute_variances <- function(ICC, prop=0.3) {
  denom <- (1 / ICC) - 1 - (1 - prop) / prop
  sigma2B <- sigma2W_true / denom
  tau2 <- (1 - prop) / prop * sigma2B
  ICC_total <- (sigma2B + tau2) / (sigma2B + tau2 + sigma2W_true)
  return(data.frame(
    ICC = ICC,
    ICC_total= ICC_total,
    sigma2B_true = sigma2B,
    tau2_true = tau2 )) 
}
variance_df <- do.call(rbind, lapply(ICC, function(x) compute_variances(x, prop=0.5))); variance_df

compute_design_params <- function(icc) {
  VIF <- 1 + (m - 1) * icc
  n_per_group <- 2 * sigma2W_true * (z_1minus_beta + 1.96)^2 * VIF / desired_theta^2
  total_clusters <- ceiling(n_per_group * 2 / m)
  grid_size <- floor(sqrt(total_clusters))
  return(data.frame(VIF = VIF, ICC = icc, total_clusters = total_clusters, grid_size = grid_size))
}

design_df <- compute_design_params(ICC); design_df
design_df$grid_size <- 4

# Create full scenario grid
scenario_grid <- expand.grid(
  kernel = c("exp","matern"),
  prior_strength = c("weak","moderate","informative"),
  ICC = ICC,
  phi_true = phi_true,
  beta_true = beta_true
) %>%
  left_join(variance_df, by = "ICC") %>%
  left_join(design_df, by = "ICC") %>%
  mutate(
    m = m,
    scenario_id = row_number()
  )
head(scenario_grid)


'
m=40; grid_size = 4; sigma2B_true = 0.1; phi_true = 1; beta_true = 0.3; tau2_true = 0.1
'
SimulateData <- function(m, grid_size, sigma2B_true, phi_true, tau2_true, beta_true) {
  
  cluster_grid <- expand.grid(x = 1:grid_size, y = 1:grid_size)
  n_clusters <- nrow(cluster_grid)
  cluster_grid$cluster_id <- 1:n_clusters
  
  # Checkerboard treatment assignment
  #cluster_grid$treatment <- (cluster_grid$x + cluster_grid$y + 1) %% 2
  treatment_assign <- c(rep(0, floor(n_clusters/2)), rep(1, ceiling(n_clusters/2)))
  treatment_assign <- sample(treatment_assign)
  cluster_grid$treatment <- treatment_assign
  
  # Simulate data for each cluster
  data <- do.call(rbind, lapply(1:n_clusters, function(j) {
    base_locs <- matrix(runif(m * 2, min = -0.5, max = 0.5), ncol = 2)
    x_coords <- cluster_grid$x[j] + base_locs[,1]
    y_coords <- cluster_grid$y[j] + base_locs[,2]
    treatment <- rep(cluster_grid$treatment[j], m)
    
    data.frame(cluster_id = cluster_grid$cluster_id[j],
               x = x_coords,
               y = y_coords,
               treatment = treatment) }))
  # # --- SPDE Mesh and simulation ---
  coords <- as.matrix(data[, c("x", "y")])
  mesh <- inla.mesh.2d(loc = coords, max.edge = c(0.2, 0.5), cutoff = 0.05)
  #plot(mesh)
  
  spde_model <- inla.spde2.matern(mesh = mesh, alpha = 1.5) #exponential covariance: alpha=1.5
  A_eta <- inla.spde.make.A(mesh = mesh, loc = coords)
  tau <- sqrt(1/tau2_true); kappa <- 2/phi_true
  Q_eta <- inla.spde2.precision(spde_model, theta = c(log(tau), log(kappa)))
  spde_model$param.inla
  # Simulate spatial field w ~ N(0, Q^-1)
  cholQ <- Cholesky(Q_eta, LDL = FALSE)
  w_field <- as.vector(Matrix::solve(cholQ, rnorm(nrow(Q_eta))))
  w_ij <- as.vector(A_eta %*% w_field)
  
  # Covariates and noise
  data$x_ij <- rnorm(nrow(data), mean = 0, sd = 1)
  cluster_reff <- rnorm(n_clusters, mean = 0, sd = sqrt(sigma2B_true))
  u_i <- cluster_reff[data$cluster_id]
  e_ij <- rnorm(nrow(data), mean = 0, sd = sqrt(sigma2W_true))
  
  mu_ij <- beta_true * data$treatment +
    gamma_true * data$x_ij +
    delta_true * data$treatment * data$x_ij
  
  data$y_data <- mu_ij + u_i + w_ij + e_ij
  return(data)
}


FitModels <- function(data, Delta = 0, ICC_true, beta_true, kernel, prior_strength) {
  
  df <- data %>% group_by(cluster_id) %>% mutate(x_ij = x_ij - mean(x_ij)) %>% ungroup()
  df$cluster_id <- relevel(factor(df$cluster_id), ref = "1") 
  
  # ---- Mesh ----
  non_convex_bdry <- inla.nonconvex.hull(points = cbind(df$x, df$y), -0.03, -0.05, resolution = c(100, 100))
  mesh <- inla.mesh.2d(boundary = non_convex_bdry, max.edge = c(0.2, 0.5), cutoff = 0.05)
  # Project continuous spatial field onto the observed data locatiosns
  A_df <- inla.spde.make.A(mesh, loc = as.matrix(df[, c("x", "y")]))
  spatial_index <- inla.spde.make.index(name = "s", n.spde = mesh$n)
  
  df_stack_SMM <- inla.stack(
    data = list(y_data = df$y_data),
    A = list(1, A_df),
    effects = list(df[, c("cluster_id", "treatment", "x_ij")], spatial_index),
    tag = "smm_model" )
  
  # ---- SPDE helper (PC-Matérn) ----
  mk_spde <- function(mesh, alpha, range_pair, sigma_pair) {
    inla.spde2.pcmatern( mesh = mesh, alpha = alpha, prior.range = range_pair, prior.sigma = sigma_pair ) }
  
  alpha <- switch(as.character(kernel), exp = 1.5, matern = 2.0 )
  
  prior_pairs <- switch(
    as.character(prior_strength),
    "weak" = list(range = c(9.8995, 0.5), sigma = c(3, 0.1)),
    "moderate" = list(range = c(4.95, 0.5), sigma = c(0.5, 0.5)),
    "informative" = list(range = c(4.95, 0.8), sigma = c(0.5, 0.2)) )
  
  model_spde <- mk_spde(mesh, alpha, prior_pairs$range, prior_pairs$sigma)
  
  # ---- Priors ----
  # Likelihood (residual) precision: PC prior with P(sigma_W > 10) = 0.1  (very weak on standardized scale)
  hyperparameters <- list(
    prec.family  = list(prior = "pc.prec", param = c(10, 0.1)),
    prec.cluster = list(prior = "pc.prec", param = c(3, 0.1))   # iid cluster effect precision
  )
  
  # Fixed effects: Normal(0, Var = 1000) ⇒ prec = 1/1000
  control.fixed <- list(mean = list(default = 0), prec = list(default = 1/1000))
  
  # ---- Fit CRT-SMM ----
  CRT_SMM <- inla(
    formula = y_data ~ treatment + x_ij + treatment:x_ij +
      f(cluster_id, model = "iid", constr = TRUE, hyper = list(prec = hyperparameters$prec.cluster)) +
      f(s, model = model_spde),
    data = inla.stack.data(df_stack_SMM),
    family = "gaussian",
    control.predictor = list(A = inla.stack.A(df_stack_SMM), compute = FALSE),
    control.fixed = control.fixed,
    control.family = list(hyper = list(prec = hyperparameters$prec.family)),
    control.inla = list(strategy = "adaptive", int.strategy = "eb"),
    control.compute = list(config = TRUE)
  )
  
  # ---- Posterior measures ----
  # Pr(theta > Delta)
  prob_rej.SMM <- 1 - inla.pmarginal(Delta, CRT_SMM$marginals.fixed$treatment)
  
  # Coverage for theta
  CI_theta <- inla.qmarginal(c(0.025, 0.975), CRT_SMM$marginals.fixed$treatment)
  theta.covered.SMM <- (beta_true >= CI_theta[1]) & (beta_true <= CI_theta[2])
  
  # Coverage for ICC
  sigma2B_samples <- 1 / inla.rmarginal(1000, CRT_SMM$marginals.hyperpar[["Precision for cluster_id"]])
  tau2_samples <- inla.rmarginal(1000, CRT_SMM$marginals.hyperpar[["Stdev for s"]])^2
  sigma2W_samples <- 1 / inla.rmarginal(1000, CRT_SMM$marginals.hyperpar[["Precision for the Gaussian observations"]])
  
  ICC_samples <- (sigma2B_samples) /
    (sigma2B_samples + tau2_samples + sigma2W_samples)
  ICC_CI <- quantile(ICC_samples, probs = c(0.025, 0.975))
  ICC_covered.SMM <- ICC_true >= ICC_CI[1] & ICC_true <= ICC_CI[2]
  
  # Bias/MSE components
  Ebeta.SMM <- CRT_SMM$summary.fixed["treatment", "mean"]
  SE.SMM    <- CRT_SMM$summary.fixed["treatment", "sd"]
  
  data.frame(
    prob_rej.SMM = prob_rej.SMM,
    theta.covered.SMM = theta.covered.SMM,
    ICC_covered.SMM = ICC_covered.SMM,
    Ebeta.SMM = Ebeta.SMM,
    SE.SMM = SE.SMM
  )
}