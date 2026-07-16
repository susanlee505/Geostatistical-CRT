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

#phi_true <- c(5) # Decay parameter for exponential covariance. should be consider the grid size.
phi_true <- c(1.5, 3.5)
ICC <- c(0.05, 0.15, 0.25)

# For sample size calculation
m <- 40 # cluster size
desired_power <- 0.85
desired_theta <- 0.6 # 0.4*sqrt(sigma2W_true)
sigma2W_true <- (desired_theta/0.4)^2; sigma2W_true
z_1minus_beta <- qnorm(desired_power)
beta_true <- seq(0, 1.4, by=0.1); beta_true 

#'@prop sigma2B/(sigma2B + tau2)
compute_variances <- function(ICC, prop) {
  denom <- (1 / ICC) - 1 - (1 - prop) / prop
  sigma2B <- sigma2W_true / denom
  tau2 <- (1 - prop) / prop * sigma2B
  ICC_total <- (sigma2B + tau2) / (sigma2B + tau2 + sigma2W_true)
  return(data.frame(
    ICC = ICC,
    sigma2B_true = sigma2B,
    tau2_true = tau2 )) 
}
variance_df <- do.call(rbind, lapply(ICC, function(x) compute_variances(x, prop=0.3))); variance_df

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


FitModels <- function(data, Delta = 0, ICC_true, beta_true) {
  df <- data %>% group_by(cluster_id) %>% mutate(x_ij = x_ij - mean(x_ij)) %>% ungroup()
  df$cluster_id <- relevel(factor(df$cluster_id), ref = "1") 
  df_cl <- data %>% group_by(cluster_id) %>% 
    summarise(y_bar = mean(y_data), x_bar = mean(x_ij), treatment = first(treatment))
  
  # Build spatial mesh (non-convex hull approach to avoid adding many small triangles outside the domain of interest)
  non_convex_bdry <- inla.nonconvex.hull(points = cbind(df$x, df$y), -0.03, -0.05, resolution = c(100, 100))
  mesh <- inla.mesh.2d(boundary = non_convex_bdry, max.edge = c(0.2, 0.5), cutoff = 0.05)
  #plot(mesh)
  #points(df$x, df$y, pch = 1, col = "green")
  #--------------
  # Build SPDE
  #--------------
  # Build SPDE model for computational efficiency
  model_spde <- inla.spde2.pcmatern(mesh = mesh, alpha = 1.5, 
                                    # Prior for range for rho = 1.41*phi_true: p(rho < 20) = 0.5
                                    prior.range = c(9.8995, 0.5), 
                                    # Prior for spatial variance: p(tau > 3) = 0.1
                                    prior.sigma = c(3, 0.1))
  
  # must be created to include spatial effects
  spatial_index <- inla.spde.make.index(name = "s", n.spde = model_spde$n.spde)
  # Project continuous spatial field onto the observed data locatiosns
  A_df <- inla.spde.make.A(mesh, loc = as.matrix(df[, c("x", "y")]))
  
  dim(A_df) # num_obs, df.spde$n.spde * n_clusters
  any(is.na(A_df))  # Should be FALSE
  any(A_df < 0)  # Should be FALSE
  #-------------
  # Build INLA
  #-------------
  df_stack_SMM <- inla.stack(
    data = list(y_data = df$y_data),
    A = list(1, A_df),
    effects = list(df[, c("cluster_id", "treatment", "x_ij")], spatial_index),
    tag = "smm_model"
  )
  dim(inla.stack.A(df_stack_SMM)) # simplified projector matrix, where each column holds one effect block
  names(inla.stack.data(df_stack_SMM))
  
  # Define hyperparameters
  hyperparameters <- list(
    # PC prior: p(sigma > 1)=0.01 => c(1, 0.01)--explicitly for random effect
    # sqrt(scenario_grid$sigma2W_true)
    prec.family = list(prior = "pc.prec", param = c(10, 0.1)), 
    # sqrt(scenario_grid$sigma2B_true)
    prec.cluster = list(prior = "pc.prec", param = c(3, 0.1))
  )
  control.fixed <- list(
    mean = list(default = 0),    # normal priors for fixed effects
    prec = list(default = 0.01)
  )
  #-------------------------------------------------------------
  # Fit CRT-cluster model
  CRT_cluster <- inla(
    formula = y_bar ~ treatment + x_bar + treatment:x_bar,
    data = df_cl,
    family = "gaussian",
    control.predictor = list(compute = FALSE),
    control.fixed = control.fixed,
    control.family = list(hyper = list(prec = hyperparameters$prec.family)),
    control.inla = list(strategy = "adaptive", int.strategy = "eb"),
    control.compute = list(config = FALSE)
  )
  # Fit CRT-FM_naive
  CRT_FM_naive <- inla(
    formula = y_data ~ treatment + x_ij + treatment:x_ij,
    data = df,
    family = "gaussian",
    control.fixed = control.fixed,
    control.family = list(hyper = list(prec = hyperparameters$prec.family)),
    control.inla = list(strategy = "adaptive", int.strategy = "eb"),
    control.compute = list(config = FALSE)
  )
  
  # Fit CRT-FM
  CRT_FM <- inla(
    formula = y_data ~ treatment + x_ij + treatment:x_ij + factor(cluster_id),
    data = df,
    family = "gaussian",
    control.fixed = list(
      mean = list(default = 0),    # normal priors for fixed effects
      prec = list(default = 1)
    ),
    control.family = list(hyper = list(prec = hyperparameters$prec.family)),
    control.inla = list(strategy = "adaptive", int.strategy = "eb"),
    control.compute = list(config = FALSE)
  )
  summary(CRT_FM)
  
  # Fit CRT-MM
  CRT_MM <- inla(
    formula = y_data ~ treatment  + x_ij + treatment:x_ij + f(cluster_id, model = "iid", constr = FALSE, hyper = list(prec = hyperparameters$prec.cluster)),
    #constr = TRUE: constriants of sum(u_j)=0
    data = df,
    family = "gaussian",
    control.fixed = control.fixed,
    control.family = list(hyper = list(prec = hyperparameters$prec.family)),
    control.inla = list(strategy = "adaptive", int.strategy = "eb"),
    control.compute = list(config = FALSE)
  )
  # Fit CRT-SMM
  CRT_SMM <- inla(
    formula = y_data ~ treatment + x_ij + treatment:x_ij + 
      f(cluster_id, model = "iid", constr = TRUE, hyper = list(prec = hyperparameters$prec.cluster)) + f(s, model = model_spde),
    data = inla.stack.data(df_stack_SMM, spde = model_spde),
    family = "gaussian",
    control.predictor = list(A = inla.stack.A(df_stack_SMM), compute = FALSE),
    control.fixed = control.fixed,
    control.family = list(hyper = list(prec = hyperparameters$prec.family)),
    control.inla = list(strategy = "adaptive", int.strategy = "eb"),
    control.compute = list(config = TRUE)
  )
  # ==== plot ====  
  '
png("images\\density_plot.png", width = 800, height = 600)

plot(CRT_FM_naive$marginals.fixed$treatment, col = "skyblue", type = "l",
     main = "Posterior of Treatment Effect",
     xlab = expression("Treatment Effect (" * hat(theta) * ")"),
     ylab = "Density", xlim = c(-3, 3))
lines(CRT_FM$marginals.fixed$treatment, col = "purple")  
lines(CRT_MM$marginals.fixed$treatment, col = "blue")
lines(CRT_SMM$marginals.fixed$treatment, col = "red", lty = 2,lwd = 3)
lines(CRT_cluster$marginals.fixed$treatment, col = "green", lty = 3, lwd = 2) 
legend("topright",
       legend = c("FM-naive", "FM", "MM", "SMM", "Cluster"),
       col = c("skyblue","purple", "blue", "red",  "green"),
       lty = c(1, 1, 1, 2, 3))
       
dev.off()
'
  # ==== Compute Pr(theta>Delta) ====
  prob_rej.cluster <- 1 - inla.pmarginal( Delta, CRT_cluster$marginals.fixed$treatment )
  prob_rej.FM_naive <- 1 - inla.pmarginal( Delta, CRT_FM_naive$marginals.fixed$treatment )
  prob_rej.FM <- 1 - inla.pmarginal( Delta, CRT_FM$marginals.fixed$treatment )
  prob_rej.MM <- 1 - inla.pmarginal( Delta, CRT_MM$marginals.fixed$treatment ) # e.g., beta_i<=delta
  prob_rej.SMM <- 1 - inla.pmarginal( Delta, CRT_SMM$marginals.fixed$treatment) 
  
  # ==== coverage probability for theta ====
  CI_theta <-inla.qmarginal(c(0.025, 0.975), CRT_cluster$marginals.fixed$treatment)
  theta.covered.cluster <- beta_true >= CI_theta[1] & beta_true <= CI_theta[2]
  
  CI_theta <-inla.qmarginal(c(0.025, 0.975), CRT_FM_naive$marginals.fixed$treatment)
  theta.covered.FM_naive <- beta_true >= CI_theta[1] & beta_true <= CI_theta[2]
  
  CI_theta <-inla.qmarginal(c(0.025, 0.975), CRT_FM$marginals.fixed$treatment)
  theta.covered.FM <- beta_true >= CI_theta[1] & beta_true <= CI_theta[2]
  
  CI_theta <- inla.qmarginal(c(0.025, 0.975), CRT_MM$marginals.fixed$treatment)
  theta.covered.MM <- beta_true >= CI_theta[1] & beta_true <= CI_theta[2]
  
  CI_theta <- inla.qmarginal(c(0.025, 0.975), CRT_SMM$marginals.fixed$treatment)
  theta.covered.SMM <- beta_true >= CI_theta[1] & beta_true <= CI_theta[2]
  
  # ==== coverage probability for ICC ====
  sigma2B_samples <- 1 / inla.rmarginal(1000, CRT_MM$marginals.hyperpar[["Precision for cluster_id"]])
  sigma2W_samples <- 1 / inla.rmarginal(1000, CRT_MM$marginals.hyperpar[["Precision for the Gaussian observations"]])
  
  ICC_samples <- sigma2B_samples / (sigma2B_samples + sigma2W_samples)
  ICC_CI <- quantile(ICC_samples, probs = c(0.025, 0.975))
  ICC_covered.MM <- ICC_true >= ICC_CI[1] & ICC_true <= ICC_CI[2]
  #--------------------------------------------
  sigma2B_samples <- 1 / inla.rmarginal(1000, CRT_SMM$marginals.hyperpar[["Precision for cluster_id"]])
  tau2_samples <- inla.rmarginal(1000, CRT_SMM$marginals.hyperpar[["Stdev for s"]])^2
  sigma2W_samples <- 1 / inla.rmarginal(1000, CRT_SMM$marginals.hyperpar[["Precision for the Gaussian observations"]])
  
  ICC_samples <- (sigma2B_samples) /
    (sigma2B_samples + tau2_samples + sigma2W_samples)
  ICC_CI <- quantile(ICC_samples, probs = c(0.025, 0.975))
  ICC_covered.SMM <- ICC_true >= ICC_CI[1] & ICC_true <= ICC_CI[2]
  #--------------------------------------------
  ## For Bias, MSE
  Ebeta.cluster <- CRT_cluster$summary.fixed["treatment", "mean"]
  Ebeta.FM_naive <- CRT_FM_naive$summary.fixed["treatment", "mean"]
  Ebeta.FM <- CRT_FM$summary.fixed["treatment", "mean"]
  Ebeta.MM <- CRT_MM$summary.fixed["treatment", "mean"]
  Ebeta.SMM <- CRT_SMM$summary.fixed["treatment", "mean"]
  
  SE.cluster <- CRT_cluster$summary.fixed["treatment", "sd"]
  SE.FM_naive <- CRT_FM_naive$summary.fixed["treatment", "sd"]
  SE.FM <- CRT_FM$summary.fixed["treatment", "sd"]
  SE.MM <- CRT_MM$summary.fixed["treatment", "sd"]
  SE.SMM <- CRT_SMM$summary.fixed["treatment", "sd"]
  
  result <- data.frame(
    prob_rej.cluster = prob_rej.cluster,
    prob_rej.FM_naive = prob_rej.FM_naive,
    prob_rej.FM = prob_rej.FM,
    prob_rej.MM = prob_rej.MM,
    prob_rej.SMM = prob_rej.SMM,
    theta.covered.cluster = theta.covered.cluster,
    theta.covered.FM_naive = theta.covered.FM_naive,
    theta.covered.FM = theta.covered.FM,
    theta.covered.MM = theta.covered.MM,
    theta.covered.SMM = theta.covered.SMM,
    ICC_covered.MM = ICC_covered.MM,
    ICC_covered.SMM = ICC_covered.SMM,
    Ebeta.cluster = Ebeta.cluster,
    Ebeta.FM_naive = Ebeta.FM_naive,
    Ebeta.FM = Ebeta.FM,
    Ebeta.MM = Ebeta.MM,
    Ebeta.SMM = Ebeta.SMM,
    SE.cluster = SE.cluster,
    SE.FM_naive = SE.FM_naive,
    SE.FM = SE.FM,
    SE.MM = SE.MM,
    SE.SMM = SE.SMM
  )
  return(result)
}
