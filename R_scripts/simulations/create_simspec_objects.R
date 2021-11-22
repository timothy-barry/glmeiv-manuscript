library(magrittr)
library(glmeiv)
library(simulatr)

args <- commandArgs(trailingOnly = TRUE)
overwrite <- if (is.na(args[1])) TRUE else as.logical(args[1])
save_obj <- function(obj, file_path, overwrite) {
  if (!file.exists(file_path)) { # if file does not exist, save
    saveRDS(obj, file_path)
  } else { # if file does exist, save only if overwrite true
    if (overwrite) {
      saveRDS(obj, file_path)
    }
  }
}
sim_dir <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/simulations/spec_objects")

one_rep_times <- list(generate_data_function = NA_integer_,
                      thresholding = NA_integer_,
                      glmeiv_fast = NA_integer_,
                      glmeiv_slow = NA_integer_)

##############################
# Experiment 0
# Vary g_pert
# Fix distribution to Poisson
# All three methods
# One covariate (batch)
##############################
set.seed(4)
theta <- 20
n <- 150000
g_perturbation_grid <- log(seq(1, 4, 0.5))
param_grid <- expand.grid(g_perturbation = g_perturbation_grid)
param_grid$grid_id <- seq(1, nrow(param_grid))

fixed_params <- list(
  m_fam = poisson() %>% augment_family_object(),
  g_fam = poisson() %>% augment_family_object(),
  seed = 4,
  n = n,
  B = 500,
  m_intercept = log(0.01),
  m_perturbation = log(0.25),
  g_intercept = log(0.005),
  covariate_matrix = data.frame(batch = rbinom(n = n, size = 1, prob = 0.5)),
  m_covariate_coefs = log(0.9),
  g_covariate_coefs = log(1.1),
  n_processors = 20,
  alpha = 0.95,
  n_em_rep = 15,
  save_membership_probs_mult = 500,
  pi = 0.02,
  m_offset = log(rpois(n = n, lambda = 10000)),
  g_offset = log(rpois(n = n, lambda = 5000)),
  pi_guess_range = c(1e-5, 0.03),
  m_perturbation_guess_range = log(c(0.1, 1.5)),
  g_perturbation_guess_range = log(c(0.5, 10)),
  m_intercept_guess_range = log(c(1e-4, 1e-1)),
  g_intercept_guess_range = log(c(1e-4, 1e-1)),
  m_covariate_coefs_guess_range = log(c(0.25, 2)),
  g_covariate_coefs_guess_range = log(c(0.25, 2)),
  run_unknown_theta_precomputation = FALSE)

sim_spec_0 <- create_simulatr_specifier_object(param_grid = param_grid,
                                               fixed_params = fixed_params,
                                               one_rep_times = one_rep_times,
                                               methods = c("glmeiv_fast", "glmeiv_slow", "thresholding"))
# check <- simulatr::check_simulatr_specifier_object(simulatr_spec = sim_spec_0, B_in = 2, parallel = TRUE)
save_obj(obj = sim_spec_0, file_path = paste0(sim_dir, "/sim_spec_0.rds"), overwrite = overwrite)


####################################################################
# Experiment 1
# Vary g_pert, fix all other parameters
# Vary distribution (Poisson, NB (known theta), NB (estimated theta)
# One covariate: batch (bernoulli variable)
####################################################################
set.seed(4)
theta <- 20
n <- 150000
g_perturbation_grid <- log(seq(1, 4, 0.5))
param_grid <- expand.grid(g_perturbation = g_perturbation_grid,
                          fam_str = c("nb_theta_unknown", "nb_theta_known"))
param_grid$grid_id <- seq(1, nrow(param_grid))
fam_obj <- lapply(as.character(param_grid$fam_str), function(str) {
  switch(EXPR = str,
         "nb_theta_unknown" = MASS::negative.binomial(theta),
         "nb_theta_known" = MASS::negative.binomial(theta)) %>% glmeiv::augment_family_object()
})
param_grid$m_fam <- param_grid$g_fam <- fam_obj
param_grid$run_unknown_theta_precomputation <- as.character(param_grid$fam_str) == "nb_theta_unknown"

fixed_params <- list(
  seed = 4,
  n = n,
  B = 500,
  m_intercept = log(0.01),
  m_perturbation = log(0.25),
  g_intercept = log(0.005),
  covariate_matrix = data.frame(batch = rbinom(n = n, size = 1, prob = 0.5)),
  m_covariate_coefs = log(0.9),
  g_covariate_coefs = log(1.1),
  n_processors = 20,
  alpha = 0.95,
  n_em_rep = 15,
  save_membership_probs_mult = 500,
  pi = 0.02,
  m_offset = log(rpois(n = n, lambda = 10000)),
  g_offset = log(rpois(n = n, lambda = 5000)),
  pi_guess_range = c(1e-5, 0.03),
  m_perturbation_guess_range = log(c(0.1, 1.5)),
  g_perturbation_guess_range = log(c(0.5, 10)),
  m_intercept_guess_range = log(c(1e-4, 1e-1)),
  g_intercept_guess_range = log(c(1e-4, 1e-1)),
  m_covariate_coefs_guess_range = log(c(0.25, 2)),
  g_covariate_coefs_guess_range = log(c(0.25, 2)))

sim_spec_1 <- create_simulatr_specifier_object(param_grid = param_grid,
                                               fixed_params = fixed_params,
                                               one_rep_times = one_rep_times,
                                               methods = c("glmeiv_fast", "glmeiv_slow", "thresholding"))
# check <- simulatr::check_simulatr_specifier_object(simulatr_spec = sim_spec_1, B_in = 2, parallel = TRUE)
save_obj(obj = sim_spec_1, file_path = paste0(sim_dir, "/sim_spec_1.rds"), overwrite = overwrite)


########################################
# Experiment 2
# Gaussian response distribution
# Vary g_perturbation
# Two covariates: batch, library size
########################################
set.seed(4)
n <- 150000
g_perturbation_grid <- seq(0, 7)
param_grid <- data.frame(g_perturbation = g_perturbation_grid,
                         grid_id = seq(1, length(g_perturbation_grid)))
fixed_params <- list(
  m_fam = gaussian() %>% augment_family_object(),
  g_fam = gaussian() %>% augment_family_object(),
  seed = 4,
  n = n,
  B = 1000,
  m_intercept = 3,
  m_perturbation = -4,
  g_intercept = 1,
  covariate_matrix = data.frame(lib_size = rpois(n = n, lambda = 10000),
                                batch = rbinom(n = n, size = 1, prob = 0.5)),
  m_covariate_coefs = c(0.0025, 0.1),
  g_covariate_coefs = c(-0.005, 0.2),
  n_processors = 20,
  alpha = 0.95,
  n_em_rep = 25,
  save_membership_probs_mult = 500,
  pi = 0.05,
  m_offset = NULL,
  g_offset = NULL,
  pi_guess_range = c(0.0, 0.1),
  m_perturbation_guess_range = c(-6, -2),
  g_perturbation_guess_range = c(0, 8),
  m_intercept_guess_range = c(3, 3),
  g_intercept_guess_range = c(1, 1),
  m_covariate_coefs_guess_range = c(0.0025, 0.0025),
  g_covariate_coefs_guess_range = c(-0.005, -0.005),
  run_unknown_theta_precomputation = FALSE,
  exponentiate_coefs = FALSE,
  ep_tol = 1e-7)

sim_spec_2 <- create_simulatr_specifier_object(param_grid = param_grid,
                                               fixed_params = fixed_params,
                                               one_rep_times = one_rep_times,
                                               methods = c("glmeiv_fast", "thresholding"))

# check <- simulatr::check_simulatr_specifier_object(simulatr_spec = sim_spec_2, B_in = 2, parallel = TRUE)
save_obj(obj = sim_spec_2, file_path = paste0(sim_dir, "/sim_spec_2.rds"), overwrite = overwrite)
