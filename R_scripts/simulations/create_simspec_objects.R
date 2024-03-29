library(magrittr)
library(glmeiv)
library(simulatr)

print("Running create sim_spec objects")

args <- commandArgs(trailingOnly = TRUE)
overwrite <- TRUE # if (is.na(args[1])) TRUE else as.logical(args[1])
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

####################################################################
# Experiment 1 (main text simulation)
# Vary g_pert, fix all other parameters
# Vary distribution: Poisson, NB (known theta), NB (estimated theta)
# One covariate: batch (bernoulli variable)
####################################################################
set.seed(4)
m_perturbation <- log(0.25)
theta <- 20
n <- 50000
g_perturbation_grid <- log(seq(1, 4, 0.5))

param_grid <- expand.grid(g_perturbation = g_perturbation_grid,
                          fam_str = c("nb_theta_unknown", "nb_theta_known", "poisson"))
param_grid$grid_id <- seq(1, nrow(param_grid))
param_grid$ground_truth <- m_perturbation
fam_obj <- lapply(as.character(param_grid$fam_str), function(str) {
  switch(EXPR = str,
         "nb_theta_unknown" = MASS::negative.binomial(theta) |> glmeiv::augment_family_object(),
         "nb_theta_known" = MASS::negative.binomial(theta) |> glmeiv::augment_family_object(),
         "poisson" = poisson() |> glmeiv::augment_family_object())
})
param_grid$m_fam <- fam_obj
param_grid$run_mrna_unknown_theta_precomputation <- as.character(param_grid$fam_str) == "nb_theta_unknown"

fixed_params <- list(
  grna_duplet_rate = 0.0,
  mrna_duplet_rate = 0.0,
  g_fam = poisson() |> augment_family_object(),
  run_grna_unknown_theta_precomputation = FALSE,
  seed = 4,
  n = n,
  B = 500,
  m_intercept = log(0.01),
  m_perturbation = log(0.25),
  g_intercept = log(0.005),
  covariate_matrix = data.frame(batch = rbinom(n = n, size = 1, prob = 0.5)),
  rm_covariate = "",
  m_covariate_coefs = log(0.9),
  g_covariate_coefs = log(1.1),
  alpha = 0.95,
  n_em_rep = 15,
  save_membership_probs_mult = 1000L,
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
  exponentiate_coefs = FALSE,
  ep_tol = 1e-4)

sim_spec_1 <- create_simulatr_specifier_object(param_grid = param_grid,
                                               fixed_params = fixed_params,
                                               methods = c("glmeiv_slow", "glmeiv_fast", "thresholding", "unimodal_mixture"))
#sim_spec_1 <- create_simulatr_specifier_object(param_grid = param_grid,
#                                               fixed_params = fixed_params,
#                                               methods = c("run_replogle_method_simulatr", "thresholding"))
# check <- simulatr::check_simulatr_specifier_object(simulatr_spec = sim_spec_1, B_in = 2, parallel = TRUE)
save_obj(obj = sim_spec_1, file_path = paste0(sim_dir, "/sim_spec_1.rds"), overwrite = overwrite)

########################################
# Experiment 2
# Gaussian response distribution
# Vary g_perturbation
# Two covariates: batch, library size
########################################
set.seed(4)
n <- 50000
g_perturbation_grid <- seq(0, 7)
param_grid <- data.frame(g_perturbation = g_perturbation_grid,
                         grid_id = seq(1, length(g_perturbation_grid)))
m_perturbation <- -4
param_grid$ground_truth <- m_perturbation

fixed_params <- list(
  grna_duplet_rate = 0.0,
  mrna_duplet_rate = 0.0,
  run_mrna_unknown_theta_precomputation = FALSE,
  run_grna_unknown_theta_precomputation = FALSE,
  m_fam = gaussian() |> augment_family_object(),
  g_fam = gaussian() |> augment_family_object(),
  seed = 4,
  n = n,
  B = 500,
  m_intercept = 3,
  g_intercept = 1,
  m_perturbation = m_perturbation,
  covariate_matrix = data.frame(lib_size = rpois(n = n, lambda = 10000),
                                batch = rbinom(n = n, size = 1, prob = 0.5)),
  rm_covariate = "",
  m_covariate_coefs = c(0.0025, 0.1),
  g_covariate_coefs = c(-0.005, 0.2),
  alpha = 0.95,
  n_em_rep = 25,
  save_membership_probs_mult = 1000L,
  pi = 0.05,
  m_offset = NULL,
  g_offset = NULL,
  pi_guess_range = c(0.0, 0.1),
  m_perturbation_guess_range = c(-6, -2),
  g_perturbation_guess_range = c(0, 8),
  run_unknown_theta_precomputation = FALSE,
  exponentiate_coefs = FALSE,
  ep_tol = 1e-7)

sim_spec_2 <- create_simulatr_specifier_object(param_grid = param_grid,
                                               fixed_params = fixed_params,
                                               methods = c("glmeiv_fast", "thresholding", "unimodal_mixture"))
# check <- check_simulatr_specifier_object(simulatr_spec = sim_spec_2, B_in = 2)
save_obj(obj = sim_spec_2, file_path = paste0(sim_dir, "/sim_spec_2.rds"), overwrite = overwrite)

##########################################################################################
# Experiment 3: varying the size parameter
# We vary the size parameter over the grid 10^(seq(log(1, base = 10), 2, length.out = 10))
# We hold fixed other parameters. We apply NB regression (known theta)
# and NB regression (estimated theta) using both thresholding method and GLM-EIV fast
##########################################################################################
n <- 50000
m_perturbation <- log(0.25) # what happens when we set this to a different value, e.g. a value closer to 0?
thetas <- 10^(seq(log(1, base = 10), 2, length.out = 10))
m_fams <- lapply(thetas, function(theta) MASS::negative.binomial(theta) |> augment_family_object())
param_grid <- expand.grid(m_fam = m_fams,
                          run_mrna_unknown_theta_precomputation = c(TRUE, FALSE))
param_grid$theta <- sapply(param_grid$m_fam, function(fam) fam$theta)
param_grid$ground_truth <- m_perturbation
param_grid$grid_id <- seq(1L, nrow(param_grid))

fixed_params <- list(
  grna_duplet_rate = 0.0,
  mrna_duplet_rate = 0.0,
  run_grna_unknown_theta_precomputation = FALSE,
  g_fam = poisson() |> augment_family_object(),
  seed = 4,
  n = n,
  B = 500,
  g_perturbation = log(2.5),
  m_intercept = log(0.01),
  g_intercept = log(0.005),
  m_perturbation = m_perturbation,
  covariate_matrix = data.frame(batch = rbinom(n = n, size = 1, prob = 0.5)),
  rm_covariate = "",
  m_covariate_coefs = log(0.9),
  g_covariate_coefs = log(1.1),
  alpha = 0.95,
  n_em_rep = 25,
  save_membership_probs_mult = 1000L,
  pi = 0.02,
  m_offset = log(rpois(n = n, lambda = 10000)),
  g_offset = log(rpois(n = n, lambda = 5000)),
  pi_guess_range = c(1e-5, 0.1),
  m_perturbation_guess_range = log(c(0.1, 1.5)),
  g_perturbation_guess_range = log(c(0.5, 10)),
  exponentiate_coefs = FALSE,
  ep_tol = 1e-4)
sim_spec_3 <- create_simulatr_specifier_object(param_grid = param_grid,
                                               fixed_params = fixed_params,
                                               methods = c("glmeiv_fast", "thresholding", "unimodal_mixture"))
# check <- check_simulatr_specifier_object(simulatr_spec = sim_spec_3, B_in = 2)
save_obj(obj = sim_spec_3, file_path = paste0(sim_dir, "/sim_spec_3.rds"), overwrite = overwrite)

#########################################################
# Experiment 4: varying m_pert while keeping g_pert fixed
# We vary the parameter m_pert over log(.2, 1)
# We keep g_pert fixed at log(2)
##########################################################
n <- 50000
m_perturbations <- log(seq(0.2, 1, length.out = 9))
param_grid <- expand.grid(m_perturbation = m_perturbations,
                          fam_str = c("nb_theta_unknown", "nb_theta_known", "poisson"))
param_grid$grid_id <- seq(1, nrow(param_grid))
param_grid$ground_truth <- m_perturbations
fam_obj <- lapply(as.character(param_grid$fam_str), function(str) {
  switch(EXPR = str,
         "nb_theta_unknown" = MASS::negative.binomial(20) |> glmeiv::augment_family_object(),
         "nb_theta_known" = MASS::negative.binomial(20) |> glmeiv::augment_family_object(),
         "poisson" = poisson() |> glmeiv::augment_family_object())
})
param_grid$run_mrna_unknown_theta_precomputation <- as.character(param_grid$fam_str) == "nb_theta_unknown"
param_grid$m_fam <- fam_obj

fixed_params <- list(
  grna_duplet_rate = 0.0,
  mrna_duplet_rate = 0.0,
  g_fam = poisson() |> glmeiv::augment_family_object(),
  run_grna_unknown_theta_precomputation = FALSE,
  seed = 4,
  n = n,
  B = 500,
  g_perturbation = log(2.5),
  m_intercept = log(0.01),
  g_intercept = log(0.005),
  covariate_matrix = data.frame(batch = rbinom(n = n, size = 1, prob = 0.5)),
  rm_covariate = "",
  m_covariate_coefs = log(0.9),
  g_covariate_coefs = log(1.1),
  alpha = 0.95,
  n_em_rep = 25,
  save_membership_probs_mult = 1000L,
  pi = 0.02,
  m_offset = log(rpois(n = n, lambda = 10000)),
  g_offset = log(rpois(n = n, lambda = 5000)),
  pi_guess_range = c(1e-5, 0.05),
  m_perturbation_guess_range = log(c(0.1, 1.5)),
  g_perturbation_guess_range = log(c(0.5, 10)),
  exponentiate_coefs = FALSE,
  ep_tol = 1e-4)
sim_spec_4 <- create_simulatr_specifier_object(param_grid = param_grid,
                                               fixed_params = fixed_params,
                                               methods = c("glmeiv_fast", "thresholding", "unimodal_mixture"))
# check <- check_simulatr_specifier_object(simulatr_spec = sim_spec_4, B_in = 2)
save_obj(obj = sim_spec_4, file_path = paste0(sim_dir, "/sim_spec_4.rds"), overwrite = overwrite)


###########################################################################################
# Experiment 5: misspecified gRNA model
# We vary g_pert, keeping all other parameters fixed; the gRNA distribution is misspecified
# due to missing covariate and duplets
###########################################################################################
set.seed(4)
theta <- 20
n <- 50000
g_perturbation_grid <- log(seq(1, 7, 1))
m_perturbation_grid <- log(seq(1, 0.25, by = -0.25))

param_grid <- expand.grid(g_perturbation = g_perturbation_grid,
                          m_perturbation = m_perturbation_grid)
param_grid$grid_id <- seq(1, nrow(param_grid))
param_grid$ground_truth <- param_grid$m_perturbation

fixed_params <- list(
  g_fam = poisson() |> augment_family_object(),
  m_fam = MASS::negative.binomial(20) |> augment_family_object(),
  run_grna_unknown_theta_precomputation = FALSE,
  run_mrna_unknown_theta_precomputation = FALSE,
  seed = 4,
  n = n,
  B = 500,
  m_intercept = log(0.01),
  g_intercept = log(0.005),
  covariate_matrix = data.frame(batch = rbinom(n = n, size = 1, prob = 0.5),
                                cell_cycle = runif(n = n, min = 0, max = 1)),
  grna_duplet_rate = 0.01,
  mrna_duplet_rate = 0.0,
  rm_covariate = "cell_cycle",
  m_covariate_coefs = log(c(0.9, 1)),
  g_covariate_coefs = log(c(0.8, 1.25)),
  alpha = 0.95,
  n_em_rep = 15,
  save_membership_probs_mult = 1000L,
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
  exponentiate_coefs = FALSE,
  ep_tol = 1e-4)

sim_spec_5 <- create_simulatr_specifier_object(param_grid = param_grid,
                                               fixed_params = fixed_params,
                                               methods = c("glmeiv_fast", "thresholding", "unimodal_mixture"))
# check <- simulatr::check_simulatr_specifier_object(simulatr_spec = sim_spec_5, B_in = 2, parallel = TRUE)
save_obj(obj = sim_spec_5, file_path = paste0(sim_dir, "/sim_spec_5.rds"), overwrite = overwrite)

#######################################
# Experiment 6: misspecified mRNA model
#######################################
set.seed(4)
theta <- 20
n <- 50000
g_perturbation_grid <- log(seq(1.5, 4, 0.5))
m_perturbation_grid <- log(seq(1, 0.25, by = -0.25))

param_grid <- expand.grid(g_perturbation = g_perturbation_grid,
                          m_perturbation = m_perturbation_grid)
param_grid$grid_id <- seq(1, nrow(param_grid))
param_grid$ground_truth <- param_grid$m_perturbation

fixed_params <- list(
  g_fam = poisson() |> augment_family_object(),
  m_fam = MASS::negative.binomial(20) |> augment_family_object(),
  run_grna_unknown_theta_precomputation = FALSE,
  run_mrna_unknown_theta_precomputation = FALSE,
  seed = 4,
  n = n,
  B = 500,
  m_intercept = log(0.01),
  g_intercept = log(0.005),
  covariate_matrix = data.frame(batch = rbinom(n = n, size = 1, prob = 0.5),
                                cell_cycle = runif(n = n, min = 0, max = 1)),
  grna_duplet_rate = 0.00,
  mrna_duplet_rate = 0.01,
  rm_covariate = "cell_cycle",
  m_covariate_coefs = log(c(0.9, 1.25)),
  g_covariate_coefs = log(c(0.8, 1.0)),
  alpha = 0.95,
  n_em_rep = 15,
  save_membership_probs_mult = 1000L,
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
  exponentiate_coefs = FALSE,
  ep_tol = 1e-4)

sim_spec_6 <- create_simulatr_specifier_object(param_grid = param_grid,
                                               fixed_params = fixed_params,
                                               methods = c("glmeiv_fast", "thresholding", "unimodal_mixture"))
# check <- simulatr::check_simulatr_specifier_object(simulatr_spec = sim_spec_6, B_in = 2, parallel = TRUE)
save_obj(obj = sim_spec_6, file_path = paste0(sim_dir, "/sim_spec_6.rds"), overwrite = overwrite)
