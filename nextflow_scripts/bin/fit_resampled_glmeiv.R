#!/usr/bin/env Rscript

##############################
# 0. Set a few hyperparameters
##############################
pi_guess_range <- c(0.001, 0.03)
m_perturbation_guess_range <- log(c(0.1, 1.5))
g_perturbation_guess_range <- log(c(3, 15))
alpha <- 0.95
n_em_rep <- 15

########################################
# 1. Load packages and command-line args
########################################
library(magrittr)
library(ondisc)
library(glmeiv)
args <- commandArgs(trailingOnly = TRUE)
n_args <- length(args)
covariate_matrix_fp <- args[1L]
gene_odm_fp <- args[2L]
gene_metadata_fp <- args[3L]
m_offsets_fp <- args[4L]
g_offsets_fp <- args[5L] 
base_fit_fp <- args[6L] 
contam_level <- as.integer(args[7L])
B <- as.integer(args[8L])

#########################################
# 2. Load ODMs, covariate matrix, offsets
#########################################
gene_odm <- read_odm(gene_odm_fp, gene_metadata_fp)
m_offset <- readRDS(m_offsets_fp)
g_offset <- readRDS(g_offsets_fp)
covariate_matrix <- readRDS(covariate_matrix_fp)
raw_fit <- readRDS(base_fit_fp)

##############
# 3. Load data
##############
gene_id <- raw_fit$gene_id
gRNA_id <- raw_fit$gRNA_id
m <- gene_odm[[gene_id,]] %>% as.numeric()
m_fam <- raw_fit$m_precomp$fam
m_precomp <- raw_fit$m_precomp
g_fam <- poisson() %>% augment_family_object()
param_ests <- raw_fit$tbl %>% dplyr::filter(parameter %in% c("g_intercept" , "g_perturbation", "g_batch", "g_p_mito"), target == "estimate") %>%
  dplyr::select(-target)
param_ests_vect <- set_names(param_ests$value, param_ests$parameter) %>% log()
posterior_pert_probs <- raw_fit$tbl %>% dplyr::filter(target == "membership_probability") %>% dplyr::pull(value)
pi_est <- raw_fit$tbl %>% dplyr::filter(parameter == "pi", target == "estimate") %>% dplyr::pull(value)
n <- length(posterior_pert_probs)

############################
# 4. Generate synthetic data
############################
# get new_intercept and new_pert values, given background_contamination
get_new_coefs <- function(background_contamination, g_int, g_pert) {
  if (background_contamination == 0) {
    out <- c(new_g_int = g_int, new_g_pert = g_pert)
  } else {
    mean_expression_1 <- exp(g_int + g_pert)
    new_beta_0 <- log(background_contamination * mean_expression_1)
    new_beta_1 <- log(mean_expression_1) - new_beta_0
    out <- c(new_g_int = new_beta_0, new_g_pert = new_beta_1)
  }
  return(out)
}
new_coefs <- get_new_coefs(background_contamination = contam_level,
                           g_int = param_ests_vect[["g_intercept"]],
                           g_pert = param_ests_vect[["g_perturbation"]])

# generate indicators
p_mat <- sapply(X = seq(1, n), FUN = function(i) {
  rbinom(n = B, size = 1, prob = posterior_pert_probs[i])
}) %>% t()

# generate gRNA counts
g_resample <- generate_glm_data_sim(intercept = new_coefs[["new_g_int"]],
                                    perturbation_coef = new_coefs[["new_g_pert"]],
                                    perturbation_indicators = p_mat,
                                    fam = g_fam, covariate_matrix = covariate_matrix,
                                    covariate_coefs = param_ests_vect[c("g_batch", "g_p_mito")],
                                    offset = g_offset, n = n, B = B)

# set the bayes-optimal decision boundary
bdy <- get_optimal_threshold(new_coefs[["new_g_int"]], new_coefs[["new_g_pert"]], g_fam, pi_est, covariate_matrix,
                             param_ests_vect[c("g_batch", "g_p_mito")], g_offset)

###########################
# 5. fit models to the data
###########################
out_l <- vector(mode = "list", length = B)
for (i in seq(1, B)) {
  print(i)
  g <- g_resample[,i]
  # first, glmeiv
  g_precomp <- run_glmeiv_precomputation(y = g, covariate_matrix = covariate_matrix, offset = g_offset, fam = g_fam)
  fit_glmeiv <- run_glmeiv_given_precomputations(m = m, g = g, m_precomp = m_precomp, g_precomp = g_precomp,
                                                covariate_matrix = covariate_matrix, m_offset = m_offset,
                                                g_offset = g_offset, n_em_rep = n_em_rep, pi_guess_range = pi_guess_range,
                                                m_perturbation_guess_range = m_perturbation_guess_range,
                                                g_perturbation_guess_range = g_perturbation_guess_range)
  s_glmeiv <- run_inference_on_em_fit(fit_glmeiv, alpha)
  tbl_glmeiv <- wrangle_glmeiv_result(s_glmeiv, 0, fit_glmeiv, TRUE, save_membership_probs_mult = 250, 1) %>%
    dplyr::mutate(method = "glmeiv")
  
  # second, thresholding
  phat <- as.integer(g >= bdy)
  tbl_thresh <- run_thresholding_method(phat = phat, m = m, m_fam = m_fam, m_offset = m_offset, covariate_matrix = covariate_matrix,
                          n_examples_per_param = 5, alpha = alpha, exponentiate_coefs = TRUE) %>% dplyr::mutate(method = "thresholding")
  
  # combine outputs
  tbl_out <- rbind(tbl_glmeiv, tbl_thresh) %>% dplyr::mutate(run_id = i)
  out_l[[i]] <- tbl_out
}

to_save <- do.call(what = rbind, args = out_l) %>% dplyr::mutate(pair_id = paste0(gene_id, ":", gRNA_id), contam_level = contam_level) %>%
  dplyr::mutate_at(.tbl = ., .vars = c("parameter", "target", "method", "run_id", "pair_id", "contam_level"), .funs = factor)
saveRDS(object = to_save, file = "raw_result.rds")
