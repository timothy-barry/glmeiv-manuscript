#!/usr/bin/env Rscript

##############################
# 0. Set a few hyperparameters
##############################
n_em_rep <- 15
pi_guess_range <- c(0.001, 0.03)
m_perturbation_guess_range <- log(c(0.1, 1.5))
g_perturbation_guess_range <- log(c(3, 15))
alpha <- 0.95

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
gRNA_odm_fp <- args[5L]
gRNA_metadata_fp <- args[6L]
g_offsets_fp <- args[7L] 
pair_id <- args[8L] # "ENSG00000061794:MRPS35_TSS"

#########################################
# 2. Load ODMs, covariate matrix, offsets
#########################################
gene_odm <- read_odm(gene_odm_fp, gene_metadata_fp)
gRNA_odm <- read_odm(gRNA_odm_fp, gRNA_metadata_fp)
m_offset <- readRDS(m_offsets_fp)
g_offset <- readRDS(g_offsets_fp)
covariate_matrix <- readRDS(covariate_matrix_fp)
pair_id_sep <- strsplit(x = pair_id, split = ":")[[1]]
gene_id <- pair_id_sep[1]
gRNA_id <- pair_id_sep[2]

###############
# 4. Fit models
###############
m <- gene_odm[[gene_id,]] %>% as.numeric()
g <- gRNA_odm[[gRNA_id,]] %>% as.numeric()
m_fam <- MASS::negative.binomial(NA) %>% augment_family_object()
g_fam <- poisson() %>% augment_family_object()

# run precomputations
m_precomp <- run_glmeiv_precomputation(y = m, covariate_matrix = covariate_matrix, offset = m_offset, fam = m_fam)
g_precomp <- run_glmeiv_precomputation(y = g, covariate_matrix = covariate_matrix, offset = g_offset, fam = g_fam)

# fit glmeiv model
fit_glmeiv <- run_glmeiv_given_precomputations(m = m, g = g, m_precomp = m_precomp, g_precomp = g_precomp,
                                        covariate_matrix = covariate_matrix, m_offset = m_offset,
                                        g_offset = g_offset, n_em_rep = n_em_rep, pi_guess_range = pi_guess_range,
                                        m_perturbation_guess_range = m_perturbation_guess_range,
                                        g_perturbation_guess_range = g_perturbation_guess_range)
s <- run_inference_on_em_fit(fit_glmeiv, alpha)
tbl_glmeiv <- wrangle_glmeiv_result(s, 0, fit_glmeiv, TRUE, 500, 1) %>% dplyr::mutate(method = "glmeiv")

# fit thresholding method, using Bayes-optimal decision boundary as threshold
g_coef <- coef(fit_glmeiv$fit_g)
bdy <- get_optimal_threshold(g_coef[["(Intercept)"]], g_coef[["perturbation"]], g_fam, fit_glmeiv$fit_pi, covariate_matrix,
                             g_coef[c("batch", "p_mito")], g_offset)
phat <- as.integer(g >= bdy)
tbl_thresh <- run_thresholding_method(phat = phat, m = m, m_fam = m_precomp$fam, m_offset = m_offset, covariate_matrix = covariate_matrix,
                                      n_examples_per_param = 5, alpha = alpha, exponentiate_coefs = TRUE) %>% dplyr::mutate(method = "thresholding")

#################
# 5. Save outputs
#################
out_tbl <- rbind(tbl_glmeiv, tbl_thresh) %>% dplyr::mutate(run_id = -1, pair_id = pair_id, contam_level = -1) %>%
  dplyr::mutate_at(.tbl = ., .vars = c("parameter", "target", "method", "pair_id"), .funs = factor)

to_save <- list(m_precomp = m_precomp, tbl = out_tbl, gene_id = gene_id, gRNA_id = gRNA_id,
                posterior_perturbation_probs = fit_glmeiv$posterior_perturbation_probs, g_coef = coef(fit_glmeiv$fit_g), fit_pi = fit_glmeiv$fit_pi)

saveRDS(object = to_save, file = "baseline_fit.rds")