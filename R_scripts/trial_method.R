library(ondisc)
library(glmeiv)
library(magrittr)
# load mRNA, gRNA, and covariate data

# directories
gasp_offsite_dir <- paste0(.get_config_path("LOCAL_GASPERINI_2019_DATA_DIR"), "at-scale/processed/")
glmeiv_data_dir <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/gasperini/data/")
# genes
gene_odm_fp <- paste0(gasp_offsite_dir, "gene/gasp_scale_gene_expressions.odm")
gene_metadata_fp <- paste0(glmeiv_data_dir, "gene_qc_metadata.rds")
# gRNAs
gRNA_odm_fp <- paste0(gasp_offsite_dir, "gRNA_grouped/gasp_scale_gRNA_counts_grouped.odm")
gRNA_metadata_fp <- paste0(glmeiv_data_dir, "gRNA_qc_metadata.rds")
# covariate matrix, pairs
covariate_matrix_fp <- paste0(glmeiv_data_dir, "covariate_matrix.rds")
gRNA_gene_pairs_fp <- paste0(glmeiv_data_dir, "gRNA_gene_pairs.rds")

# Load ODMs, covariate matrix, pairs
gene_odm <- read_odm(odm_fp = gene_odm_fp, metadata_fp = gene_metadata_fp)
gRNA_odm <- read_odm(odm_fp = gRNA_odm_fp, metadata_fp = gRNA_metadata_fp)
global_covariate_matrix <- readRDS(covariate_matrix_fp) %>% dplyr::mutate(batch = ifelse(batch == "prep_batch_1", 0, 1))
gRNA_gene_pairs <- readRDS(gRNA_gene_pairs_fp)
gRNA_gene_pairs %>% dplyr::filter(site_type == "DHS")

# get gene and gRNA data
# set.seed(5)
# my_gene <- gene_odm %>% mutate_feature_covariates(gene_id = get_feature_ids(gene_odm)) %>%
#  get_feature_covariates() %>% dplyr::filter(mean_expression > 15) %>% dplyr::sample_n(1) %>% dplyr::pull(gene_id)
# my_gRNA <- gRNA_odm %>% mutate_feature_covariates(gRNA_id = get_feature_ids(gRNA_odm)) %>%
#  get_feature_covariates() %>% dplyr::filter(coef_of_variation < 15) %>% dplyr::sample_n(1) %>% dplyr::pull(gRNA_id)

m <- as.numeric(gene_odm[["ENSG00000169813",]])
g <- as.numeric(gRNA_odm[["chr10.1340_top_two",]])

m_lib_size <- exp(global_covariate_matrix$lg_mRNA_lib_size)
g_lib_size <-  exp(global_covariate_matrix$lg_gRNA_lib_size)

m_precomp <- glmeiv::run_glmeiv_precomputation(y = m,
                                               covariate_matrix = dplyr::select(global_covariate_matrix, batch, p_mito),
                                               offset = log(m_lib_size),
                                               fam = augment_family_object( MASS::negative.binomial(NA) ))
g_precomp <- glmeiv::run_glmeiv_precomputation(y = g,
                                               covariate_matrix = dplyr::select(global_covariate_matrix, batch, p_mito),
                                               offset = log(g_lib_size),
                                               fam = augment_family_object( poisson() ))
if (FALSE) {
fit <- run_glmeiv_given_precomputations(m = m, g = g, m_precomp = m_precomp, g_precomp = g_precomp,
                                        covariate_matrix = dplyr::select(global_covariate_matrix, batch, p_mito),
                                        m_offset = dplyr::pull(global_covariate_matrix, lg_mRNA_lib_size),
                                        g_offset = dplyr::pull(global_covariate_matrix, lg_gRNA_lib_size),
                                        n_em_rep = 10, pi_guess_range = c(1e-5, 0.03),
                                        m_perturbation_guess_range = log(c(0.1, 1.5)),
                                        g_perturbation_guess_range = log(c(0.5, 10)))
} else {
  m_fam <- m_precomp$fam; g_fam <- g_precomp$fam
  covariate_matrix = dplyr::select(global_covariate_matrix, batch, p_mito);
  m_offset = dplyr::pull(global_covariate_matrix, lg_mRNA_lib_size);
  g_offset = dplyr::pull(global_covariate_matrix, lg_gRNA_lib_size);
  n_em_rep = 10; pi_guess_range = c(1e-5, 0.03);
  m_perturbation_guess_range = log(c(0.1, 1.5));
  g_perturbation_guess_range = log(c(5, 10));
  
  # get random starting guesses
  guesses <- lapply(list(pi = pi_guess_range,
                         m_perturbation = m_perturbation_guess_range,
                         g_perturbation = g_perturbation_guess_range), function(r) {
                           stats::runif(n = n_em_rep, min = r[1], max = r[2])
                         })
  
  # fit the reduced GLM-EIV model over the random starting vectors
  exp_m_offset <- m_precomp$fitted_values; exp_g_offset <- g_precomp$fitted_values
  reduced_fits <- lapply(seq(1L, n_em_rep), function(i) {
    run_reduced_em_algo(m = m, g = g, exp_m_offset = exp_m_offset, exp_g_offset = exp_g_offset,
                        m_pert_guess = guesses$m_perturbation[i],
                        g_pert_guess = guesses$g_perturbation[i],
                        pi_guess = guesses$pi[i], m_fam = m_fam, g_fam = g_fam)
  })
  
  # obtain the best run according to log-likelihood
  best_run <- select_best_em_run(reduced_fits)
  
  # Finally, run full GLM-EIV using pilot estimates
  fit <- run_full_glmeiv_given_pilot_params(m = m, g = g, m_fam = m_fam, g_fam = g_fam,
                                            pi_guess = best_run$pi, m_intercept_guess = m_precomp$fitted_intercept,
                                            m_perturbation_guess = best_run$m_perturbation, m_covariate_coefs_guess = m_precomp$covariate_coefs,
                                            g_intercept_guess = g_precomp$fitted_intercept, g_perturbation_guess = best_run$g_perturbation,
                                            g_covariate_coefs_guess = g_precomp$covariate_coefs, covariate_matrix = covariate_matrix,
                                            m_offset = m_offset, g_offset = g_offset)
}

s <- run_inference_on_em_fit(fit)
run_delta_method(s$estimate[s$variable == "m_perturbation"],
                 s$std_error[s$variable == "m_perturbation"])
table(phat = g >= 5, em = fit$posterior_perturbation_probs > 0.5)
sum(g >= 5); sum(fit$posterior_perturbation_probs > 0.9)
