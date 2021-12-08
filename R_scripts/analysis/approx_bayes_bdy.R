library(magrittr)
load_all("~/research_code/glmeiv/")
library(dplyr)

########
# 1. XIE
########
# load the results
xie_result_dir <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/xie/results/")
glmeiv_res <- paste0(xie_result_dir, "glmeiv_result.rds") %>% readRDS()
med_est <- glmeiv_res %>% filter(type == "cis", target == "estimate") %>% select(parameter, value, pair_id) %>%
  tidyr::pivot_wider(data = ., id_cols = "pair_id", names_from = "parameter", values_from = "value") %>% 
  select(-pair_id) %>% summarize_all(median) %>% select(., starts_with("g"), "pi")

# Load the data
xie_data_dir <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/xie/data/")
covariate_matrix <- paste0(xie_data_dir, "covariate_matrix.rds") %>% readRDS()
g_offset <- paste0(xie_data_dir, "g_offset.rds") %>% readRDS()
g_fam <- poisson() %>% augment_family_object()

# set arguments to Bayes-opt function
g_intercept <- med_est$g_intercept %>% log()
g_pert <- med_est$g_perturbation %>% log()
g_coef <- med_est[, c("g_batchbatch_2", "g_batchbatch_3", "g_batchbatch_4", "g_batchbatch_5", "g_p_mito")] %>% as.numeric()
pi <- med_est$pi

# compute median bayes-optimal decision boundary
bdy <- get_optimal_threshold(g_intercept = g_intercept, g_perturbation = g_pert, g_fam = g_fam, pi = pi,
                      covariate_matrix = NULL, g_covariate_coefs = NULL, g_offset = median(g_offset)) %>% round()
f_name <- paste0(xie_data_dir, "dec_bdy.txt")
if (file.exists(f_name)) file.remove(f_name)
file_con <- file(f_name)
writeLines(text = as.character(bdy), con = file_con)
close(file_con)

##############
# 2. GASPERINI
##############
gasp_result_dir <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/gasperini/results/")
glmeiv_res <- paste0(gasp_result_dir, "result_glmeiv.rds") %>% readRDS()
med_est <- glmeiv_res %>% filter(site_type == "selfTSS", target == "estimate") %>% select(parameter, value, pair_id) %>%
  tidyr::pivot_wider(data = ., id_cols = "pair_id", names_from = "parameter", values_from = "value") %>% 
  select(-pair_id) %>% summarize_all(median) %>% select(., starts_with("g"), "pi")

# Load the data
gasp_data_dir <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/gasperini/data/")
covariate_matrix <- paste0(gasp_data_dir, "covariate_matrix.rds") %>% readRDS()
g_offset <- paste0(gasp_data_dir, "g_offsets.rds") %>% readRDS()
g_fam <- poisson() %>% augment_family_object()

# set arguments to Bayes-opt function
g_intercept <- med_est$g_intercept %>% log()
g_pert <- med_est$g_perturbation %>% log()
pi <- med_est$pi

# compute bayes-opt dec bdy and save
bdy <- get_optimal_threshold(g_intercept = g_intercept, g_perturbation = g_pert, g_fam = g_fam, pi = pi,
                             covariate_matrix = NULL, g_covariate_coefs = NULL, g_offset = median(g_offset)) %>% round()

f_name <- paste0(gasp_data_dir, "dec_bdy.txt")
if (file.exists(f_name)) file.remove(f_name)
file_con <- file(f_name)
writeLines(text = as.character(bdy), con = file_con)
close(file_con)
