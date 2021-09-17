library(magrittr)

# load the results
gasp_result_dir <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/gasperini/results/")
glmeiv_res <- readRDS(paste0(gasp_result_dir, "result_glmeiv.rds"))
thresh_res <- readRDS(paste0(gasp_result_dir, "result_thresholding.rds"))

# first, examine g_perturbation for glmeiv in the different categories
glmeiv_res %>% dplyr::filter(parameter == "g_perturbation",
                             target %in% c("confint_lower", "confint_upper", "estimate"),
                             site_type %in% c("DHS", "NTC", "positive_ctrl", "selfTSS")) %>% 
  dplyr::select(target, value, pair_id, site_type) %>%
  tidyr::pivot_wider(id_cols = c("pair_id" , "site_type"), names_from = "target", values_from = "value") %>%
  dplyr::group_by(site_type) %>% dplyr::summarize(m_est = mean(estimate),
                                                  m_confint_lower = mean(confint_lower),
                                                  m_confint_upper = mean(confint_upper))

# examine CI coverage rate for negative control pairs for both methods
get_negctrl_ci_coverage_rate <- function(df) {
  df %>% dplyr::filter(parameter == "m_perturbation", site_type == "NTC",
                               target %in% c("confint_lower", "confint_upper", "estimate")) %>%
    dplyr::select(target, value, pair_id) %>%
    tidyr::pivot_wider(id_cols = "pair_id", names_from = "target", values_from = "value") %>%
    dplyr::mutate(covered = confint_lower < 1 & confint_upper > 1) %>%
    dplyr::summarize(mean(covered, na.rm = TRUE))
}

get_negctrl_ci_coverage_rate(glmeiv_res)
get_negctrl_ci_coverage_rate(thresh_res %>% dplyr::filter(threshold == 5))
