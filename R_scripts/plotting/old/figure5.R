# set figure 3 dir
fig_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "glmeiv-manuscript/figures/")
fig5_dir <- paste0(fig_dir, "fig5/")
if (!dir.exists(fig5_dir)) dir.create(fig5_dir, recursive = TRUE)

# load packages and results
xie_result_dir <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/xie/results/")
glmeiv_res <- readRDS(paste0(xie_result_dir, "result_glmeiv.rds")) %>% dplyr::distinct() # ensure ahead of time rows are distinct; check gRNA_gene_pairs data frame.
thresh_res <- readRDS(paste0(xie_result_dir, "result_thresholding.rds")) %>% dplyr::distinct()
library(ggplot2)
library(magrittr)

# examine g_perturbation for glmeiv in the different categories
glmeiv_res %>% dplyr::filter(parameter == "g_perturbation",
                             target %in% c("confint_lower", "confint_upper", "estimate"),
                             type %in% c("neg_control", "cis")) %>%
  dplyr::select(target, value, pair_id, type) %>%
  tidyr::pivot_wider(id_cols = c("pair_id", "type"),
                     names_from = "target", values_from = "value") %>%
  dplyr::group_by(type) %>% dplyr::summarize(m_est = mean(estimate),
                                             m_confint_lower = mean(confint_lower, na.rm = TRUE),
                                             m_confint_upper = mean(confint_upper, na.rm = TRUE))

# examine CI coverage rate for negative control pairs for both methods
get_negctrl_ci_coverage_rate <- function(df) {
  df %>% dplyr::filter(parameter == "m_perturbation",
                       target %in% c("confint_lower", "confint_upper", "estimate"),
                       type %in% c("neg_control")) %>%
    dplyr::select(target, value, pair_id, type) %>%
    tidyr::pivot_wider(id_cols = c("pair_id", "type"),
                       names_from = "target", values_from = "value") %>%
    dplyr::mutate(covered = confint_upper > 1 & confint_lower < 1) %>%
    dplyr::summarize(coverage_rate = mean(covered, na.rm = TRUE) * 100) 
}

get_negctrl_ci_coverage_rate(glmeiv_res)
get_negctrl_ci_coverage_rate(thresh_res)

# Compare estimates on candidate cis pairs
thresh_df <- thresh_res %>% dplyr::filter(parameter == "m_perturbation",
                                          target %in% c("estimate", "confint_lower", "confint_upper"), type == "cis") %>%
  dplyr::select(target, value, pair_id)
glmeiv_df <- glmeiv_res %>% dplyr::filter(parameter == "m_perturbation",
                                          target %in% c("estimate", "confint_lower", "confint_upper"), type == "cis") %>%
  dplyr::select(target, value, pair_id)
compare_df <- dplyr::left_join(x = glmeiv_df, y = thresh_df,
                               by = c("pair_id", "target"), suffix = c("_glmeiv", "_thresholding")) %>%
  na.omit() %>% dplyr::mutate(target = factor(labels = c("Estimate", "Lower CI", "Upper CI"),
                                              levels = c("estimate", "confint_lower", "confint_upper"),
                                              x = target))

p1 <- ggplot(data = compare_df, mapping = aes(x = value_glmeiv, y = value_thresholding)) + 
  geom_point(alpha = 0.7) + facet_grid(.~target) + xlab("GLM-EIV") + ylab("Thresholding") + theme_bw() +
  geom_abline(slope = 1, intercept = 0, col = "blue") + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = paste0(fig5_dir, "comparison.pdf"), plot = p1, device = "pdf", scale = 1, width = 7, height = 3)
