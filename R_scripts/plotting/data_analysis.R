library(dplyr)
library(ggplot2)

my_cols <- c("firebrick3", "dodgerblue3", "orchid4")
fig_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "glmeiv-manuscript/figures/data_analysis")
if (!dir.exists(fig_dir)) dir.create(fig_dir)
result_dir <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/gasperini/results")

########################
# 1. Resampling analysis
########################
resampling_df <- readRDS(paste0(result_dir, "/resampling_result.rds"))
# a function to plot the result given a given pair
f <- function(resampling_df, pair_id, x_max = 0.5) {
  df <- dplyr::filter(resampling_df, pair_id == !!pair_id, contam_level <= x_max) %>%
    dplyr::filter(parameter == "m_perturbation", target %in% c("estimate", "confint_lower", "confint_upper")) %>%
    dplyr::select(-parameter, -pair_id) %>%
    tidyr::pivot_wider(data = ., names_from = "target", values_from = "value") %>%
    dplyr::mutate(method = factor(method, c("glmeiv", "thresholding"), c("GLM-EIV", "Thresholding")))
  df_synth <- df %>% dplyr::filter(contam_level >= 0) %>% dplyr::group_by(method, contam_level) %>%
    dplyr::summarize(m_est = mean(estimate, na.rm = TRUE), m_ci_upper = mean(confint_upper, na.rm = TRUE), m_ci_lower = mean(confint_lower, na.rm = TRUE))
  df_gt <- df %>% dplyr::filter(contam_level == -1)
  p <- ggplot(data = df_synth, mapping = aes(x = contam_level, y = m_est)) + facet_grid(.~method) +
    geom_hline(data = df_gt, mapping = aes(yintercept = estimate), color = my_cols[2], lwd = 0.75) +
    theme_bw() + xlab("Excess background contamination") + 
    scale_color_manual(values=c("grey70", "grey70", "black")) + ylab("Estiamte") + geom_ribbon(aes(ymin = m_ci_lower, ymax = m_ci_upper), fill = "grey70", alpha = 0.5) + geom_line(lwd = 0.75) +
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing.x = unit(4, "mm"), panel.border = element_blank(), axis.line = element_line(colour = "black")) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))
  return(p)
}
pair_ids <- as.character(unique(resampling_df$pair_id))
f(resampling_df, pair_ids[10], x_max = 0.4) %>% plot()

# an "aggregate" version of the above analysis
aggregate_df <- resampling_df %>% filter(contam_level >= 0, parameter == "m_perturbation", target == "estimate") %>%
  select(-parameter, -target) %>% group_by(method, pair_id, contam_level) %>%
  summarize(m = mean(value)) %>% ungroup() %>% group_by(method, pair_id) %>%
  group_modify(.f = function(tbl, key) {
    baseline_est <- tbl$m[tbl$contam_level == 0]
    p_change <- abs(tbl$m - baseline_est)/baseline_est
    mutate(tbl, p_change = p_change)
  }) %>% ungroup() %>% group_by(contam_level, method) %>% summarize(median_p_change = median(p_change),
                                                                    mean_p_change = mean(p_change),
                                                                    lower_median_ci = sort(p_change)[qbinom(.025, length(p_change), 0.5)],
                                                                    upper_median_ci = sort(p_change)[qbinom(.975, length(p_change), 0.5)]) %>%
  mutate(method = factor(method, c("glmeiv", "thresholding"), c("GLM-EIV", "Thresholding")))

p2 <- ggplot(data = aggregate_df %>% filter(contam_level <= 0.4),
             mapping = aes(x = contam_level, y = median_p_change, col = method)) + 
  geom_hline(yintercept = 0, col = "black", lwd = 0.5) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.21, 0.82), legend.title=element_blank()) +
  scale_color_manual(values = c(my_cols[1], my_cols[3])) +
  xlab("Excess background contamination") + scale_x_continuous(expand = expansion(mult = c(0,0.01))) + geom_ribbon(aes(ymin = lower_median_ci, ymax = upper_median_ci, fill = method), alpha = 0.5) + geom_line(lwd = 0.85) +
  ylab("Median REC across pairs") + scale_fill_manual(values = c(my_cols[1], my_cols[3]))

# difference between resampled 0 and fitted 0.