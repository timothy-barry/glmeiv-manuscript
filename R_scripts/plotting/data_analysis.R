library(dplyr)
library(ggplot2)
library(ggpmisc)
library(tidyr)

my_theme <- theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_rect("white"))


my_cols <- c("firebrick3", "dodgerblue3", "orchid4")
fig_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "glmeiv-manuscript/figures/data_analysis")
if (!dir.exists(fig_dir)) dir.create(fig_dir)
result_dir_gasp <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/gasperini/results")

########################
# 1. Resampling analysis
########################
n_pairs <- resampling_df$pair_id %>% unique() %>% length()
resampling_df <- readRDS(paste0(result_dir_gasp, "/resampling_result.rds"))
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
p1 <- f(resampling_df, pair_ids[37], x_max = 0.4) # good choices: 37

# an "aggregate" version of the above analysis
aggregate_df <- resampling_df %>% filter(contam_level >= 0, parameter == "m_perturbation", target %in% c("estimate", "confint_lower", "confint_upper")) %>%
  select(-parameter) %>% tidyr::pivot_wider(names_from = "target", values_from = "value") %>% group_by(method, pair_id, contam_level) %>%
  summarize(m = mean(estimate), m_lower_ci = mean(confint_lower, na.rm = TRUE), m_upper_ci = mean(confint_upper, na.rm = TRUE)) %>% ungroup() %>% group_by(method, pair_id) %>%
  group_modify(.f = function(tbl, key) {
    baseline_est <- tbl$m[tbl$contam_level == 0]
    p_change <- abs(tbl$m - baseline_est)/baseline_est
    ci_cover <- (tbl$m_lower_ci < baseline_est) & (tbl$m_upper_ci > baseline_est)
    mutate(tbl, p_change = p_change, ci_cover = ci_cover)
  }) %>% ungroup() %>% group_by(contam_level, method) %>% summarize(median_p_change = median(p_change),
                                                                    lower_median_ci = sort(p_change)[qbinom(.025, length(p_change), 0.5)],
                                                                    upper_median_ci = sort(p_change)[qbinom(.975, length(p_change), 0.5)],
                                                                    mean_coverage = mean(ci_cover),
                                                                    lower_cover_ci = mean_coverage - 1.96 * sqrt((1/n_pairs) * mean_coverage * (1 - mean_coverage)),
                                                                    upper_cover_ci = mean_coverage + 1.96 * sqrt((1/n_pairs) * mean_coverage * (1 - mean_coverage))) %>%
  mutate(method = factor(method, c("glmeiv", "thresholding"), c("GLM-EIV", "Thresholding")))

p2 <- ggplot(data = aggregate_df %>% filter(contam_level <= 0.4),
             mapping = aes(x = contam_level, y = median_p_change, col = method)) + 
  geom_hline(yintercept = 0, col = "black", lwd = 0.5) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.21, 0.82), legend.title=element_blank()) +
  scale_color_manual(values = c(my_cols[1], my_cols[3])) +
  xlab("Excess background contamination") + scale_x_continuous(expand = expansion(mult = c(0,0.01))) + geom_ribbon(aes(ymin = lower_median_ci, ymax = upper_median_ci, fill = method), alpha = 0.5) + geom_line(lwd = 0.85) +
  ylab("Median REC across pairs") + scale_fill_manual(values = c(my_cols[1], my_cols[3]))

p3 <- ggplot(data = aggregate_df %>% filter(contam_level <= 0.4),
             mapping = aes(x = contam_level, y = mean_coverage, col = method)) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "blank", legend.title=element_blank()) +
  scale_color_manual(values = c(my_cols[1], my_cols[3])) +
  xlab("Excess background contamination") + scale_x_continuous(expand = expansion(mult = c(0, 0.01))) + geom_ribbon(aes(ymin = lower_cover_ci, ymax = upper_cover_ci, fill = method), alpha = 0.5) + geom_line(lwd = 0.85) +
  ylab("CI stability across pairs") + scale_fill_manual(values = c(my_cols[1], my_cols[3]))

# difference between resampled 0 and fitted 0.
raw_est_vs_fitted_model_est <- resampling_df %>% filter(contam_level %in% c(-1, 0), parameter == "m_perturbation", target == "estimate") %>%
  select(-parameter, -target) %>% group_by(method, pair_id, contam_level) %>%
  dplyr::summarize(m_est = mean(value)) %>% group_modify(.f = function(tbl, key) {
    raw_est <- tbl$m_est[tbl$contam_level < 0]
    resampled_est <- tbl$m_est[tbl$contam_level == 0]
    delta <- abs(raw_est - resampled_est)/raw_est
    tbl %>% mutate(delta = delta)
  }) %>% ungroup() %>% group_by(method) %>% summarize(median_delta = median(delta))
raw_est_vs_fitted_model_est

####################
# Raw data analyses
####################
# Define functions
# 1. join data frames for comparison
join_results <- function(glmeiv_res, thresh_res, site_type) {
  left_join(x = filter(.data = glmeiv_res, target == "estimate", parameter == "m_perturbation", type == !!site_type) %>% select(pair_id, value),
                           y = filter(thresh_res, target == "estimate", parameter == "m_perturbation", type == !!site_type) %>% select(pair_id, value),
                           by = "pair_id", suffix = c("_glmeiv", "_thresh"))
}

# 2. Compute CI coverage rate (assuming ground truth of 1) and width
get_ci_coverage_rate <- function(df, site_type = NULL) {
  df %>% dplyr::filter(parameter == "m_perturbation",
                       target %in% c("confint_lower", "confint_upper", "estimate"),
                       type %in% if (is.null(site_type)) as.character(unique(df$type)) else site_type) %>%
    dplyr::select(target, value, pair_id, type) %>%
    tidyr::pivot_wider(id_cols = c("pair_id", "type"),
                       names_from = "target", values_from = "value") %>%
    dplyr::mutate(covered = confint_upper > 1 & confint_lower < 1,
                  width = confint_upper - confint_lower) %>%
    dplyr::summarize(coverage_rate = mean(covered, na.rm = TRUE) * 100,
                     mean_width = mean(width, na.rm = TRUE)) 
}

####################
# Gasperini analysis
####################
#  1. load the results
glmeiv_res <- paste0(result_dir_gasp, "/result_glmeiv.rds") %>% readRDS() %>% mutate(type = site_type)
thresh_res <- paste0(result_dir_gasp, "/result_thresholding.rds") %>% readRDS() %>% mutate(type = site_type)

# 2. Ensure we are examining the same pairs (Gasp missing ~100 due to probable node failure)
glmeiv_pairs <- glmeiv_res %>% filter(parameter == "m_perturbation", target == "estimate") %>% pull(pair_id) %>% as.character() %>% sort()
thresh_pairs <- thresh_res %>% filter(parameter == "m_perturbation", target == "estimate") %>% pull(pair_id) %>% as.character() %>% sort()

glmeiv_res <- glmeiv_res %>% filter(pair_id %in% ok_glmeiv_pairs)
thresh_res <- thresh_res %>% filter(pair_id %in% ok_glmeiv_pairs)

# 4. get CI coverage for NTCs
get_ci_coverage_rate(glmeiv_res, "NTC")
get_ci_coverage_rate(thresh_res, "NTC")

# 5. Plot results on NTC pairs
comparison_df <- join_results(glmeiv_res, thresh_res, "NTC")
p3 <- ggplot(data = comparison_df, mapping = aes(x = value_glmeiv, y = value_thresh)) + geom_point(alpha = 0.5) + 
  geom_abline(slope = 1, intercept = 0, col = my_cols[2], lwd = 1) + xlab("GLM-EIV estimate") + ylab("Thresholding estimate") + my_theme

# Plot positive control pairs
comparison_df <- join_results(glmeiv_res, thresh_res, "selfTSS")
p4 <- ggplot(data = comparison_df, mapping = aes(x = value_glmeiv, y = value_thresh)) + geom_point(alpha = 0.5) + 
  geom_abline(slope = 1, intercept = 0, col = my_cols[2], lwd = 1) + xlab("GLM-EIV estimate") + ylab("Thresholding estimate") + my_theme

##############
# Xie analysis
##############
# 1. load the results
result_dir_xie <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/xie/results")
glmeiv_res <- paste0(result_dir_xie, "/glmeiv_result.rds") %>% readRDS()
thresh_res <- paste0(result_dir_xie, "/result_thresholding.rds") %>% readRDS()

# 2. confirm that the results are on the same pairs
glmeiv_pairs <- glmeiv_res %>% filter(parameter == "m_perturbation", target == "estimate") %>% pull(pair_id) %>% as.character() %>% sort()
thresh_pairs <- thresh_res %>% filter(parameter == "m_perturbation", target == "estimate") %>% pull(pair_id) %>% as.character() %>% sort()
identical(glmeiv_pairs, thresh_pairs)

# 3. Comipute mean g_perturbation
glmeiv_res %>% filter(parameter == "g_perturbation") %>% pivot_wider(id_cols = "pair_id", names_from = "target") %>%
  summarize(m_est = mean(estimate), m_confint_lower = mean(confint_lower, na.rm = TRUE), m_confint_upper = mean(confint_upper, na.rm = TRUE)) # g perturbation fairly large

# 4. get the CI coverage table
ci_info <- rbind(get_ci_coverage_rate(glmeiv_res, "neg_control"),
      get_ci_coverage_rate(thresh_res, "neg_control")) %>%
  mutate(Method = c("GLM-EIV", "Thresh."), coverage_rate = round(coverage_rate, 1), mean_width = round(mean_width, 3)) %>%
  select(Method, "Cov. rate" = coverage_rate, "Width" = mean_width)

# 5. Plot results
comparison_df <- join_results(glmeiv_res, thresh_res, "neg_control") %>% mutate(glmeiv_ok = value_glmeiv < 1.25 & value_glmeiv > 0.75)
comparison_df_sub <- comparison_df %>% filter(glmeiv_ok)
p3 <- ggplot(data = comparison_df_sub, mapping = aes(x = value_glmeiv, y = value_thresh)) + geom_point(alpha = 0.5) + 
  geom_abline(slope = 1, intercept = 0, col = my_cols[2], lwd = 1) + xlab("GLM-EIV estimate") + ylab("Thresholding estimate") + my_theme +
  annotate(geom = "table", x = min(comparison_df_sub$value_thresh), y = max(comparison_df_sub$value_thresh), label = list(ci_info), vjust = 1, hjust = 0, fill = "white")

# Examine bad pairs
bad_pairs <- glmeiv_res %>% filter(parameter == "m_perturbation", target == "estimate", type == "neg_control", value > 3.5)
glmeiv_res %>% filter(pair_id %in% bad_pairs, target  == "estimate" | parameter == "meta" ) %>% print(n = 20)
