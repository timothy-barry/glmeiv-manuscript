library(magrittr)
library(ggplot2)
library(cowplot)

my_cols <- c("firebrick3", "dodgerblue3", "orchid4", "lightskyblue3")
fig_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "glmeiv-manuscript/figures/thresholding_empirical")
if (!dir.exists(fig_dir)) dir.create(fig_dir)
my_theme <- theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# load the results
gasp_result_dir <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/gasperini/results/")
thresh_res <- readRDS(paste0(gasp_result_dir, "thresholding_result_pc.rds"))

# study how choice of threshold influences results on Gasperini self TSS pairs
thresh_res_pc <- thresh_res %>% dplyr::filter(site_type == "selfTSS", parameter == "m_perturbation")

# Question 1. How do thresh = 1 and thresh = 5 compare on estimation?
thresh_res_pc_wide <- thresh_res_pc %>% dplyr::select(target, value, threshold, pair_id) %>%
  tidyr::pivot_wider(id_cols = c("pair_id", "target"),
                     names_from = "threshold",
                     values_from = "value", names_prefix = "threshold_")

to_plot <- thresh_res_pc_wide %>% dplyr::filter(target == "estimate") %>% dplyr::select(threshold_1, threshold_5)
p1 <- ggplot(data = to_plot, mapping = aes(x = threshold_5, y = threshold_1)) + geom_point(alpha = 0.7, col = my_cols[2]) +
  xlab("Thresh = 5 (est)") + ylab("Thresh = 1 (est)") + geom_abline(slope = 1, intercept = 0) + theme_bw(base_size = 10) + my_theme

# threshold = 1 ests are bigger than (i.e., closer to one than) threshold = 5; this is attenuation bias due to the fact that excessively many unperturbed cells have been classified as perturbed.

# a) est. vs. est.
to_plot <- thresh_res_pc_wide %>% dplyr::filter(target == "estimate") %>% dplyr::select(threshold_5, threshold_20)
p2 <- ggplot(data = to_plot, mapping = aes(x = threshold_5, y = threshold_20)) + geom_point(alpha = 0.9, col = my_cols[2]) +
  xlab("Thresh = 5 (est)") + ylab("Thresh = 20 (est)") + geom_abline(slope = 1, intercept = 0, lwd = 0.75) + theme_bw(base_size = 10) + my_theme

# b) p-value vs. p-value
to_plot <- thresh_res_pc_wide %>% dplyr::filter(target == "p_value") %>% dplyr::select(threshold_5, threshold_20)
p3 <- ggplot(data = to_plot, mapping = aes(x = -log(threshold_5, base = 10), y = -log(threshold_20, base = 10))) + geom_point(alpha = 0.7, col = my_cols[3]) +
  xlab("Thresh = 5 (-log p)") + ylab("Thresh = 20 (-log p)") + geom_abline(slope = 1, intercept = 0) + theme_bw(base_size = 10) + my_theme

# c) CI width for threshold = 5 vs threshold = 20
to_plot <- thresh_res_pc_wide %>% dplyr::filter(target == c("confint_lower", "confint_upper")) %>% dplyr::select(-threshold_1) %>%
  tidyr::pivot_longer(cols = c("threshold_5", "threshold_20"), names_to = "threshold") %>% dplyr::group_by(threshold) %>%
  dplyr::group_modify(function(tbl, key) {
    tidyr::pivot_wider(data = tbl, id_cols = "pair_id", names_from = "target", values_from = "value")
  }) %>% dplyr::mutate(ci_width = log(confint_upper) - log(confint_lower)) %>%
  dplyr::mutate(threshold = factor(x = threshold, levels = c("threshold_5", "threshold_20"), labels = c("5", "20")))

mean_ci_width <- to_plot %>% dplyr::summarize(m_width = mean(ci_width))
(mean_ci_width %>% dplyr::filter(threshold == "20") %>% dplyr::pull(m_width))/(mean_ci_width %>% dplyr::filter(threshold == "5") %>% dplyr::pull(m_width))
# Threshold = 20 CI width is 1.5 times the threshold = 5 CI width.
p4 <- ggplot(data = to_plot, mapping = aes(x = threshold, y = ci_width)) + geom_violin(fill = my_cols[1], alpha = 0.7, draw_quantiles = 0.5) +
  theme_bw(base_size = 10) + xlab("Threshold") + ylab("CI width") + geom_hline(yintercept = 0) + my_theme + scale_y_continuous(expand = c(NA, 0))

# e and f) gRNA count histograms
plotting_objects <- paste0(.get_config_path("LOCAL_CODE_DIR"), "glmeiv-manuscript/figures/plotting_objects.R")
source(plotting_objects)
set.seed(2)
gRNA_ids_gasp <- sample(x = get_feature_ids(gRNA_odm_gasp), size = 1, replace = FALSE)
gRNA_ids_xie <- sample(x = get_feature_ids(gRNA_odm_xie), size = 1, replace = FALSE)

xie_gRNA <- gRNA_odm_xie[[gRNA_ids_xie,]] %>% as.numeric()
gasp_gRNA <- gRNA_odm_gasp[[gRNA_ids_gasp,]] %>% as.numeric()

# plot histograms
p5 <- ggplot(data = xie_gRNA %>% tibble::tibble(count = .) %>%
               dplyr::filter(count >= 1, count <= 40),
             mapping = aes(x = count)) +
  geom_histogram(binwidth = 2, col = "black", fill = my_cols[4], alpha = 0.7) +
  scale_y_continuous(trans='log10', expand = c(0, NA)) + xlab("gRNA count (Xie)") + ylab("") +
  theme_bw(base_size = 10) + my_theme
  
p6 <- ggplot(data = gasp_gRNA %>% tibble::tibble(count = .) %>%
               dplyr::filter(count >= 1, count <= 40),
             mapping = aes(x = count)) +
  geom_histogram(binwidth = 2, col = "black", fill = my_cols[4], alpha = 0.7) +
  scale_y_continuous(trans='log10', expand = c(0, NA)) + xlab("gRNA count (Gasperini)") +
  ylab("") + theme_bw(base_size = 10) + my_theme

# combine plots
p_out <- plot_grid(p1, p2, p3, p4, p6, p5, labels = c("a", "b", "c", "d", "e", "f"), ncol = 2, rel_heights = c(1,1,0.8), align = "vh")

ggsave(filename = paste0(fig_dir, "/plot.jpeg"), plot = p_out, device = "jpeg", scale = 1, width = 6, height = 6)
