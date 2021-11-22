library(dplyr)
library(ggplot2)
load_all("~/research_code/simulatr/")

my_cols <- c("dodgerblue3", "orchid4")
fig_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "glmeiv-manuscript/figures/supplement_sim")
if (!dir.exists(fig_dir)) dir.create(fig_dir)

# load the results and specifier objects
sim_result_dir <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/simulations/results")
sim_spec_dir <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/simulations/spec_objects")

#############
# 1. Gaussian
#############
gaus_res <- readRDS(paste0(sim_result_dir, "/raw_result_2.rds")) %>% tibble()
gaus_res$target[gaus_res$target == "confint_higher"] <- "confint_upper"
gaus_spec <- readRDS(paste0(sim_spec_dir, "/sim_spec_2.rds"))

# summarize results
summarized_gaus_results <- summarize_results(sim_spec = gaus_spec, sim_res = gaus_res,
                                             metrics = c("bias", "mse", "coverage", "ci_width"),
                                             parameters = "m_perturbation")

to_plot_all <- summarized_gaus_results %>% filter(metric != "count") %>%
  mutate(metric_fct = factor(metric, levels = c("bias", "mse", "coverage", "ci_width", "time"),
                             labels = c("Bias", "MSE", "Coverage", "CI width", "Time (s)")),
         Method = factor(method, levels = c("glmeiv_slow", "glmeiv_fast", "thresholding"),
                         labels = c("GLM-EIV", "GLM-EIV (accelerated)", "Thresholding"))) %>%
  arrange(Method)

p <- ggplot(data = to_plot_all, mapping = aes(x = g_perturbation, y = value, col = Method)) + 
  facet_wrap(. ~ metric_fct, scales = "free_y", nrow = 1) + xlab(expression(beta[g])) + scale_color_manual(values = my_cols) +
  geom_hline(data = dplyr::filter(to_plot_all, metric_fct == "Bias"), mapping = aes(yintercept = 0), colour = "black") +
  geom_hline(data = dplyr::filter(to_plot_all, metric_fct == "MSE"), mapping = aes(yintercept = 0), colour = "black") +
  geom_hline(data = dplyr::filter(to_plot_all, metric_fct == "Coverage"), mapping = aes(yintercept = 0.95), colour = "black") +
  geom_line() + geom_point() + 
  theme_bw() + theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab("")

fp <- paste0(fig_dir, "/gaussian.pdf")
ggsave(filename = fp, plot = p, device = "pdf", scale = 1.1, width = 7, height = 2.5)
