library(dplyr)
library(ggplot2)
library(simulatr)

my_cols <- c("firebrick3", "dodgerblue3", "orchid4")

# load the results and specifier objects
sim_result_dir <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/simulations/results")
sim_res_4 <- readRDS(paste0(sim_result_dir, "/sim_res_4.rds"))[["metrics"]]

to_plot <- sim_res_4 |>
  dplyr::filter(metric %in% c("bias", "ci_width", "coverage", "mse", "time")) |>
  mutate(metric_fct = factor(metric, levels = c("bias", "mse", "coverage", "ci_width", "time"),
                             labels = c("Bias", "MSE", "Coverage", "CI width", "Time (s)")),
         Method = factor(method, levels = c("glmeiv_slow", "glmeiv_fast", "thresholding"),
                         labels = c("GLM-EIV", "GLM-EIV (accelerated)", "Thresholding"))) |>
  arrange(Method) |> mutate(exp_m_perturbation = exp(m_perturbation))

p <- ggplot(data = to_plot, mapping = aes(x = exp_m_perturbation,
                                          y = mean,
                                          ymin = mean - 2 * se, 
                                          ymax = mean + 2 * se,
                                          col = Method)) + 
  xlab(expression(exp(beta[m]))) + scale_color_manual(values = my_cols) +
  facet_grid(metric_fct ~ fam_str, scales = "free") + 
  geom_hline(data = dplyr::filter(to_plot, metric_fct == "Bias"),
             mapping = aes(yintercept = 0), colour = "black") +
  geom_hline(data = dplyr::filter(to_plot, metric_fct == "MSE"),
             mapping = aes(yintercept = 0), colour = "black") +
  geom_hline(data = dplyr::filter(to_plot, metric_fct == "Coverage"),
             mapping = aes(yintercept = 0.95), colour = "black") +
  geom_line() + geom_errorbar(width = 0.05) + geom_point() + 
  theme_bw() + theme(legend.position = "bottom", panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) + ylab("")
