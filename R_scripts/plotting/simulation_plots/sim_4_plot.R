library(dplyr)
library(ggplot2)
library(simulatr)

# load the results and specifier objects
sim_result_dir <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/simulations/results")
sim_res_4 <- readRDS(paste0(sim_result_dir, "/sim_res_4.rds"))[["metrics"]]
source(paste0(.get_config_path("LOCAL_CODE_DIR"), "glmeiv-manuscript/R_scripts/plotting/simulation_plots/shared_plotting_code.R"))

to_plot <- sim_res_4 |> transform_metric_df() |> mutate(exp_m_perturbation = exp(m_perturbation))

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

ggsave(filename = "~/Desktop/p4.pdf", plot = p, device = "pdf", scale = 1, width = 5, height = 6)
