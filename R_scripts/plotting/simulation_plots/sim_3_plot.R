library(dplyr)
library(ggplot2)
library(simulatr)

# load the results and specifier objects
sim_result_dir <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/simulations/results")
sim_res_3 <- readRDS(paste0(sim_result_dir, "/sim_res_3.rds"))[["metrics"]]
source(paste0(.get_config_path("LOCAL_CODE_DIR"), "glmeiv-manuscript/R_scripts/plotting/simulation_plots/shared_plotting_code.R"))
to_plot <- sim_res_3 |>
  mutate(fam_str = ifelse(run_mrna_unknown_theta_precomputation, "nb_theta_unknown", "nb_theta_known")) |>
  transform_metric_df()

p <- ggplot(data = to_plot, mapping = aes(x = theta,
                                          y = mean,
                                          ymin = mean - 2 * se,
                                          ymax = mean + 2 * se,
                                          col = Method)) +
  xlab(expression(theta)) + scale_color_manual(values = my_cols) +
  scale_x_continuous(trans = "log10") +
  facet_grid(metric_fct ~ fam_str, scales = "free") +
  geom_hline(data = dplyr::filter(to_plot, metric_fct == "Bias"),
             mapping = aes(yintercept = 0), colour = "black") +
  geom_hline(data = dplyr::filter(to_plot, metric_fct == "MSE"),
             mapping = aes(yintercept = 0), colour = "black") +
  geom_hline(data = dplyr::filter(to_plot, metric_fct == "Coverage"),
             mapping = aes(yintercept = 0.95), colour = "black") +
  geom_line() + geom_errorbar(width = 0.05) + geom_point() +
  theme_bw() + theme(legend.position = "bottom", panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) + ylab("") +
ggplot2::guides(color = ggplot2::guide_legend(
  nrow = 3, byrow = TRUE,
  keyheight = 0.01,
  default.unit = "inch",
  override.aes = list(size = 1.25)))

ggsave(filename = "~/Desktop/p3.pdf", plot = p, device = "pdf", scale = 1, width = 4, height = 6)
