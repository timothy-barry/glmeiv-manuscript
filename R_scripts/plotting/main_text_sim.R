library(dplyr)
library(ggplot2)
load_all("~/research_code/simulatr/")

my_cols <- c("firebrick3", "dodgerblue3", "orchid4")
fig_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "glmeiv-manuscript/figures/main_text_sim")
if (!dir.exists(fig_dir)) dir.create(fig_dir)
my_theme <- theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# load the results and specifier objects
sim_result_dir <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/simulations/results")
sim_spec_dir <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/simulations/spec_objects")

############
# 1. Poisson
############
pois_res <- readRDS(paste0(sim_result_dir, "/raw_result_0.rds")) %>% tibble()
pois_res$target[pois_res$target == "confint_higher"] <- "confint_upper"

pois_spec <- readRDS(paste0(sim_spec_dir, "/sim_spec_0.rds"))

# summarize the results
summarized_pois_results <- summarize_results(sim_spec = pois_spec, sim_res = pois_res,
                                             metrics = c("bias", "mse", "coverage", "count", "time"),
                                             parameters = "m_perturbation") %>% mutate(distribution = "Poisson")

##################
# 3. Make the plot
##################

to_plot_all <- summarized_pois_results %>% filter(metric != "count") %>%
  mutate(metric_fct = factor(metric, levels = c("bias", "mse", "coverage", "time"),
                             labels = c("Bias", "MSE", "Coverage", "Time (s)")),
         Method = factor(method, levels = c("glmeiv_slow", "glmeiv_fast", "thresholding"),
                         labels = c("GLM-EIV", "GLM-EIV (accelerated)", "Thresholding"))) %>%
  arrange(Method) %>% mutate(exp_g_perturbation = exp(g_perturbation))

ggplot(data = to_plot_all, mapping = aes(x = exp_g_perturbation, y = value, col = Method)) + 
  facet_grid(metric_fct ~ distribution, scales = "free_y") + xlab(expression(exp(beta[g]))) + scale_color_manual(values = my_cols) +
  geom_hline(data = dplyr::filter(to_plot_all, metric_fct == "Bias"), mapping = aes(yintercept = 0), colour = "black") +
  geom_hline(data = dplyr::filter(to_plot_all, metric_fct == "MSE"), mapping = aes(yintercept = 0), colour = "black") +
  geom_hline(data = dplyr::filter(to_plot_all, metric_fct == "Coverage"), mapping = aes(yintercept = 0.95), colour = "black") +
  geom_line() + geom_point() + 
  theme_bw() + theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab("")
