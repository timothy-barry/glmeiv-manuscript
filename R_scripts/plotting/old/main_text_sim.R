library(dplyr)
library(ggplot2)
library(simulatr)

my_cols <- c("firebrick3", "dodgerblue3", "orchid4")
fig_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "glmeiv-manuscript/figures/main_text_sim")
if (!dir.exists(fig_dir)) dir.create(fig_dir)

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
                                             metrics = c("bias", "mse", "coverage", "ci_width", "count", "time"),
                                             parameters = "m_perturbation") %>% mutate(distribution = "Poisson")
# remove times with 0 counts
thresh_grid_row_exclude <- summarized_pois_results %>% filter(metric == "count", value == 0) %>% dplyr::pull(grid_row_id)
summarized_pois_results_trimmed <- summarized_pois_results %>%
   filter(!(grid_id %in% thresh_grid_row_exclude & method == "thresholding"))

#######
# 2. NB
#######
nb_res <- readRDS(paste0(sim_result_dir, "/raw_result_1.rds")) %>% tibble()
nb_res$target[nb_res$target == "confint_higher"] <- "confint_upper"
nb_spec <- readRDS(paste0(sim_spec_dir, "/sim_spec_1.rds"))

# summarize the results
summarized_nb_results <- summarize_results(sim_spec = nb_spec, sim_res = nb_res,
                                           metrics = c("bias", "mse", "ci_width", "coverage", "count", "time"),
                                           parameters = "m_perturbation")
thresh_grid_row_exclude <- summarized_nb_results %>% filter(metric == "count", value == 0) %>% dplyr::pull(grid_row_id)
summarized_nb_results_trimmed <- summarized_nb_results %>% 
  filter(!(grid_id %in% thresh_grid_row_exclude & method == "thresholding")) %>%
  mutate(distribution = ifelse(fam_str == "nb_theta_known", "NB (theta known)", "NB (theta est.)")) %>%
  select(-g_fam, -m_fam, -run_unknown_theta_precomputation, -fam_str)

##################
# 3. Make the plot
##################
to_plot_all <- rbind(summarized_pois_results_trimmed, summarized_nb_results_trimmed) %>% filter(metric != "count") %>%
  mutate(metric_fct = factor(metric, levels = c("bias", "mse", "coverage", "ci_width", "time"),
                             labels = c("Bias", "MSE", "Coverage", "CI width", "Time (s)")),
         Method = factor(method, levels = c("glmeiv_slow", "glmeiv_fast", "thresholding"),
                         labels = c("GLM-EIV", "GLM-EIV (accelerated)", "Thresholding")),
         distribution = factor(distribution, levels = c("Poisson", "NB (theta known)", "NB (theta est.)"),
                               labels = c("Poisson", "NB (theta known)", "NB (theta est.)"))) %>%
  arrange(Method) %>% mutate(exp_g_perturbation = exp(g_perturbation))

p <- ggplot(data = to_plot_all, mapping = aes(x = exp_g_perturbation, y = value, col = Method)) + 
  facet_grid(metric_fct ~ distribution, scales = "free") + xlab(expression(exp(beta[g]))) + scale_color_manual(values = my_cols) +
  geom_hline(data = dplyr::filter(to_plot_all, metric_fct == "Bias"), mapping = aes(yintercept = 0), colour = "black") +
  geom_hline(data = dplyr::filter(to_plot_all, metric_fct == "MSE"), mapping = aes(yintercept = 0), colour = "black") +
  geom_hline(data = dplyr::filter(to_plot_all, metric_fct == "Coverage"), mapping = aes(yintercept = 0.95), colour = "black") +
  geom_line() + geom_point() + 
  theme_bw() + theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab("")

# save
fp <- paste0(fig_dir, "/plot.pdf")
ggsave(filename = fp, plot = p, device = "pdf", scale = 1.4, width = 4.25, height = 5.25)
