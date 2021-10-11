# Figure 3 displays the results of the "warmup" simulation -- no covariates, no offsets, no biological interpretation.

# set figure 3 dir
fig_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "glmeiv-manuscript/figures/")
fig3_dir <- paste0(fig_dir, "fig3/")
if (!dir.exists(fig3_dir)) dir.create(fig3_dir, recursive = TRUE)

# load packages and set sim_dir location
sim_dir <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "private/simulations")
library(ggplot2)
library(magrittr)

######################
# Distribution 1: Gaus
######################
sim_spec <- readRDS(paste0(sim_dir, "/sim_spec_4.rds")) # simulatr specifier object
sim_res <- readRDS(paste0(sim_dir, "/raw_result_4.rds")) # raw results; note: PSC failed for certain settings in

# possibly subset the sim_res according to some criteria
id_classifications <- glmeiv::obtain_valid_ids(sim_res = sim_res, pi_upper = 0.4)
valid_ids <- id_classifications$valid_ids
sim_res_sub <- dplyr::filter(sim_res, id %in% valid_ids)

summarized_results_gaus <- simulatr::summarize_results(sim_spec = sim_spec, sim_res = sim_res_sub,
                                                       metrics = c("bias", "mse", "coverage"),
                                                       parameters = c("m_perturbation"),
                                                       threshold = 0.1) %>% tibble::as_tibble() %>% dplyr::filter(arm_arm_pi_big)
pi <- summarized_results_gaus$pi[1]  
m_pert <- sim_spec@fixed_parameters$m_perturbation
my_funct <- function(x) glmeiv::get_tresholding_estimator_bias(m_perturbation = m_pert, g_perturbation = x, pi = pi)

density_dfs_and_thresholds <- glmeiv::get_theoretical_densities_and_thresholds(sim_spec = sim_spec,
                                                                               xgrid = seq(-10, 10, 0.1))
idx <- 1
glmeiv::plot_mixture(density_df = density_dfs_and_thresholds$g_dfs[[idx]],
                     thresh = density_dfs_and_thresholds$g_thresholds[idx],
                     x_max = 5, x_min = 0, points = FALSE)

######################
# Distribution 2: Pois
######################
sim_spec <- readRDS(paste0(sim_dir, "/sim_spec_1.rds")) # simulatr specifier object
sim_res <- readRDS(paste0(sim_dir, "/raw_result_1.rds")) # raw results

id_classifications <- glmeiv::obtain_valid_ids(sim_res)
valid_ids <- id_classifications$valid_ids
sim_res_sub <- dplyr::filter(sim_res, id %in% valid_ids)

# compute summary statistics
summarized_results_pois <- simulatr::summarize_results(sim_spec = sim_spec, sim_res = sim_res_sub,
                                                  metrics = c("coverage", "bias", "mse"),
                                                  parameters = c("m_perturbation")) %>% dplyr::filter(arm_g_perturbation)

####################
# Distribution 3: NB
####################
sim_spec <- readRDS(paste0(sim_dir, "/sim_spec_5.rds")) # simulatr specifier object
sim_res <- readRDS(paste0(sim_dir, "/raw_result_5.rds")) # raw results; note: PSC failed for certain settings in



id_classifications <- glmeiv::obtain_valid_ids(sim_res = sim_res, pi_upper = 0.4)
valid_ids <- id_classifications$valid_ids
sim_res_sub <- dplyr::filter(sim_res, id %in% valid_ids)

summarized_results_nb <- simulatr::summarize_results(sim_spec = sim_spec, sim_res = sim_res_sub,
                                                     metrics = c("bias", "mse", "coverage"),
                                                     parameters = c("m_perturbation"),
                                                     threshold = 0.1) %>% tibble::as_tibble() %>% dplyr::filter(arm_g_perturbation)

#############################
# Plot all g_perturbation arm
#############################
to_plot_all <- rbind(summarized_results_gaus %>%
                       dplyr::select(g_perturbation, method, metric, value, lower_mc_ci, upper_mc_ci) %>%
                       dplyr::mutate(distribution = "Gaussian"),
                     summarized_results_pois %>% dplyr::select(g_perturbation, method, metric, value, lower_mc_ci, upper_mc_ci) %>%
                       dplyr::mutate(distribution = "Poisson"),
                     summarized_results_nb %>% dplyr::select(g_perturbation, method, metric, value, lower_mc_ci, upper_mc_ci) %>%
                       dplyr::mutate(distribution = "Negative Binomial")) %>%
  dplyr::mutate(metric_fct = factor(metric, levels = c("bias", "mse", "coverage"),
                                    labels = c("Bias", "MSE", "Coverage")),
                distribution = factor(distribution, levels = c("Gaussian", "Poisson", "Negative Binomial")),
                Method = factor(method, levels = c("em", "thresholding"), labels = c("GLM-EIV", "Thresholding")))

p1 <- ggplot(to_plot_all, aes(x = g_perturbation, value, col = Method)) +
  facet_wrap(distribution ~ metric_fct, scales = "free") + ylab("") + xlab(expression(beta[g])) +
  scale_color_manual(values = my_cols) +
  geom_hline(data = dplyr::filter(to_plot_all, metric_fct == "Bias"), mapping = aes(yintercept = 0), colour = "black") +
  geom_hline(data = dplyr::filter(to_plot_all, metric_fct == "MSE"), mapping = aes(yintercept = 0), colour = "black") +
  geom_hline(data = dplyr::filter(to_plot_all, metric_fct == "Coverage"), mapping = aes(yintercept = 0.95), colour = "black") +
  geom_line() +
  geom_function(data = dplyr::filter(to_plot_all, metric_fct == "Bias", distribution == "Gaussian"), fun = my_funct, col = "black", lwd = 0.8) +
  geom_vline(data = dplyr::filter(to_plot_all, distribution == "Poisson"), mapping = aes(xintercept = 1.5), col = "gray60") + 
  geom_vline(data = dplyr::filter(to_plot_all, distribution == "Poisson"), mapping = aes(xintercept = 2.5), col = "gray60") + geom_point() +
  theme_bw() + theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename = paste0(fig3_dir, "simulation_warmup_arm_g.pdf"), plot = p1, device = "pdf", scale = 1, width = 7.5, height = 6.5)
