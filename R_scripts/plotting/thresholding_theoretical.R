load_all("~/research_code/glmeiv/")
my_cols <- c("firebrick3", "dodgerblue4", "orchid4")
my_theme <- theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
my_lwd <- 0.65
library(magrittr)
library(ggplot2)
library(cowplot)

fig_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "glmeiv-manuscript/figures/thresholding_theoretical")
if (!dir.exists(fig_dir)) dir.create(fig_dir)

# a) estimate vs. threshold for different values of beta_g
x_grid <- seq(0, 7, length.out = 501)
beta_gs <- c(1, qnorm(0.75) * 2, 1.5)

to_plot_a <- lapply(X = beta_gs, function(curr_beta_g) {
  y <- get_tresholding_estimator_bias(m_perturbation = 1,
                                      g_perturbation = curr_beta_g,
                                      g_intercept = 0,
                                      pi = 0.5,
                                      c = x_grid,
                                      return_bias = TRUE)
  data.frame(Bias = y, Threshold = x_grid, beta_g = curr_beta_g)
}) %>% do.call(what = "rbind", args = .) %>%
  dplyr::mutate(beta_g_disp = paste0("beta[1]^g == ", round(beta_g, 2)) %>% factor())

# subplot a) bias in pi = 1/2 model
p_a <- ggplot(data = to_plot_a, mapping = aes(x = Threshold, y = Bias)) +
  geom_hline(yintercept = 0.5, col = "darkred", lwd = my_lwd) +
  facet_wrap(.~beta_g_disp, labeller = label_parsed, scales = "free_y") +
  theme_bw(base_size = 10) + scale_x_continuous(expand = c(NA, 0)) +
  geom_vline(data = dplyr::filter(to_plot_a, beta_g_disp == "beta[1]^g == 1"),
             mapping = aes(xintercept = beta_gs[1]/2), col = my_cols[2], lwd = my_lwd) +
  geom_vline(data = dplyr::filter(to_plot_a, beta_g_disp == "beta[1]^g == 1.35"),
             mapping = aes(xintercept = beta_gs[2]/2), col = my_cols[2], lwd = my_lwd) +
  geom_vline(data = dplyr::filter(to_plot_a, beta_g_disp == "beta[1]^g == 1.5"),
             mapping = aes(xintercept = beta_gs[3]/2), col = my_cols[2], lwd = my_lwd) +
  geom_line(lwd = my_lwd) + my_theme

# subplot b) bias as c -> infinity for pi in [0,1/2]
pi <- seq(0.0, 0.5, length.out = 11)
to_plot_b <- data.frame(pi = pi, Bias = pi)

p_b <- ggplot(data = to_plot_b, mapping = aes(x = pi, y = Bias)) + geom_line(lwd = my_lwd) + theme_bw(base_size = 10) +
  xlab(expression(pi)) + scale_x_continuous(expand = c(0.005, 0.005)) + scale_y_continuous(expand = c(0.005, 0.005)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# subplot c) bias-variance tradeoff
xgrid <- seq(1, 3.9, length.out = 1001)
g_beta <- 1
pi <- 0.1
m_beta <- 1
n <- 10000

bias_sq <- (1 - get_thresholding_estimator_est_no_int(c = xgrid, g_beta = g_beta, pi = pi, m_beta = m_beta))^2
var <- get_thresholding_estimator_var_no_int(c = xgrid, g_beta = g_beta, pi = pi, m_beta = m_beta, n = n)
to_plot_c <- data.frame(Metric = rep(c("Bias^2", "Variance", "MSE"), each = length(xgrid)),
                        value = c(bias_sq, var, bias_sq + var),
                        threshold = rep(xgrid, times = 3))

p_c <- ggplot(data = to_plot_c %>% dplyr::mutate(Metric = factor(Metric, levels = c("Bias^2", "Variance", "MSE"),
                                                                 labels = c("Bias^2", "Variance", "MSE"))),
              mapping = aes(x = threshold, y = value, col = Metric)) +
  geom_line(lwd = my_lwd) + theme_bw(base_size = 10) + xlab("Threshold") + ylab("") +
  scale_color_manual(values = c(my_cols[3], my_cols[1], my_cols[2]),
                     labels = c("MSE", expression(Bias^2), "Variance"),
                     breaks = c("MSE", expression(Bias^2), "Variance")) + scale_x_continuous(expand = c(0.005, 0.005)) + 
  scale_y_continuous(expand = c(0.005, 0.005)) + theme(legend.position = c(0.14, 0.34), legend.title = element_blank(),
                                                       panel.border = element_blank(), panel.grid.major = element_blank(),
                                                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# combine
# p_bottom <- plot_grid(p_b, p_c, align = "h", ncol = 2, rel_widths = c(1, 2), labels = c("b", "c"))
# p <- plot_grid(p_a, p_bottom, ncol = 1, labels = c("a", ""))

# save
ggsave(filename = paste0(fig_dir, "/plot.jpeg"), plot = p_c, device = "jpeg", scale = 1.0, width = 4.5, height = 2.5, dpi = 320)
ggsave(filename = paste0(fig_dir, "/app_plot.jpeg"), plot = p_a, device = "jpeg", scale = 1.5, width = 4.5, height = 1.5, dpi = 320)
