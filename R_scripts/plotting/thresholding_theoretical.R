load_all("~/research_code/glmeiv/")
my_cols <- c("firebrick3", "dodgerblue3")
library(magrittr)
library(ggplot2)
library(cowplot)

# a) estimate vs. threshold for different values of beta_g
x_grid <- seq(0, 7, length.out = 501)
beta_gs <- c(1, qnorm(0.75) * 2, 1.5)

to_plot_a <- lapply(X = beta_gs, function(curr_beta_g) {
  y <- get_tresholding_estimator_bias(m_perturbation = 1,
                                      g_perturbation = curr_beta_g,
                                      pi = 0.5,
                                      c = x_grid,
                                      return_bias = FALSE)
  data.frame(Estimate = y, Threshold = x_grid, beta_g = curr_beta_g)
}) %>% do.call(what = "rbind", args = .) %>%
  dplyr::mutate(beta_g_disp = paste0("beta[g] == ", round(beta_g, 2)) %>% factor())

# subplot a) bias in pi = 1/2 model
p_a <- ggplot(data = to_plot_a, mapping = aes(x = Threshold, y = Estimate)) +
  geom_hline(yintercept = 0.5, col = "darkred") +
  facet_wrap(.~beta_g_disp, labeller = label_parsed, scales = "free_y") +
  theme_cowplot(font_size = 11) + scale_x_continuous(expand = c(NA, 0)) +
  geom_vline(data = dplyr::filter(to_plot_a, beta_g_disp == "beta[g] == 1"),
             mapping = aes(xintercept = beta_gs[1]/2), col = my_cols[2]) +
  geom_vline(data = dplyr::filter(to_plot_a, beta_g_disp == "beta[g] == 1.35"),
             mapping = aes(xintercept = beta_gs[2]/2), col = my_cols[2]) +
  geom_vline(data = dplyr::filter(to_plot_a, beta_g_disp == "beta[g] == 1.5"),
             mapping = aes(xintercept = beta_gs[3]/2), col = my_cols[2]) +
  geom_line()

# subplot b) bias as c -> infinity for pi in [0,1/2]
pi <- seq(0.0, 0.5, length.out = 11)
Estimate <- 1 - pi
to_plot_b <- data.frame(pi = pi, Estimate = Estimate)

p_b <- ggplot(data = to_plot_b, mapping = aes(x = pi, y = Estimate)) + geom_line() +
  scale_x_continuous(trans = "reverse", expand = c(NA, 0), limits = c(0.5, -0.02)) + theme_cowplot(font_size = 11) +
  geom_segment(aes(x=0.0, xend=0.5, y=1, yend=1), lwd = 0.3) + scale_y_continuous(expand = c(0, NA), limits = c(0.5, 1.02)) +
  xlab(expression(pi))

# subplot c) bias-variance tradeoff
xgrid <- seq(1, 3.9, length.out = 1001)
g_beta <- 1
pi <- 0.25
m_beta <- 1
n <- 10000

bias_sq <- (1 - get_thresholding_estimator_est_no_int(c = xgrid, g_beta = g_beta, pi = pi, m_beta = m_beta))^2
var <- get_thresholding_estimator_var_no_int(c = xgrid, g_beta = g_beta, pi = pi, m_beta = m_beta, n = n)
to_plot_c <- data.frame(Metric = rep(c("Bias^2", "Variance", "MSE"), each = length(xgrid)),
                        value = c(bias_sq, var, bias_sq + var),
                        threshold = rep(xgrid, times = 3))
p_c <- ggplot(data = to_plot_c, mapping = aes(x = threshold, y = value, col = Metric)) + geom_line() + theme_cowplot(font_size = 11) + xlab("Threshold") + ylab("") +
  scale_color_manual(values = c("black", my_cols[1], my_cols[2]), labels = c("MSE", expression(Bias^2), "Variance"), breaks = c("MSE", expression(Bias^2), "Variance"))


# combine
p_bottom <- plot_grid(p_b, p_c, align = "h", ncol = 2, rel_widths = c(1,2), labels = c("b", "c"))
p <- plot_grid(p_a, p_bottom, ncol = 1, labels = c("a", ""))

