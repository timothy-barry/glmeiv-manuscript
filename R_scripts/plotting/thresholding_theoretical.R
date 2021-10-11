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

# subplot c) variance as c -> infinity
ggplot(data = to_plot_b, mapping = aes(x = pi, y = Estimate)) + geom_line() +
  scale_x_continuous(trans = "reverse", expand = c(NA, 0), limits = c(0.5, -0.02)) + theme_cowplot(font_size = 11) +
  geom_segment(aes(x=0.0, xend=0.5, y=1, yend=1), lwd = 0.3) + scale_y_continuous(expand = c(0, NA), limits = c(0.5, 1.02)) +
  xlab(expression(pi))

# 
