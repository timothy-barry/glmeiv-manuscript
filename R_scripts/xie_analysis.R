library(magrittr)
library(dplyr)
library(tidyr)
glmeiv_offsite_dir <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/")

res <- readRDS(paste0(glmeiv_offsite_dir, "xie/results/glmeiv_result_nb10.rds"))
res_sum <- res %>% filter(type == "neg_control", parameter == "m_perturbation") %>%
  select(target, value, pair_id) %>% pivot_wider(id_cols = "pair_id", names_from = "target", values_from = "value") %>%
  mutate(cover = (confint_lower <= 1 & confint_upper >= 1))

res_sum_cis <- res %>% filter(type == "cis", parameter == "m_perturbation") %>% 
  select(target, value, pair_id) %>% pivot_wider(id_cols = "pair_id", names_from = "target", values_from = "value") %>%
  mutate(cover = (confint_lower <= 1 & confint_upper >= 1)) %>% filter(estimate < 1.5, !is.na(confint_lower))

# res sum
res_sum <- filter(res_sum, estimate < 1.5, !is.na(confint_lower))
hist(res_sum$estimate, breaks = 30, freq = FALSE, xlim = c(0.75, 1.25))
abline(v = 1, col = "red")
ests <- res_sum$estimate
mu <- mean(ests); sigma <- sd(ests)
xs <- seq(0.6, 1.4, length.out = 1000)
ys <- dnorm(x = xs, mean = mu, sd = sigma)
lines(x = xs, y = ys, col = "blue")
res_sum %>% summarize(coverage = mean(cover, na.rm = TRUE),
                      p_below_0.8 = 100 * mean(estimate < 0.8))

# res sum cis
res_sum_cis %>% summarize(coverage = mean(cover, na.rm = TRUE),
                          p_below_0.8 = 100 * mean(estimate < 0.8))
