library(dplyr)
gasp_result_dir <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/gasperini/results/")
glmeiv_res <- readRDS(paste0(gasp_result_dir, "glmeiv_result.rds"))

# n pairs of each site type
site_type_table <- glmeiv_res %>% group_by(pair_id) %>% summarize(site_type = site_type[1]) %>% pull() %>% table()
n_pairs <- sum(site_type_table)

# write summarize function
summarize_res <- function(df) {
  summarized_df <- df %>% group_by(pair_id) %>%
    summarize(converged = value[target == "converged"],
              spread = value[target == "membership_probability_spread"],
              n_approx_1 = value[target == "n_approx_1"],
              estimate = value[parameter == "m_perturbation" & target == "estimate"],
              confint_lower = value[parameter == "m_perturbation" & target == "confint_lower"],
              confint_upper = value[parameter == "m_perturbation" & target == "confint_upper"],
              gene_id = gene_id[1], gRNA_id = gRNA_id[1])%>%
    mutate(covered = confint_lower < 1 & confint_upper > 1)
  summary_stats <- summarized_df %>% summarize(m_est = mean(estimate),
                             n_ci = sum(!is.na(covered)),
                             ci_coverage = mean(covered, na.rm = TRUE),
                             m_lower_ci = mean(confint_lower, na.rm = TRUE),
                             m_upper_ci = mean(confint_upper, na.rm = TRUE))
  return(list(summarized_df = summarized_df, summary_stats = summary_stats))
}

# examine NTCs (ground truth)
ntc_res <- summarize_res(filter(glmeiv_res, site_type == "NTC"))
ntc_pairs_filtered <- filter(ntc_res$summarized_df, estimate <= 1.2, estimate >= 0.8)
m_hat <- mean(ntc_pairs_filtered$estimate)
sd_hat <- sd(ntc_pairs_filtered$estimate)
ntc_pairs_filtered %>% pull(estimate) %>% hist(breaks = 28, probability = TRUE)
abline(v = m_hat, col = "darkblue", lwd = 1.5)
xs <- seq(0.8, 1.2, length.out = 100)
ys <- dnorm(x = xs, mean = m_hat, sd = sd_hat)
lines(x = xs, y = ys, col = "darkred", lwd = 1.5)
ntc_res$summary_stats
100 * nrow(filter(ntc_res$summarized_df, estimate <= 0.80))/nrow(ntc_res$summarized_df)

# examine selfTSS (ground truth)
self_tss_res <- summarize_res(filter(glmeiv_res, site_type == "selfTSS"))
self_tss_res$summary_stats
self_tss_res$summarized_df

# examine DHS (hardest, as there is no ground truth)
dhs_res <- summarize_res(filter(glmeiv_res, site_type == "DHS"))
dhs_res_filtered <- filter(dhs_res$summarized_df, estimate <= 0.80)
# pairs with a substantial shift in expression
dhs_discoveries <- filter(dhs_res$summarized_df, estimate <= 0.80)