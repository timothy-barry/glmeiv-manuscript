# set figure 2 dir
fig_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "glmeiv-manuscript/figures/")
fig2_dir <- paste0(fig_dir, "fig2/")
if (!dir.exists(fig2_dir)) dir.create(fig2_dir, recursive = TRUE)
load_all("~/research_code/glmeiv/")

# load plotting objects
plotting_objects <- paste0(.get_config_path("LOCAL_CODE_DIR"), "glmeiv-manuscript/figures/plotting_objects.R")
source(plotting_objects)

# a. gRNA count histogram
set.seed(2)
gRNA_ids_gasp <- sample(x = get_feature_ids(gRNA_odm_gasp), size = 2, replace = FALSE)
gRNA_ids_xie <- sample(x = get_feature_ids(gRNA_odm_xie), size = 2, replace = FALSE)

to_plot <- rbind(data.frame(count = as.numeric(gRNA_odm_gasp[[gRNA_ids_gasp[1],]]), gRNA_id = gRNA_ids_gasp[1], dataset = "Gasp"),
                 data.frame(count = as.numeric(gRNA_odm_gasp[[gRNA_ids_gasp[2],]]), gRNA_id = gRNA_ids_gasp[2], dataset = "Gasp"),
                 data.frame(count = as.numeric( gRNA_odm_xie[[gRNA_ids_xie[1],]]), gRNA_id = gRNA_ids_xie[1], dataset = "Xie"),
                 data.frame(count = as.numeric( gRNA_odm_xie[[gRNA_ids_xie[2],]]), gRNA_id = gRNA_ids_xie[2], dataset = "Xie")) %>% 
  dplyr::mutate(lab = factor(x = gRNA_id,
                             levels = c(gRNA_ids_gasp, gRNA_ids_xie),
                             labels = c("Gasperini gRNA 1", "Gasperini gRNA 2", "Xie gRNA 1", "Xie gRNA 2"))) %>% dplyr::filter(count >= 1)

# compute the original thresholds
thresh_df <- to_plot %>% dplyr::group_by(lab) %>% dplyr::group_modify(.f = function(tbl, key) {
  dataset <- tbl$dataset[1]
  if (dataset == "Gasp") {
    thresh <- 5
  } else {
    curr_count <- tbl$count
    curr_count <- curr_count[curr_count >= 1]
    thresh <- mean(curr_count)
  }
  return(tibble::tibble(thresh = thresh))
})

my_cols <- c("firebrick3", "dodgerblue3")

p1 <- ggplot(data = to_plot %>% dplyr::filter(count >= 1, count <= 40), mapping = aes(x = count)) +
  facet_wrap(lab~., scales = "free_y") + geom_histogram(binwidth = 1, col = "black", fill = "gray") +
  scale_y_continuous(trans='log10', expand = c(0, NA)) + xlab("gRNA count") + ylab("") +
  geom_vline(data = dplyr::filter(to_plot, lab == "Gasperini gRNA 1"),
             mapping = aes(xintercept = dplyr::filter(thresh_df, lab == "Gasperini gRNA 1") %>% dplyr::pull(thresh)), 
             col = my_cols[1], lwd = 1) +
  geom_vline(data = dplyr::filter(to_plot, lab == "Gasperini gRNA 2"),
             mapping = aes(xintercept = dplyr::filter(thresh_df, lab == "Gasperini gRNA 2") %>% dplyr::pull(thresh)), 
             col = my_cols[1], lwd = 1) +
  geom_vline(data = dplyr::filter(to_plot, lab == "Xie gRNA 1"),
             mapping = aes(xintercept = dplyr::filter(thresh_df, lab == "Xie gRNA 1") %>% dplyr::pull(thresh) %>% ceiling()), 
             col = my_cols[2], lwd = 1) +
  geom_vline(data = dplyr::filter(to_plot, lab == "Xie gRNA 2"),
             mapping = aes(xintercept = dplyr::filter(thresh_df, lab == "Xie gRNA 2") %>% dplyr::pull(thresh) %>% ceiling()), 
             col = my_cols[2], lwd = 1) + theme_bw()
ggsave(filename = paste0(fig2_dir, "histograms.pdf"), plot = p1, device = "pdf", scale = 1, width = 6.5, height = 3.75)

# b. theoretical bias of thresholding estimator
xgrid <- seq(0.2, 6, by = 0.1)
y1 <- get_tresholding_estimator_bias(m_perturbation = 1,
                                             g_perturbation = xgrid,
                                             pi = 0.2,
                                             return_bias = FALSE)
y2 <- get_tresholding_estimator_bias(m_perturbation = 1,
                                             g_perturbation = xgrid,
                                             pi = 0.4,
                                             return_bias = FALSE)
to_plot <- rbind(data.frame(y = y1, x = xgrid, pi = 0.2),
                 data.frame(y = y2, x = xgrid, pi = 0.4)) %>% dplyr::mutate(pi = factor(pi))
p2 <- ggplot(data = to_plot, mapping = aes(x = x, y = y, col = pi)) +
  geom_hline(yintercept = 1) + scale_x_continuous(expand = c(0.01, 0.01)) +
  geom_line(lwd = 0.8) + xlab(expression(beta[g])) + ylab("Estimate") +
  theme_bw() + scale_color_manual(values = my_cols) +
  theme(legend.position = c(0.85, 0.25),
        legend.background = element_rect(fill = "white", color = "black")) +
  guides(color = guide_legend(title = bquote(pi)))
ggsave(filename = paste0(fig2_dir, "att_bias.pdf"), plot = p2, device = "pdf", scale = 1, width = 3.5, height = 2.75)
