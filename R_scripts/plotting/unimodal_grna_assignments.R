library(ggplot2)
library(cowplot)
# create plots related to unimodal gRNA assignments

to_read_dir <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/aux/")
result <- readRDS(file = paste0(to_read_dir, "papalexi_unimodal_grna_assignments.rds"))

my_cols <- c("darkorchid1", "firebrick2")
# p_a: a plot of gRNA count vs. assignment (for both mixture model and ground truth)
grna_idx <- 2
grna_id <- result[[grna_idx]]$curr_grna_id
g <- result[[grna_idx]]$g
mixture_assignments <- result[[grna_idx]]$mixture_assignments
ground_truth_assignments <- result[[grna_idx]]$ground_truth_assignments
df <- data.frame(g = g,
                 mixture_assignment = mixture_assignments,
                 ground_truth_assignment = ground_truth_assignments)
p_a <- ggplot(data = df, mapping = aes(x = mixture_assignment, y = g, col = ground_truth_assignment)) +
  geom_jitter(alpha = 0.8) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 5, 20, 50, 200, 1000)) +
  labs(color = "Ground truth\n assignment") +
  xlab("Unimodal GLM-EIV assignment") +
  ylab("gRNA count") +
  ggtitle(paste0("Perturbation assignments\nfor gRNA ", grna_id)) +
  sceptre:::get_my_theme() +
  scale_color_manual(values = my_cols)

# p_b: get the confusion matrix and confusion matrix metrics for this gRNA
confusion_mat <- result[[grna_idx]]$confusion_mat
confusion_mat$byClass[["Sensitivity"]]
confusion_mat$byClass[["Specificity"]]
confusion_mat$byClass[["Balanced Accuracy"]]

# p_c: plot the sensitivity, specificity, and accuracy across gRNAs
sensitivities <- sapply(X = result, FUN = function(l) l$confusion_mat$byClass[["Sensitivity"]]) |> unlist()
specificities <- sapply(X = result, FUN = function(l) l$confusion_mat$byClass[["Specificity"]]) |> unlist()
accuracies <- sapply(X = result, FUN = function(l) l$confusion_mat$byClass[["Balanced Accuracy"]]) |> unlist()
df_2 <- data.frame(Sensitivity = sensitivities,
                 Specificity = specificities,
                 "Balanced_accuracy" = accuracies) |>
  dplyr::rename("Balanced accuracy" = "Balanced_accuracy") |>
  tidyr::pivot_longer(cols = c("Sensitivity", "Specificity", "Balanced accuracy"),
                      names_to = "metric")
df_2 |> dplyr::group_by(metric) |> dplyr::summarize(m = mean(value))

p_c <- ggplot(data = df_2, mapping = aes(y = value)) +
  geom_boxplot(fill = "dodgerblue", outlier.size = 0.8) + facet_wrap(metric ~ ., scales = "free_y") +
  sceptre:::get_my_theme() +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle("Classification metrics across all gRNAs")


# create cowplot
p <- plot_grid(p_a, NULL, p_c, NULL, labels = c("a", "b", "c", "d"), nrow = 2, rel_widths = c(0.6, 0.3), rel_heights = c(0.6, 0.4))

# also, compute these metrics for the replogle method
sensitivities <- sapply(X = result, FUN = function(l) l$confusion_mat_simple$byClass[["Sensitivity"]]) |> unlist()
specificities <- sapply(X = result, FUN = function(l) l$confusion_mat_simple$byClass[["Specificity"]]) |> unlist()
accuracies <- sapply(X = result, FUN = function(l) l$confusion_mat_simple$byClass[["Balanced Accuracy"]]) |> unlist()

# save
fig_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "glmeiv-manuscript/figures/unimodal_glmeiv_assignments")
if (!dir.exists(fig_dir)) dir.create(fig_dir)
ggsave(filename = paste0(fig_dir, "/r_plot.png"), plot = p, device = "png", dpi = 330, width = 6, height = 5, scale = 1.1)
