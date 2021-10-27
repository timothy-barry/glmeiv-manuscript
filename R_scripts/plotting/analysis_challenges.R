library(magrittr)
library(ggplot2)
library(cowplot)

my_cols <- c("firebrick3", "dodgerblue3", "orchid4", "lightskyblue3")
fig_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "glmeiv-manuscript/figures/analysis_challenges")
if (!dir.exists(fig_dir)) dir.create(fig_dir)

# a. create mixture distributions
x_grid <- seq(-1.2, 1.7, 0.01)
to_plot <- rbind(data.frame(x = x_grid,
           y = dnorm(x = x_grid, mean = 0, sd = 0.25),
           modality = "gene",
           Status = "Unperturbed"),
data.frame(x = x_grid,
           y = dnorm(x = x_grid, mean = -0.4, sd = 0.25),
           modality = "gene",
           Status = "Perturbed"),
data.frame(x = x_grid,
           y = dnorm(x = x_grid, mean = 1, sd = 0.25),
           modality = "gRNA",
           Status = "Perturbed"),
data.frame(x = x_grid,
           y = dnorm(x = x_grid, mean = 0, sd = 0.25),
           modality = "gRNA",
           Status = "Unperturbed")) %>% dplyr::mutate(modality = factor(x = modality, levels = c("gene", "gRNA")))
p1 <- ggplot(data = to_plot, mapping = aes(x = x, y = y, col = Status)) +
  facet_grid(modality~.) + ylab("Density") + xlab("Expression") + theme_bw() +
  theme(legend.position="bottom") + theme(legend.position="bottom", legend.title=element_blank()) +
  geom_vline(data = dplyr::filter(to_plot, modality == "gene"), mapping = aes(xintercept = -0.4), col = "darkred") +
  geom_vline(data = dplyr::filter(to_plot, modality == "gene"), mapping = aes(xintercept = 0.0), col = "darkred") +
  geom_hline(yintercept = 0) + scale_y_continuous(expand = c(0, 0.05)) + 
  scale_color_manual(values = c(my_cols[2], my_cols[3])) + geom_line(lwd = 1.1) +
  theme(legend.position = c(0.77, 0.85), legend.title = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line.y = element_line(colour = "black"))

# set up cowplot
p_out <- plot_grid(NULL, p1, NULL, NULL, labels = c("a", "c", "b", "d"), nrow = 2, rel_heights = c(0.5, 0.5))
f_name <- paste0(fig_dir, "/ggplot.jpg")

ggsave(filename = f_name, plot = p_out, device = "jpg", scale = 1, width = 7, height = 6, dpi = 320)
