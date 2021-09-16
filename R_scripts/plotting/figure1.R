# set figure 1 dir
fig_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "glmeiv-manuscript/figures/")
fig1_dir <- paste0(fig_dir, "fig1/")
if (!dir.exists(fig1_dir)) dir.create(fig1_dir, recursive = TRUE)

# load plotting objects
plotting_objects <- paste0(.get_config_path("LOCAL_CODE_DIR"), "glmeiv-manuscript/figures/plotting_objects.R")
source(plotting_objects)

# a. create mixture distributions
x_grid <- seq(-1.2, 1.7, 0.01)
to_plot <- rbind(data.frame(x = x_grid,
           y = dnorm(x = x_grid, mean = 0, sd = 0.25),
           modality = "mRNA",
           Status = "Unperturbed"),
data.frame(x = x_grid,
           y = dnorm(x = x_grid, mean = -0.4, sd = 0.25),
           modality = "mRNA",
           Status = "Perturbed"),
data.frame(x = x_grid,
           y = dnorm(x = x_grid, mean = 1, sd = 0.25),
           modality = "gRNA",
           Status = "Perturbed"),
data.frame(x = x_grid,
           y = dnorm(x = x_grid, mean = 0, sd = 0.25),
           modality = "gRNA",
           Status = "Unperturbed")) %>% dplyr::mutate(modality = factor(x = modality, levels = c("mRNA", "gRNA")))
my_cols <- c("royalblue4", "seagreen4")
p1 <- ggplot(data = to_plot, mapping = aes(x = x, y = y, col = Status)) +
  facet_grid(modality~.) + ylab("Density") + xlab("Expression") + theme_bw() +
  theme(legend.position="bottom") + theme(legend.position="bottom", legend.title=element_blank()) +
  geom_vline(data = dplyr::filter(to_plot, modality == "mRNA"), mapping = aes(xintercept = -0.4), col = "darkred") +
  geom_vline(data = dplyr::filter(to_plot, modality == "mRNA"), mapping = aes(xintercept = 0.0), col = "darkred") +
  scale_color_manual(values = my_cols) + geom_line(lwd = 1.1)
ggsave(filename = paste0(fig1_dir, "density_plot.pdf"), plot = p1, device = "pdf", scale = 1, width = 4, height = 4)


# b. plot histogram of gene counts
to_plot <- rbind(data.frame(dataset = "Gasperini", counts = as.numeric(gene_odm_gasp[["ENSG00000188290",]])),
                 data.frame(dataset = "Xie", counts = as.numeric(gene_odm_xie[["ENSG00000240409.1",]]))) # %>% dplyr::filter(counts >= 1)
p2 <- ggplot(data = to_plot, mapping = aes(x = counts)) + geom_histogram(binwidth = 0.5, col = my_cols[1], fill = "black", lwd = 0.7) + facet_grid(dataset~.) + ylab("") + xlab("N transcripts") + ylab("N cells") + xlim(0, 15) + scale_y_continuous(trans='log10') + theme_bw()
ggsave(filename = paste0(fig1_dir, "mRNA_count_hist.pdf"), plot = p2, device = "pdf", scale = 1, width = 4.5, height = 3.5)
