load_all("~/research_code/glmeiv")
library(caret)

# prepare the data
papalexi_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/papalexi/eccite_screen/")
gene_odm_fp <- paste0(papalexi_dir, "gene/matrix.odm")
grna_odm_fp <- paste0(papalexi_dir, "grna_expression/matrix.odm")
mm_metadata_fp <- paste0(papalexi_dir, "multimodal_metadata.rds")
mm_odm <- ondisc::read_multimodal_odm(odm_fps = c(gene_odm_fp, grna_odm_fp),
                                      multimodal_metadata_fp = mm_metadata_fp)
covariate_mat <- mm_odm |> ondisc::get_cell_covariates()
covariate_mat_trans <- covariate_mat |>
  dplyr::select(gene_n_nonzero, gene_n_umis, grna_expression_n_nonzero, grna_expression_n_umis, bio_rep) |>
  dplyr::mutate(gene_n_nonzero = log(gene_n_nonzero),
                gene_n_umis = log(gene_n_umis),
                grna_expression_n_nonzero = log(grna_expression_n_nonzero),
                grna_expression_n_umis = log(grna_expression_n_umis))

odm <- mm_odm |> ondisc::get_modality("grna_expression")
grna_mat <- odm[[1:nrow(odm),]] |> as.matrix()
grna_ids_str <- ondisc::get_feature_ids(odm)

# 0. identify the top 95% most highly expressed gRNAs -- analyze these
n_nonzero_per_grna <- rowSums(grna_mat >= 1)
cutoff <- quantile(n_nonzero_per_grna, probs = 0.05)
grnas_to_keep <- which(n_nonzero_per_grna > cutoff)
cellwise_grna_counts <- colSums(grna_mat)

# 1. loop over the grnas
output <- sapply(grnas_to_keep, function(grna_id) {
  print(grna_id)
  g <- grna_mat[grna_id,]
  curr_grna_id <- grna_ids_str[grna_id]
  
  # first, assign based on UMI counts; if g makes up more than 25% of the reads, assign 
  ratios <- g/cellwise_grna_counts
  ground_truth_assignments <- ratios >= 0.25
  
  if (all(!ground_truth_assignments)) {
    return(NULL)
  }
  
  # second, assign based on the mixture model
  res <- assign_grnas_unimodal_mixture(g = g, g_fam = poisson() |> augment_family_object(),
                                       covariate_matrix = covariate_mat_trans, g_offset = NULL)
  mixture_assignments <- res$phat
  # third assign according to the simple two-component gaussian
  g_trans <- log2(g + 1)
  simple_assignments <- tryCatch({
    fit <- flexmix::flexmix(g_trans ~ 1, k = 2)
    ifelse(fit@cluster == 2, "Perturbed", "Unperturbed") |> factor(levels = c("Perturbed", "Unperturbed"))
  }, error = function(e) {
    rep("Unperturbed", length(g_trans)) |> factor(levels = c("Perturbed", "Unperturbed"))
  })
  
  # compute the confusion matrix
  ground_truth_assignments <- ifelse(ground_truth_assignments, "Perturbed", "Unperturbed") |> factor(levels = c("Perturbed", "Unperturbed"))
  mixture_assignments <- ifelse(mixture_assignments == 0, "Unperturbed", "Perturbed") |> factor(levels = c("Perturbed", "Unperturbed"))

  confusion_mat <- confusionMatrix(data = mixture_assignments, reference = ground_truth_assignments) 
  confusion_mat_simple <- confusionMatrix(data = simple_assignments, reference = ground_truth_assignments)
  
  return(list(confusion_mat = confusion_mat, confusion_mat_simple = confusion_mat_simple, g = g,
              mixture_assignments = mixture_assignments,
              ground_truth_assignments = ground_truth_assignments,
              curr_grna_id = curr_grna_id))
}, simplify = FALSE)

is_null <- sapply(output, is.null)
output <- output[!is_null]

# save the output
to_save_dir <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/aux/")
saveRDS(object = output, file = paste0(to_save_dir, "papalexi_unimodal_grna_assignments.rds"))
