####################################
# 0. Load data and packages; set fps
####################################
glmeiv_offsite_dir_gasp_data <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/gasperini/data/")
gasp_offsite_dir <- .get_config_path("LOCAL_GASPERINI_2019_DATA_DIR")
gasp_gene_dir <- paste0(gasp_offsite_dir, "at-scale/processed/gene/")
gasp_gRNA_dir <- paste0(gasp_offsite_dir, "at-scale/processed/gRNA_grouped/")

library(magrittr)
library(ondisc)

# genes
gene_exp_fp <- paste0(gasp_gene_dir, "gasp_scale_gene_expressions.odm")
gene_exp_meta <- paste0(gasp_gene_dir, "gasp_scale_gene_metadata.rds")

# gRNAs
gRNA_counts_fp <- paste0(gasp_gRNA_dir, "gasp_scale_gRNA_counts_grouped.odm")
gRNA_counts_meta <- paste0(gasp_gRNA_dir, "gasp_scale_gRNA_metadata_grouped.rds")

# pairs
all_pairs <- readRDS(paste0(gasp_gRNA_dir, "pairs_grouped.rds"))

# read ODMs
gene_odm <- read_odm(gene_exp_fp, gene_exp_meta)
gRNA_odm <- read_odm(gRNA_counts_fp, gRNA_counts_meta)

################
# 1. QC on cells
################
# construct a multimodal ODM for cell QC.
multimodal_odm <- multimodal_ondisc_matrix(list(gene = gene_odm, gRNA = gRNA_odm))

# compute quantiles of mRNA and gRNA counts
count_quants <- multimodal_odm %>%
  get_cell_covariates() %>% dplyr::select(gene_n_umis, gRNA_n_umis) %>%
  dplyr::summarize(gene_quant = quantile(gene_n_umis, c(0.05, 0.95)),
                   gRNA_quant = quantile(gRNA_n_umis, c(0.05, 0.95)))

# determine which cells are OK to use:
# (i) mRNA count in [0.05, 0.95] percentile
# (ii) gRNA count in [0.05, 0.95] percentile
# (iii) p_mito < 0.08 
ok_cells <- multimodal_odm %>% get_cell_covariates() %>%
  dplyr::mutate(gene_count_ok = (gene_n_umis > count_quants$gene_quant[1] & gene_n_umis < count_quants$gene_quant[2]),
                gRNA_count_ok = (gRNA_n_umis > count_quants$gRNA_quant[1] & gRNA_n_umis < count_quants$gRNA_quant[2]),
                p_mito_ok = gene_p_mito < 0.08,
                ok = gene_count_ok & gRNA_count_ok & p_mito_ok) %>% dplyr::pull(ok)

# subset the multimodal ODM according to QC metrics
multimodal_odm_sub <- multimodal_odm[,ok_cells]

# obtain global covariate matrix
global_covariate_matrix <- multimodal_odm_sub %>%
  get_cell_covariates() %>% dplyr::mutate(lg_mRNA_lib_size = log(gene_n_umis),
                                          lg_gRNA_lib_size = log(gRNA_n_umis),
                                          batch = gene_batch,
                                          p_mito = gene_p_mito) %>%
  dplyr::select(lg_mRNA_lib_size, lg_gRNA_lib_size, batch, p_mito) %>%
  dplyr::mutate(batch = ifelse(batch == "prep_batch_1", 0L, 1L))
covariate_matrix_to_save <- global_covariate_matrix %>% dplyr::select(batch, p_mito)
m_offsets_to_save <- global_covariate_matrix %>% dplyr::pull(lg_mRNA_lib_size)
g_offsets_to_save <- global_covariate_matrix %>% dplyr::pull(lg_gRNA_lib_size)

################
# 2. QC on genes
################
# obtain the modalities
gene_odm_sub <- get_modality(multimodal_odm_sub, "gene")

# keep only genes with mean expression >= 1 and p_expressed > 0.1.
ok_genes <- gene_odm_sub %>% mutate_feature_covariates(p_expressed = n_nonzero/ncol(multimodal_odm), ok = (mean_expression >= 1 & p_expressed > 0.1)) %>% get_feature_covariates() %>% dplyr::pull(ok)

gene_odm_qc <- gene_odm_sub[ok_genes,]
gene_odm_qc_to_save <- gene_odm_qc %>% mutate_cell_covariates(n_nonzero = NULL, n_umis = NULL,
                                                              p_mito = NULL, batch = NULL)

################
# 3. QC on gRNAs
################
gRNA_odm_sub <- get_modality(multimodal_odm_sub, "gRNA")
gRNA_odm_qc_to_save <- gRNA_odm_sub %>% mutate_cell_covariates(n_nonzero = NULL, n_umis = NULL)

###########################
# 4. Subset gene-gRNA pairs
###########################
subsetted_pairs <- all_pairs %>% dplyr::filter(gene_id %in% get_feature_ids(gene_odm_qc_to_save),
                                               gRNA_id %in% get_feature_ids(gRNA_odm_qc_to_save))
set.seed(11)
sample_pairs <- subsetted_pairs %>% dplyr::filter(gene_id %in% sample(x = subsetted_pairs$gene_id, size = 4, replace = FALSE) &
                                  gRNA_id %in% sample(x = subsetted_pairs$gRNA_id, size = 4, replace = FALSE)) %>%
  dplyr::slice_sample(n = 15)
pc_pairs <- subsetted_pairs %>% dplyr::filter(site_type == "selfTSS")
sample_pairs_pc <- subsetted_pairs %>% dplyr::filter(site_type == "selfTSS") %>% dplyr::slice_sample(n = 15)

########################################################
# 5. Save ODMs, subsetted pairs, global covariate matrix
########################################################
# gene metadata
save_odm(odm = gene_odm_qc_to_save,
         paste0(glmeiv_offsite_dir_gasp_data, "gene_qc_metadata.rds"))
# gRNA metadata
save_odm(odm = gRNA_odm_qc_to_save,
         paste0(glmeiv_offsite_dir_gasp_data, "gRNA_qc_metadata.rds"))
# covariate matrix
saveRDS(object = covariate_matrix_to_save,
        paste0(glmeiv_offsite_dir_gasp_data, "covariate_matrix.rds"))
# subsetted pairs
saveRDS(object = subsetted_pairs,
        paste0(glmeiv_offsite_dir_gasp_data, "gRNA_gene_pairs.rds"))
# subsetted pc pairs
saveRDS(object = pc_pairs,
        paste0(glmeiv_offsite_dir_gasp_data, "gRNA_gene_pairs_pc.rds"))
# sample pairs
saveRDS(sample_pairs,
        paste0(glmeiv_offsite_dir_gasp_data, "gRNA_gene_pairs_sample.rds"))
# sample pairs pc
saveRDS(sample_pairs_pc,
        paste0(glmeiv_offsite_dir_gasp_data, "gRNA_gene_pairs_sample_pc.rds"))
# m offsets
saveRDS(object = global_covariate_matrix %>% dplyr::pull(lg_mRNA_lib_size),
        file = paste0(glmeiv_offsite_dir_gasp_data, "m_offsets.rds"))
# g offsets
saveRDS(object = global_covariate_matrix %>% dplyr::pull(lg_gRNA_lib_size),
        file = paste0(glmeiv_offsite_dir_gasp_data, "g_offsets.rds"))
