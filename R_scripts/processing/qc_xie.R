library(magrittr)
library(ondisc)
glmeiv_offsite_dir_xie_data <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/xie/data/")
xie_offsite <- paste0(.get_config_path("LOCAL_XIE_2019_DATA_DIR"), "processed/")

##############
# 0. Load data
##############
gene_odm <- read_odm(odm_fp = paste0(xie_offsite, "gene/expression_matrix.odm"),
                     metadata_fp = paste0(xie_offsite, "gene/metadata.rds"))
gRNA_odm <- read_odm(odm_fp = paste0(xie_offsite, "gRNA/raw_grouped.odm"),
                     metadata_fp = paste0(xie_offsite, "gRNA/raw_grouped_metadata.rds"))

##################
# 1. QC on cells
##################
multimodal_odm <- multimodal_ondisc_matrix(list(gene = gene_odm, gRNA = gRNA_odm))

# determine which cells contain fewer than 0.08 p_mito
ok_cells <- multimodal_odm %>% get_cell_covariates() %>%
  dplyr::summarize(ok = gene_p_mito < 0.08) %>% dplyr::pull(ok)

# subset the multimodal ODM according to QC metrics
multimodal_odm_sub <- multimodal_odm[,ok_cells]

# obtain global covariate matrix
global_covariate_matrix <- multimodal_odm_sub %>%
  get_cell_covariates() %>% dplyr::mutate(lg_mRNA_lib_size = log(gene_n_umis),
                                          lg_gRNA_lib_size = log(gRNA_n_umis),
                                          batch = gene_batch,
                                          p_mito = gene_p_mito) %>%
  dplyr::select(lg_mRNA_lib_size, lg_gRNA_lib_size, batch, p_mito) %>%
  dplyr::mutate(batch = factor(batch))

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
gene_odm_qc_to_save <- gene_odm_qc %>% mutate_cell_covariates(n_nonzero = NULL, n_umis = NULL, p_mito = NULL, batch = NULL)

################
# 3. QC on gRNAs
################
gRNA_odm_sub <- get_modality(multimodal_odm_sub, "gRNA")
gRNA_odm_qc_to_save <- gRNA_odm_sub %>% mutate_cell_covariates(n_nonzero = NULL, n_umis = NULL)

#############################
# 4. Get the pairs to analyze
#############################
gRNA_gene_pairs <- readRDS(paste0(xie_offsite, "aux/pairs_grouped.rds"))
set.seed(4)
gRNA_gene_pairs_sub <- dplyr::filter(gRNA_gene_pairs, gene_id %in% get_feature_ids(gene_odm_qc_to_save))
# check for duplication of pair id
gRNA_gene_pairs_sub %>% dplyr::summarize(pair_id = paste0(gene_id, gRNA_id)) %>%
  dplyr::pull() %>% duplicated() %>% any()
# combine the cis pairs with 50,000, randomly-selected negative control pairs
gRNA_gene_sample <- rbind(gRNA_gene_pairs_sub %>% dplyr::filter(type == "cis"),
                          gRNA_gene_pairs_sub %>% dplyr::filter(type == "neg_control") %>% dplyr::slice_sample(n = 50000))

############################################################
# 5. Save the metadata RDS, offsets, global covariate matrix
############################################################
save_odm(odm = gene_odm_qc_to_save, metadata_fp = paste0(glmeiv_offsite_dir_xie_data, "gene_metadata"))
save_odm(odm = gRNA_odm_qc_to_save, metadata_fp = paste0(glmeiv_offsite_dir_xie_data, "gRNA_metadata"))
saveRDS(object = m_offsets_to_save, file = paste0(glmeiv_offsite_dir_xie_data, "m_offset.rds"))
saveRDS(object = g_offsets_to_save, file = paste0(glmeiv_offsite_dir_xie_data, "g_offset.rds"))
saveRDS(object = covariate_matrix_to_save, file = paste0(glmeiv_offsite_dir_xie_data, "covariate_matrix.rds"))
saveRDS(object = gRNA_gene_sample, file = paste0(glmeiv_offsite_dir_xie_data, "gRNA_gene_pairs.rds"))
saveRDS(object = gRNA_gene_sample %>% dplyr::slice_sample(n = 15), file = paste0(glmeiv_offsite_dir_xie_data, "gRNA_gene_pairs_trial.rds"))
