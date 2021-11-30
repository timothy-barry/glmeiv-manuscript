library(magrittr)
library(ondisc)
glmeiv_offsite_dir_xie_data <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/xie/data/")
xie_offsite <- paste0(.get_config_path("LOCAL_XIE_2019_DATA_DIR"), "processed/")

##############
# 1. Load data
##############
gene_odm <- read_odm(odm_fp = paste0(xie_offsite, "gene/expression_matrix.odm"),
                     metadata_fp = paste0(xie_offsite, "gene/metadata.rds"))
gRNA_odm <- read_odm(odm_fp = paste0(xie_offsite, "gRNA/raw_grouped.odm"),
                     metadata_fp = paste0(xie_offsite, "gRNA/raw_grouped_metadata.rds"))

#############################
# 2. QC on genes (no gRNA qc)
#############################
highly_expressed_genes <- gene_odm %>% get_feature_covariates() %>%
  dplyr::mutate(p_exp = n_nonzero / ncol(gene_odm)) %>%
  dplyr::filter(p_exp >= 0.08) %>% row.names()
gene_odm <- gene_odm[highly_expressed_genes,]

##########################################
# 3. Global cell covariate matrix, offsets
##########################################
multimodal_odm <- multimodal_ondisc_matrix(list(gene = gene_odm, gRNA = gRNA_odm))
global_covariate_matrix <- multimodal_odm %>% get_cell_covariates() %>%
  dplyr::select(batch = gene_batch)
m_offset <- multimodal_odm %>% get_cell_covariates() %>% dplyr::pull(gene_n_umis) %>% log()
g_offset <- multimodal_odm %>% get_cell_covariates() %>% dplyr::pull(gRNA_n_umis) %>% log()

##################
# 4. Get the pairs
##################
gRNA_gene_pairs <- readRDS(paste0(xie_offsite, "aux/pairs_grouped.rds"))
set.seed(4)
gRNA_gene_pairs_sub <- dplyr::filter(gRNA_gene_pairs, gene_id %in% get_feature_ids(gene_odm))
# check for duplication of pair id
gRNA_gene_pairs_sub %>% dplyr::summarize(pair_id = paste0(gene_id, gRNA_id)) %>%
  dplyr::pull() %>% duplicated() %>% any()
# combine the cis pairs with 50,000, randomly-selected negative control pairs
gRNA_gene_sample <- rbind(gRNA_gene_pairs_sub %>% dplyr::filter(type == "cis"),
                          gRNA_gene_pairs_sub %>% dplyr::filter(type == "neg_control") %>% dplyr::slice_sample(n = 50000))

############################################################
# 5. Save the metadata RDS, offsets, global covariate matrix
############################################################
save_odm(odm = gene_odm,
         metadata_fp = paste0(glmeiv_offsite_dir_xie_data, "gene_metadata"))
save_odm(odm = gRNA_odm,
         metadata_fp = paste0(glmeiv_offsite_dir_xie_data, "gRNA_metadata"))
saveRDS(object = m_offset, file = paste0(glmeiv_offsite_dir_xie_data, "m_offset.rds"))
saveRDS(object = g_offset, file = paste0(glmeiv_offsite_dir_xie_data, "g_offset.rds"))
saveRDS(object = global_covariate_matrix, file = paste0(glmeiv_offsite_dir_xie_data, "covariate_matrix.rds"))
saveRDS(object = gRNA_gene_sample, file = paste0(glmeiv_offsite_dir_xie_data, "gRNA_gene_pairs.rds"))
saveRDS(object = gRNA_gene_sample %>% dplyr::slice_sample(n = 15), file = paste0(glmeiv_offsite_dir_xie_data, "gRNA_gene_pairs_trial.rds"))
