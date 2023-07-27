# plotting objects

# 1. packages
library(ondisc)
library(ggplot2)
library(magrittr)

# 2. Load the ondisc matrices of gasperini and xie gene/gRNA data
xie_2019_offsite <- paste0(.get_config_path("LOCAL_XIE_2019_DATA_DIR"), "processed/")
xie_glmeiv_offsite <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/xie/data/")
gasp_2019_offsite <- paste0(.get_config_path("LOCAL_GASPERINI_2019_DATA_DIR"), "at-scale/processed/")
gasp_glmeiv_offsite <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/gasperini/data/")

gene_odm_xie <- read_odm(odm_fp = paste0(xie_2019_offsite, "gene/expression_matrix.odm"),
                         metadata_fp = paste0(xie_glmeiv_offsite, "gene_metadata.rds"))
gRNA_odm_xie <- read_odm(odm_fp = paste0(xie_2019_offsite, "gRNA/raw_grouped.odm"),
                         metadata_fp = paste0(xie_2019_offsite, "gRNA/raw_grouped_metadata.rds"))
gene_odm_gasp <- gene_odm <- read_odm(odm_fp = paste0(gasp_2019_offsite, "gene/gasp_scale_gene_expressions.odm"),
                                      metadata_fp = paste0(gasp_glmeiv_offsite, "gene_qc_metadata.rds"))
gRNA_odm_gasp <- read_odm(odm_fp = paste0(gasp_2019_offsite, "gRNA_grouped/gasp_scale_gRNA_counts_grouped.odm"),
                          metadata_fp = paste0(gasp_glmeiv_offsite, "gRNA_qc_metadata.rds"))
