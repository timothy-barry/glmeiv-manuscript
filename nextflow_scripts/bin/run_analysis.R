##############################
# 0. Set a few hyperparameters
##############################
threshold <- "adaptive"

########################################
# 1. Load packages and command-line args
########################################
library(magrittr)
library(ondisc)
args <- commandArgs(trailingOnly = TRUE)
n_args <- length(args)
covariate_matrix_fp <- args[1L]
gene_odm_fp <- args[2L]
gene_metadata_fp <- args[3L]
m_offsets_fp <- args[4L]
gRNA_odm_fp <- args[5L]
gRNA_metadata_fp <- args[6L]
other_args <- args[seq(7L, n_args)]

#########################################
# 2. Load ODMs, covariate matrix, offsets
#########################################
gene_odm <- read_odm(gene_odm_fp, gene_metadata_fp)
gRNA_odm <- read_odm(gRNA_odm_fp, gRNA_metadata_fp)
m_offset <- readRDS(m_offsets_fp)
covariate_matrix <- readRDS(covariate_matrix_fp)

#############################################################
# 3. Set vectors of gene IDs, gRNA IDs, and precomp locations
#############################################################
gene_ids <- other_args[seq(from = 1L, by = 3, to = length(other_args))]
gRNA_ids <- other_args[seq(from = 2L, by = 3, to = length(other_args))]
gene_precomp_fps <- other_args[seq(from = 3L, by = 3, to = length(other_args))]
n_pairs <- length(gene_ids)

############################################################
# 4. Loop through pairs, running method given precomputation
############################################################
out_l <- vector(mode = "list", length = n_pairs)
for (i in seq(1L, n_pairs)) {
  gene <- gene_ids[i]
  if (i == 1 || gene_ids[i] != gene_ids[i - 1]) { # only load gene data if necessary
    m <- as.numeric(gene_odm[[gene,]]); m_fam <- readRDS(gene_precomp_fps[i])
  }
  gRNA <- gRNA_ids[i]
  if (i == 1 || gRNA_ids[i] != gRNA_ids[i - 1]) { # likewise for gRNAs
    g <- as.numeric(gRNA_odm[[gRNA,]])
  }
  # loop through different thresholds, fitting model to each
  # perform threshold
  if (is.numeric(threshold)) {
    curr_threshold <- threshold
  } else {
    g_prime <- g[g >= 1]
    curr_threshold <- mean(g_prime)
  }
  phat <- as.integer(g >= curr_threshold)
  # run thresholding method
  fit <- glmeiv::run_thresholding_method(phat = phat, m = m, m_fam = m_fam,
                                         m_offset = m_offset, covariate_matrix = covariate_matrix,
                                         n_examples_per_param = 5) %>% dplyr::mutate(gene_id = gene, gRNA_id = gRNA)
  out_l[[i]] <- fit
}

out <- do.call(rbind, out_l) %>% dplyr::mutate_at(c("parameter", "target", "gene_id", "gRNA_id"), factor)
saveRDS(object = out, file = "raw_result.rds")
