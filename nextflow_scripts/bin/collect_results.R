if (!("/jet/home/timbar/R/x86_64-redhat-linux-gnu-library/4.0" %in% .libPaths())) .libPaths("/jet/home/timbar/R/x86_64-redhat-linux-gnu-library/4.0")

#######################
# 1. Set fps, load data
#######################
library(magrittr)
args <- commandArgs(trailingOnly = TRUE)
n_args <- length(args)
pairs_fp <- args[1L]
fp_name <- args[2L]
raw_result_fps <- args[seq(3L, n_args)]

# load data
pairs_df <- readRDS(pairs_fp)
combined_result <- do.call(rbind, lapply(raw_result_fps, function(fp) readRDS(fp)))

########################################
# 2. Append site type to combined_result
########################################
if (ncol(pairs_df) >= 3) {
  combined_result <- dplyr::mutate(combined_result, pair_id = factor(paste0(gene_id, ":", gRNA_id)))
  pairs_df <- dplyr::mutate(pairs_df, pair_id = factor(paste0(gene_id, ":", gRNA_id))) %>%
    dplyr::mutate(gene_id = NULL, gRNA_id = NULL)
  out <- dplyr::left_join(x = combined_result, y = pairs_df, by = "pair_id")
} else {
  out <- combined_result
}

saveRDS(out, fp_name)
