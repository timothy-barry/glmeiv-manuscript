if (!("/jet/home/timbar/R/x86_64-redhat-linux-gnu-library/4.0" %in% .libPaths())) .libPaths("/jet/home/timbar/R/x86_64-redhat-linux-gnu-library/4.0")

########################################
# 1. Load packages and command-line args
########################################
library(ondisc)
args <- commandArgs(trailingOnly = TRUE)
n_args <- length(args)
odm_fp <- args[1L]
metadata_fp <- args[2L]
covariate_matrix_fp <- args[3L]
m_offsets_fp <- args[4L]
ids <- args[seq(5L, n_args)]

####################################################
# 2. Load ODM, covariate mat, and offset; set family
####################################################
odm <- read_odm(odm_fp = odm_fp, metadata_fp = metadata_fp)
covariate_matrix <- readRDS(covariate_matrix_fp)
offset <- readRDS(m_offsets_fp)
fam_obj <- MASS::negative.binomial(NA)

####################################################
# 3. Iterate through the IDs, running precomputation
####################################################
for (id in ids) {
  counts <- as.numeric(odm[[id,]])
  precomp <- glmeiv::run_glmeiv_precomputation(y = counts,
                                               covariate_matrix = covariate_matrix,
                                               offset = offset,
                                               fam = glmeiv::augment_family_object(fam_obj))
  saveRDS(precomp$fam, paste0(id, ".rds"))
}