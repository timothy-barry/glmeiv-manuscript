/*********************
* gene precomputations
*********************/
// Obtain the positive control pair IDs.
process obtain_pair_ids {
  time "60s"

  output:
  stdout pair_id_ch_raw

  input:
  path pairs_fp from params.pairs

  """
  Rscript -e '
  library(magrittr);
  pairs <- readRDS("$pairs_fp");
  pair_names <- pairs %>% dplyr::filter(site_type == "selfTSS") %>%
  dplyr::summarize(paste0(gene_id, ":", gRNA_id)) %>% dplyr::slice(seq(1, 2)) %>% dplyr::pull()
  cat(paste0(pair_names, collapse = "\n"))'
  """
}
pair_id_ch = pair_id_ch_raw.splitText().map{it.trim()}


// fit the baseline glmeiv models to the real data.
process fit_base_glmeiv_model {
  time 2.m
  echo true

  output:
  file 'baseline_fit.rds' into baseline_fit_ch

  input:
  path covariate_matrix_fp from params.covariate_matrix
  path gene_odm_fp from params.gene_odm
  path gene_metadata_fp from params.gene_metadata
  path m_offsets_fp from params.m_offsets
  path gRNA_odm_fp from params.gRNA_odm
  path gRNA_metadata_fp from params.gRNA_metadata
  path g_offsets_fp from params.g_offsets
  val pair_id from pair_id_ch

  """
  run_baseline_fit.R $covariate_matrix_fp $gene_odm_fp $gene_metadata_fp $m_offsets_fp $gRNA_odm_fp $gRNA_metadata_fp $g_offsets_fp $pair_id
  """
}

// create channel of contimination levels
my_arr = []
counter = params.seq_start
while (counter <= params.seq_end) {
  my_arr.push(counter)
  counter += params.seq_by
}
contamination_level = Channel.fromList(my_arr)
pair_contamination_product_ch = baseline_fit_ch.combine(contamination_level)

// fit glmeiv on resampled data for each (contamination level, pair) pair.
process fit_resampled_glmeiv {
  time {2.m * params.B}

  output:
  file 'raw_result.rds' into raw_results_ch

  input:
  path covariate_matrix_fp from params.covariate_matrix
  path gene_odm_fp from params.gene_odm
  path gene_metadata_fp from params.gene_metadata
  path m_offsets_fp from params.m_offsets
  path g_offsets_fp from params.g_offsets
  tuple file('base_fit.rds'), val(contamination_level) from pair_contamination_product_ch

  """
  fit_resampled_glmeiv.R $covariate_matrix_fp $gene_odm_fp $gene_metadata_fp $m_offsets_fp $g_offsets_fp base_fit.rds $contamination_level ${params.B}
  """
}

// combine results
