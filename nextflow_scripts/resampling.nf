/*********************
* gene precomputations
*********************/
// Obtain the positive control pair IDs.
process obtain_pair_ids {
  time {60.s}

  output:
  stdout pair_id_ch_raw

  input:
  path pairs_fp from params.pairs

  """
  Rscript -e '
  library(magrittr);
  pairs <- readRDS("$pairs_fp") %>% dplyr::filter(site_type == "selfTSS");
  if ("${params.trial}" == "true") {
      pairs <- pairs %>% dplyr::slice(1:2);
  } else {
      set.seed(4);
      pairs <- pairs %>% dplyr::sample_n(100);
  }
  pair_names <- pairs %>% dplyr::summarize(paste0(gene_id, ":", gRNA_id)) %>% dplyr::pull();
  cat(paste0(pair_names, collapse = "\n"));'
  """
}
pair_id_ch = pair_id_ch_raw.splitText().map{it.trim()}
pair_id_ch.view()


// fit the baseline glmeiv models to the real data.
process fit_base_glmeiv_model {
  errorStrategy  { task.attempt <= 3  ? 'retry' : 'finish' }
  time { 4.m * task.attempt }

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
baseline_fit_ch.into{ baseline_fit_ch_a; baseline_fit_ch_b }

// create channel of contimination levels
my_arr = []
counter = params.seq_start
while (counter <= params.seq_end) {
  my_arr.push(counter)
  counter += params.seq_by
}
contamination_level = Channel.fromList(my_arr)
pair_contamination_product_ch = baseline_fit_ch_a.combine(contamination_level)


// fit glmeiv on resampled data for each (contamination level, pair) pair.
process fit_resampled_glmeiv {
  errorStrategy  { task.attempt <= 3  ? 'retry' : 'finish' }
  time { 4.m * params.B * task.attempt }

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


process combine_results_1 {
  time {10.m}

  input:
  file 'raw_result' from baseline_fit_ch_b.collect()

  output:
  file 'combined_result_real_data.rds' into combined_result_real_data_ch

  """
  Rscript -e 'library(magrittr)
  args <- commandArgs(trailingOnly = TRUE)
  n_args <- length(args)
  raw_result_fps <- args[seq(1L, n_args)]
  combined_result <- do.call(rbind, lapply(raw_result_fps, function(fp) readRDS(fp)[["tbl"]]))
  saveRDS(object = combined_result, file = "combined_result_real_data.rds")' raw_result*
  """
}


// combine results
process combine_results_2 {
  time {10.m}
  publishDir params.result_dir, mode: "copy"

  input:
  file 'raw_result' from raw_results_ch.collect()
  file 'combined_result_real_data' from combined_result_real_data_ch

  output:
  file "resampling_result.rds" into collected_results_ch

  """
  Rscript -e '
  library(magrittr)
  args <- commandArgs(trailingOnly = TRUE)
  n_args <- length(args)
  raw_result_fps <- args[seq(1L, n_args)]
  combined_result <- do.call(rbind, lapply(raw_result_fps, function(fp) readRDS(fp)))
  saveRDS(object = combined_result, file = "resampling_result.rds")
  ' combined_result_real_data raw_result*
  """
}
