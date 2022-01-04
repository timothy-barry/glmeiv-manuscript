source set_gasp_args.sh

#######################
# Run Nextflow pipeline
#######################
if [ $trial = true ]
then
  gRNA_gene_pairs_pc=$processed_data_dir"gRNA_gene_pairs_sample_pc.rds"
else
  gRNA_gene_pairs_pc=$processed_data_dir"gRNA_gene_pairs_pc.rds"
fi

rm -f trace.txt
nextflow run $thresholding_nf_pipeline --pairs $gRNA_gene_pairs_pc \
--covariate_matrix $covariate_matrix \
--pair_pod_size 50 \
--result_dir $result_dir \
--out_file_name "thresholding_result_pc.rds" \
--gene_pod_size 20 \
--gene_odm $gene_odm \
--gene_metadata $gene_metadata \
--m_offsets $m_offsets \
--gRNA_odm $gRNA_odm \
--gRNA_metadata $gRNA_metadata \
--threshold "1,5,20" \
-bg > $PWD/log -ansi-log false -with-trace -w $work_dir
