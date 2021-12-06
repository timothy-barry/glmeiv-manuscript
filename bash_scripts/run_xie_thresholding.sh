source set_xie_args.sh

#######################
# Run Nextflow pipeline
#######################
rm -f trace.txt
nextflow run $thresholding_nf_pipeline --pairs $gRNA_gene_pairs \
--covariate_matrix $covariate_matrix \
--pair_pod_size $pair_pod_size \
--result_dir $result_dir \
--gene_pod_size $gene_pod_size \
--gene_odm $gene_odm \
--gene_metadata $gene_metadata \
--m_offsets $m_offsets \
--gRNA_odm $gRNA_odm \
--gRNA_metadata $gRNA_metadata \
--threshold $thresh \
-bg > $PWD/log -ansi-log false -with-trace -w $work_dir
