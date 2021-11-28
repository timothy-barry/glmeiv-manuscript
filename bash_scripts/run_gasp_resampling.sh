source set_gasp_args.sh

#######################
# Run Nextflow pipeline
#######################
rm -f trace.txt

nextflow run $resampling_nf_pipeline \
--pairs $gRNA_gene_pairs \
--covariate_matrix $covariate_matrix \
--result_dir $result_dir \
--gene_pod_size $gene_pod_size \
--gene_odm $gene_odm \
--gene_metadata $gene_metadata \
--m_offsets $m_offsets \
--m_fam_str $m_fam \
--gRNA_pod_size $gRNA_pod_size \
--gRNA_odm $gRNA_odm \
--gRNA_metadata $gRNA_metadata \
--g_offsets $g_offsets \
--g_fam_str $g_fam \
--B 5 \
--seq_start 0 \
--seq_end 0.3 \
--seq_by 0.05 \
-ansi-log false -bg > $PWD/log -with-trace -w $work_dir -resume
