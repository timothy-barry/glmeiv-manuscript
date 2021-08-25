source set_xie_args.sh

#######################
# Run Nextflow pipeline
#######################
rm -f trace.txt
# nextflow run $glmeiv_nf_pipeline \
nextflow run https://github.com/timothy-barry/glmeiv-pipeline -r main \
--pairs $gRNA_gene_pairs \
--covariate_matrix $covariate_matrix \
--pair_pod_size $pair_pod_size \
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
-bg > $PWD/log -ansi-log false -with-trace
