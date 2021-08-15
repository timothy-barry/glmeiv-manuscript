source ~/.research_config

# check for trial flag -t on command line
trial=false
while getopts t OPT
do
    case "$OPT" in
        t) trial=true ;;
    esac
done


#############################
# 1. Locate nextflow pipeline
#############################
# location of Nextflow pipeline
nf_pipeline=$LOCAL_CODE_DIR"glmeiv-pipeline/main.nf"

######################
# 2. Set all arguments
######################
# i. Location of backing .odm files
backing_files_dir=$LOCAL_GASPERINI_2019_DATA_DIR"at-scale/processed/"
gene_odm=$backing_files_dir"gene/gasp_scale_gene_expressions.odm"
gRNA_odm=$backing_files_dir"gRNA_grouped/gasp_scale_gRNA_counts_grouped.odm"
# ii. Location of processed data
processed_data_dir=$LOCAL_GLMEIV_DATA_DIR"public/gasperini/data/"
# iii. metadata files
gene_metadata=$processed_data_dir"gene_qc_metadata.rds"
gRNA_metadata=$processed_data_dir"gRNA_qc_metadata.rds"
# iv. covariate matrix and offsets
covariate_matrix=$processed_data_dir"covariate_matrix.rds"
m_offsets=$processed_data_dir"m_offsets.rds"
g_offsets=$processed_data_dir"g_offsets.rds"
# v. family strings
m_fam="nb"
g_fam="poisson"
# vi. pairs to analyze and pod sizes
if [ $trial = true ]
then
  gRNA_gene_pairs=$processed_data_dir"gRNA_gene_pairs_sample.rds"
  gene_pod_size=3 #500
  gRNA_pod_size=3 #500
  pair_pod_size=3 #500
else
  gRNA_gene_pairs=$processed_data_dir"gRNA_gene_pairs.rds"
  gene_pod_size=500
  gRNA_pod_size=200
  pair_pod_size=500
fi

# viii. results directory
result_dir=$LOCAL_GLMEIV_DATA_DIR"public/gasperini/results"

##########################
# 3. Run Nextflow pipeline
##########################
rm -f trace.txt
nextflow run $nf_pipeline --pairs $gRNA_gene_pairs \
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
