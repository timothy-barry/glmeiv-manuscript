source ~/.research_config

# use -t flag to set "trial" parameters; otherwise, sets "at-scale" parameters
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
glmeiv_nf_pipeline=$LOCAL_CODE_DIR"glmeiv-pipeline/main.nf"
thresholding_nf_pipeline=$PWD"/../nextflow_scripts/thresholding.nf"

######################
# 2. Set all arguments
######################
work_dir=$LOCAL_GLMEIV_DATA_DIR"work"
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
  gRNA_gene_pairs=$processed_data_dir"gene_gRNA_pairs_sample_problem.rds"
  gene_pod_size=5 #500
  gRNA_pod_size=5 #500
  pair_pod_size=5 #500
else
  gRNA_gene_pairs=$processed_data_dir"gRNA_gene_pairs.rds"
  gene_pod_size=200
  gRNA_pod_size=500
  pair_pod_size=250
fi

# vii. results directory
result_dir=$LOCAL_GLMEIV_DATA_DIR"public/gasperini/results"
result_file_name="glmeiv_result.rds"
