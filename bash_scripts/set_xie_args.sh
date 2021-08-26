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
thresholding_nf_pipeline=$PWD"/../nextflow_scripts/thresholding.nf"

######################
# 2. Set all arguments
######################
# i. Location of backing .odm files
backing_files_dir=$LOCAL_XIE_2019_DATA_DIR"processed/"
gene_odm=$backing_files_dir"gene/expression_matrix.odm"
gRNA_odm=$backing_files_dir"gRNA/raw_grouped.odm"
# ii. Location of processed data
processed_data_dir=$LOCAL_GLMEIV_DATA_DIR"public/xie/data/"
# iii. metadata files
gene_metadata=$processed_data_dir"gene_metadata.rds"
gRNA_metadata=$processed_data_dir"gRNA_metadata.rds"
# iv. covariate matrix and offsets
covariate_matrix=$processed_data_dir"covariate_matrix.rds"
m_offsets=$processed_data_dir"m_offset.rds"
g_offsets=$processed_data_dir"g_offset.rds"
# v. family strings
m_fam="nb"
g_fam="nb"
g_theta=10
# vi. pairs to analyze and pod sizes
if [ $trial = true ]
then
  gRNA_gene_pairs=$processed_data_dir"gRNA_gene_pairs_trial.rds"
  gene_pod_size=3 #500
  gRNA_pod_size=3 #500
  pair_pod_size=3 #500
else
  gRNA_gene_pairs=$processed_data_dir"gRNA_gene_pairs.rds"
  gene_pod_size=200
  gRNA_pod_size=500
  pair_pod_size=500
fi

# vii. results directory
result_dir=$LOCAL_GLMEIV_DATA_DIR"public/xie/results"
result_file_name="glmeiv_result_nb10.rds"
