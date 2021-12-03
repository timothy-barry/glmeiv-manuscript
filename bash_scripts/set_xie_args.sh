source ~/.research_config

# use -t flag to set "trial" parameters; otherwise, sets "at-scale" parameters
trial=false
azure=false
while getopts ta OPT
do
    case "$OPT" in
        t) trial=true ;;
        a) azure=true ;;
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
# i. Local vs Azure
if [ $azure = true ]
then
  glmeiv_dir="az://main/glmeiv/"
  xie_2019_dir="az://main/xie-2019/"
  work_dir="az://main/work"
  prof="-profile az"
else
  glmeiv_dir=$LOCAL_GLMEIV_DATA_DIR
  xie_2019_dir=$LOCAL_XIE_2019_DATA_DIR
  work_dir=$LOCAL_GLMEIV_DATA_DIR"work"
  prof=""
fi
backing_files_dir=$xie_2019_dir"processed/"
processed_data_dir=$glmeiv_dir"public/xie/data/"

# ii. Location of backing .odm files
gene_odm=$backing_files_dir"gene/expression_matrix.odm"
gRNA_odm=$backing_files_dir"gRNA/raw_grouped.odm"
# iii. metadata files
gene_metadata=$processed_data_dir"gene_metadata.rds"
gRNA_metadata=$processed_data_dir"gRNA_metadata.rds"
# iv. covariate matrix and offsets
covariate_matrix=$processed_data_dir"covariate_matrix.rds"
m_offsets=$processed_data_dir"m_offset.rds"
g_offsets=$processed_data_dir"g_offset.rds"
# v. family strings
m_fam="nb"
g_fam="poisson"
# vi. pairs to analyze and pod sizes
if [ $trial = true ]
then
  gRNA_gene_pairs=$processed_data_dir"gRNA_gene_pairs_problem.rds"
  gene_pod_size=3 #500
  gRNA_pod_size=3 #500
  pair_pod_size=3 #500
  result_file_name="glmeiv_result_trial.rds"
else
  gRNA_gene_pairs=$processed_data_dir"gRNA_gene_pairs.rds"
  gene_pod_size=200
  gRNA_pod_size=500
  pair_pod_size=200
  result_file_name="glmeiv_result.rds"
fi

# vii. results directory
result_dir=$LOCAL_GLMEIV_DATA_DIR"public/xie/results"
