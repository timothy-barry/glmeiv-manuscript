# get the input simulation id
sim_no=$1

# source file paths
source ~/.research_config
sim_dir=$LOCAL_GLMEIV_DATA_DIR"public/simulations/"
result_dir=$sim_dir"results/"
spec_obj=$sim_dir"spec_objects/sim_spec_$sim_no.rds"

# ensure all up to date
Rscript ../../../R_scripts/simulations/create_simspec_objects.R
nextflow pull timothy-barry/simulatr-pipeline

# run
nextflow run timothy-barry/simulatr-pipeline \
 --simulatr_specifier_fp $spec_obj \
 --result_dir $result_dir \
 --result_file_name "sim_res_$sim_no.rds" \
 --max_gb 3
 