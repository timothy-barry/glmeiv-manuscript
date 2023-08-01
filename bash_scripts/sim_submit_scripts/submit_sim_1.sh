# source file paths
source ~/.research_config
sim_dir=$LOCAL_GLMEIV_DATA_DIR"public/simulations/"
result_dir=$sim_dir"results/"
spec_obj_1=$sim_dir"spec_objects/sim_spec_1.rds"

nextflow pull timothy-barry/simulatr-pipeline
# nextflow run timothy-barry/simulatr-pipeline \
nextflow run /Users/timbarry/research_code/simulatr-pipeline/main.nf \
 --simulatr_specifier_fp $spec_obj_1 \
 --result_dir $result_dir \
 --result_file_name "sim_res_1.rds" \
 --max_gb 3
 