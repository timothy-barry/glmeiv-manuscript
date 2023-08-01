# source file paths
source ~/.research_config
sim_dir=$LOCAL_GLMEIV_DATA_DIR"public/simulations/"
result_dir=$sim_dir"results/"
spec_obj_3=$sim_dir"spec_objects/sim_spec_3.rds"

nextflow pull timothy-barry/simulatr-pipeline
nextflow run timothy-barry/simulatr-pipeline \
 --simulatr_specifier_fp $spec_obj_3 \
 --result_dir $result_dir \
 --result_file_name "sim_res_3.rds" \
 --max_gb 3
 