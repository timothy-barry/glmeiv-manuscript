# source file paths
bash setup.sh
source ~/.research_config
sim_dir=$LOCAL_GLMEIV_DATA_DIR"public/simulations/"
result_dir=$sim_dir"results/"

# sim 0
spec_obj_0=$sim_dir"spec_objects/sim_spec_0.rds"
nextflow pull timothy-barry/simulatr-pipeline
nextflow run timothy-barry/simulatr-pipeline \
 --simulatr_specifier_fp $spec_obj_0 \
 --result_dir $result_dir \
 --result_file_name "sim_res_0.rds" \
 --B 5
