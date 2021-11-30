# Runs the entire analysis from start to finish.
# It is assumed that this script is being executed from the "bash_scripts" subdirectory.

##########
# 0. setup
##########
# source file paths
bash setup.sh
source ~/.research_config

# set a couple variables
sim_dir=$LOCAL_GLMEIV_DATA_DIR"public/simulations/"

#######################
# 1. Gasperini analysis
#######################
# QC Gasperini data
Rscript ../R_scripts/processing/qc_gasperini.R
# run glmeiv at scale on Gasperini data
bash run_gasp_glmeiv.sh
# run thresholding method at scale on Gasperini data
bash run_gasp_thresh.sh

#################
# 2. Xie analysis
#################
# QC Xie data
Rscript ../R_scripts/processing/qc_xie.R
# run glmeiv at scale on Xie data
bash run_xie_glmeiv.sh
# run thresholding method at scale on Xie data
bash run_xie_thresh.sh

################
# 3. Simulations
################
Rscript ../R_scripts/simulations/create_simspec_objects.R
$SIMULATR -f $sim_dir"spec_objects/sim_spec_0.rds" -r $sim_dir"results/raw_result_0.rds"
$SIMULATR -f $sim_dir"spec_objects/sim_spec_1.rds" -r $sim_dir"results/raw_result_1.rds"
$SIMULATR -f $sim_dir"spec_objects/sim_spec_2.rds" -r $sim_dir"results/raw_result_2.rds"

#############
# 4. Plotting
#############
