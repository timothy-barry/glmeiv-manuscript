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
Rscript ../R_scripts/analysis/qc_gasperini.R
# run glmeiv at scale on Gasperini data
bash run_gasp_glmeiv.sh
# run thresholding method at scale on Gasperini data
bash run_gasp_thresh.sh

#################
# 2. Xie analysis
#################
# QC Xie data

# run glmeiv at scale on Xie data
Rscript ../R_scripts/analysis/qc_xie.R
Rscript run_xie_glmeiv.sh

################
# 3. Simulations
################
Rscript ../R_scripts/simulations/create_simspec_objects.R
$SIMULATR -f $sim_dir"spec_objects/sim_spec_0.rds" -r $sim_dir"results/raw_result_1.rds" -b 2
