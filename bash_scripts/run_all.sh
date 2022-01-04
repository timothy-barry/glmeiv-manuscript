# Runs the entire analysis from start to finish.
# It is assumed that this script is being executed from the "bash_scripts" subdirectory.

##########
# 0. setup
##########
# source file paths
bash setup.sh
source ~/.research_config

################
# 1. QC datasets
################
# QC Gasperini data
Rscript ../R_scripts/processing/qc_gasperini.R
# QC Xie data
Rscript ../R_scripts/processing/qc_xie.R

#####################
# 2. GLM-EIV analysis
#####################
# run glmeiv at scale on Gasperini data
bash run_gasp_glmeiv.sh
# run glmeiv at scale on Xie data
bash run_xie_glmeiv.sh

##########################
# 3. Thresholding analysis
##########################
# Get the optimal threshold from the GLM-EIV results
Rscript ../R_scripts/analysis/approx_bayes_bdy.R
# run thresholding on gasperini data
bash run_gasp_thresholding.sh
# run thresholding on xie data
bash run_xie_thresh.sh
# run thresholding on positive control Gasperini pairs, setting thresh to 1,5,20.
bash run_gasp_pc_thresholding.sh

################
# 4. Simulations
################
sim_dir=$LOCAL_GLMEIV_DATA_DIR"public/simulations/"
Rscript ../R_scripts/simulations/create_simspec_objects.R
$SIMULATR -f $sim_dir"spec_objects/sim_spec_0.rds" -r $sim_dir"results/raw_result_0.rds"
$SIMULATR -f $sim_dir"spec_objects/sim_spec_1.rds" -r $sim_dir"results/raw_result_1.rds"
$SIMULATR -f $sim_dir"spec_objects/sim_spec_2.rds" -r $sim_dir"results/raw_result_2.rds"

##################################
# 5. Gasperini resampling analysis
##################################
bash run_gasp_resampling.sh

#############
# 6. Plotting
#############
Rscript ../R_scripts/plotting/analysis_challenges.R
Rscript ../R_scripts/plotting/data_analysis.R
Rscript ../R_scripts/plotting/main_text_sim.R
Rscript ../R_scripts/plotting/thresholding_empirical.R
Rscript ../R_scripts/plotting/thresholding_theoretical.R
Rscript ../R_scripts/plotting/supplement_sim.R
