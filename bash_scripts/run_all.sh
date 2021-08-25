# Runs the entire analysis from start to finish.
# It is assumed that this script is being executed from the "bash_scripts" subdirectory.

##########
# 0. setup
##########
bash setup.sh
source ~/.research_config

#######################
# 1. Gasperini analysis
#######################
# QC Gasperini data
Rscript ../R_scripts/qc_gasperini.R
# run glmeiv at scale on Gasperini data
bash run_gasp_glmeiv.sh
# run thresholding method at scale on Gasperini data
bash run_gasp_thresh.sh

#################
# 2. Xie analysis
#################
# QC Xie data

# run glmeiv at scale on Xie data
Rscript ../R_scripts/qc_xie.R
Rscript run_xie_glmeiv.sh
