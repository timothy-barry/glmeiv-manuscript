# Code for replicating "Exponential family measurement error models for single-cell CRISPR screens"

This repository contains code to reproduce the analyses reported in the paper "Exponential family measurement error models for single-cell CRISPR screens" by Tim Barry, Eugene Katsevich, and Kathryn Roeder (2022). The paper introduces the "GLM-EIV" (generalized linear model-based errors-in-variables) model and associated method. The `glmeiv` package is available [here](https://github.com/timothy-barry/glmeiv), and the at-scale GLM-EIV pipeline is available [here](https://github.com/timothy-barry/glmeiv-pipeline).

# 1. Dependencies

The code in this repository has several dependences.

a. **Nextflow**. `Nextflow` is a programming language that faciliates building data-intensive pipelines for deployment on HPC and cloud. Download and install `Nextflow`using the instructions [here](https://github.com/timothy-barry/glmeiv).

b. **R language**.

c. **The following `R` packages**.
```
# Our packages
install.packages("devtools")
devtools::install_github("timothy-barry/ondisc")
devtools::install_github("timothy-barry/simulatr")
devtools::install_github("timothy-barry/glmeiv")

# Other peoples' packages
install.packages("tidyverse")
install.packages("cowplot")
install.packages("MASS")
```

d. **simulatr command line utility**. `simulatr` is a (currently in-house) R package and associated command-line utility for running simulations on HPC and cloud. Install the `simulatr` command line utility via the following command:
```
git clone https://github.com/timothy-barry/simulatr-command-line.git
```

# 2. Data

Download the the following three data directories from Box: https://upenn.box.com/s/34bi9u90sgsajbvwuacqcd8pog1twq6z (Gasperini data), https://upenn.box.com/s/qgbr6dmhr1dkviupibu7ef7bymde37of (Xie data), https://upenn.box.com/s/9gv48i342liv7zz8ykedch2m2t9ogisw (GLM-EIV results and metadata). The first two directories contain the Gasperini 2019 and Xie 2019 data, respectively, and the third contains some metadata specific to the current project. The large files in the first two directories are single-cell expression matrices stored as [ondisc](https://github.com/timothy-barry/ondisc) objects.

# 3. Config file

We use config files to increase code portability. You will need to create a config file of your own to replicate our analyses. Use the following steps.

i. Open the terminal.

ii. `cd` to the home directory and create a file called `.research_config`:
```
cd; touch .research_config
```
iii. Edit the .research_config file. Define the variables `LOCAL_GLMEIV_DATA_DIR`, `LOCAL_GASPERINI_2019_DATA_DIR`, `LOCAL_XIE_2019_DATA_DIR`, and `SIMULATR` to be the file paths to the GLM-EIV results directory, Gasperini data directory, Xie data directory, and simulatr directory on your machine. For example, I (Tim) have the above directories stored in the following locations on my machine: "/Users/timbarry/research_offsite/glmeiv", "/Users/timbarry/research_offsite/gasperini-2019", "/Users/timbarry/research_offsite/xie-2019", and "/Users/timbarry/research_code/simulatr-command-line", respectively. Therefore, my config file is as follows:  

```
LOCAL_GLMEIV_DATA_DIR="/Users/timbarry/research_offsite/glmeiv"
LOCAL_GASPERINI_2019_DATA_DIR="/Users/timbarry/research_offsite/gasperini-2019"
LOCAL_XIE_2019_DATA_DIR="/Users/timbarry/research_offsite/xie-2019"
SIMULATR="/Users/timbarry/research_code/simulatr-command-line"
```
