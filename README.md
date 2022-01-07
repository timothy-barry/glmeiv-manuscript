# Code for replicating "Exponential family measurement error models for single-cell CRISPR screens"

This repository contains code to reproduce the analyses reported in the paper "Exponential family measurement error models for single-cell CRISPR screens" by T Barry, E Katsevich, and K Roeder (2022). This paper introduces the "GLM-EIV" (GLM-based errors-in-variables) method and applies the method to analyze two large-scale single-cell CRISPR screen datasets. The `glmeiv` package is available [here](https://github.com/timothy-barry/glmeiv), and the at-scale GLM-EIV pipeline is available [here](https://github.com/timothy-barry/glmeiv-pipeline).

# 1. Dependencies

The code in this repository has several dependences.

a. **Nextflow**. `Nextflow` is a programming language that facilitates building data-intensive pipelines for deployment on HPC and cloud. Download and install `Nextflow`using the instructions [here](https://www.nextflow.io).

b. **R language**.

c. **The following `R` packages**.
```
# Other developers' packages
install.packages("tidyverse")
install.packages("cowplot")
install.packages("MASS")
install.packages("devtools")

# Our packages
devtools::install_github("timothy-barry/ondisc")
devtools::install_github("timothy-barry/simulatr")
devtools::install_github("timothy-barry/glmeiv")
```

d. **simulatr command line utility**. `simulatr` is a (currently in-house) R package and associated command-line utility for running simulation experiments on HPC and cloud. Install the `simulatr` command line utility via the following command:
```
git clone https://github.com/timothy-barry/simulatr-command-line.git
```

# 2. Data

Download the the following three data directories from Box: https://upenn.box.com/v/gasperini-data-2019 (Gasperini data), https://upenn.box.com/v/xie-data-2019 (Xie data), https://upenn.box.com/v/glmeiv-files-v1 (GLM-EIV results and metadata). The first two directories contain the Gasperini 2019 and Xie 2019 data, respectively, and the third contains metadata, results, and other files specific to the GLM-EIV project.

The large files in the first two directories are single-cell expression matrices stored as [ondisc](https://github.com/timothy-barry/ondisc) objects. `ondisc` is an R package for large-scale single-cell computing.

# 3. Config file

We use a config file in our projects to increase code portability. You will need to create a config file of your own to replicate our analyses. Use the following steps.

i. Open the terminal.

ii. `cd` to the home directory and create a file called `.research_config`:
```
cd; touch .research_config
```
iii. Edit the .research_config file. Set the variables `LOCAL_GLMEIV_DATA_DIR`, `LOCAL_GASPERINI_2019_DATA_DIR`, `LOCAL_XIE_2019_DATA_DIR`, and `SIMULATR` to the file paths of the GLM-EIV results directory, Gasperini data directory, Xie data directory, and simulatr directory, respectively, on your machine. For example, I (TB) have the above directories stored in the following locations on my machine: "/Users/timbarry/research_offsite/glmeiv", "/Users/timbarry/research_offsite/gasperini-2019", "/Users/timbarry/research_offsite/xie-2019", and "/Users/timbarry/research_code/simulatr-command-line". Therefore, my config file is as follows:  

```
LOCAL_GLMEIV_DATA_DIR="/Users/timbarry/research_offsite/glmeiv"
LOCAL_GASPERINI_2019_DATA_DIR="/Users/timbarry/research_offsite/gasperini-2019"
LOCAL_XIE_2019_DATA_DIR="/Users/timbarry/research_offsite/xie-2019"
SIMULATR="/Users/timbarry/research_code/simulatr-command-line"
```

Edit the config file above so that the variables point to the correct locations on your machine.

# 4. Manuscript analysis code
Clone the `glmeiv-manuscript` repository (which contains all the analysis replication code) to your machine:
```
git clone https://github.com/timothy-barry/glmeiv-manuscript.git
```
# 5. Replicate

We are now ready to replicate the analayses! Change directories to the `bash_scripts` subdirectory of the `glmeiv-manuscript` directory.
```
cd glmeiv-manuscript/bash_scripts
```
The script `run_all.sh` is the master bash script that replicates all analayses reported in the paper. In theory you could run this script by calling `bash run_all.sh`. However, executing this script line-by-line is a better option, as some commands within this script take a month or longer to execute on a laptop. Note that the data directories (downloaded in step 2) contain all intermediate and results files required to reproduce the analyses. Therefore, you can execute _any_ line of the run_all.sh script in any order.

# 6. A note on Nextflow
The master `run_all.sh` script calls a sequence of bash scripts, each available in the `glmeiv-manuscript/bash_scripts` directory. Several of these bash scripts -- namely, `bash run_gasp_gmeiv.sh`, `bash run_xie_glmeiv.sh`, `bash run_gasp_thresholding.sh`, `bash run_xie_thresh.sh`, `bash run_gasp_pc_thresholding.sh`, and `bash run_gasp_resampling.sh` -- call Nextflow pipelines. You can execute all Nextflow-calling bash scripts on your local machine, provided that you have downloaded and installed the Nextflow executable.

We note that two of the Nextflow pipelines -- namely, those contained in `bash run_gasp_gmeiv.sh` and `bash run_xie_glmeiv.sh` -- may take weeks or longer to run on a laptop. You instead can execute these pipelines on an HPC or cloud to reduce the compute time to a few hours. The details of how to use Nextflow are beyond the scope of this brief note (see the [docs](https://www.nextflow.io/docs/latest/index.html)). However, if you are working on an HPC with a SLURM scheduler (a common choice on many university HPCs), you can use the following quick approach.

1. Login to your HPC. Install the dependences and data and set up your config file.
2. Create a file called `~/.nextflow/config` with the following contents:

```
process {
   executor = 'slurm'
}

executor {
   queueSize = 200
   submitRateLimit = '30 sec'
   pollInterval = '1 min'
   queueStatInterval = '5 min'
}
```
3. Execute the bash script by submitting it as a job to the scheduler using `sbatch`. For example, `sbatch run_gasp_gmeiv.sh`.

If you are interested in running the pipeline on an HPC that uses a different scheduler (e.g., Sun Grid) or on a cloud (e.g., Azure or AWS), please contact the authors after briefly familiarizing yourself with the Nextflow docs.
