# Code for replicating "Exponential family measurement error models for single-cell CRISPR screens"

This repository contains code to reproduce the analyses reported in the paper "Exponential family measurement error models for single-cell CRISPR screens" by T Barry, E Katsevich, and K Roeder (2022). The paper introduces the "GLM-EIV" (generalized linear model-based errors-in-variables) model and associated method. The `glmeiv` package is available [here](https://github.com/timothy-barry/glmeiv), and the at-scale GLM-EIV pipeline is available [here](https://github.com/timothy-barry/glmeiv-pipeline).

# 1. Dependencies

The code in this repository has several dependences.

1. **Nextflow**. `Nextflow` is a programming language that faciliates building data-intensive pipelines for deployment on HPC and cloud. Download and install `Nextflow`using the instructions [here](https://github.com/timothy-barry/glmeiv).

2. **R language**.

3. **The following `R` packages**.
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

4. **simulatr command line utility**. `simulatr` is a (currently in-house) R package and associated command-line utility for running simulations on HPC and cloud. Install the `simulatr` command line utility via the following command:
```
git clone https://github.com/timothy-barry/simulatr-command-line.git
```

# 2. Data

Download the the following three data directories from Box: https://upenn.box.com/s/34bi9u90sgsajbvwuacqcd8pog1twq6z, https://upenn.box.com/s/qgbr6dmhr1dkviupibu7ef7bymde37of, https://upenn.box.com/s/9gv48i342liv7zz8ykedch2m2t9ogisw. The first two directories contain the Gasperini 2019 and Xie 2019 data, and the third contains some metadata specific to the GLM-EIV paper. The data are single-cell expression matrices stored as [ondisc](https://github.com/timothy-barry/ondisc) objects.
