# Code for _Exponential family measurement error models for single-cell CRISPR screens_

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
