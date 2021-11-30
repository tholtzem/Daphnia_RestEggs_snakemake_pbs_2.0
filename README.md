# Daphnia_RestEggs_snakemake_pbs

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.6.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

======================================================

Snakemake workflow for the analysis of DNA from *Daphnia*-resting eggs on the mach2 HPC cluster with the PBS-torque batch job submission system. 

The snakemake workflow was modified from the [daphnia_snakemake_pbs Tutorial](https://github.com/tholtzem/daphnia_snakemake_pbs), which was initially based on the [ta_dna_snakemake_pbs Tutorial](https://github.com/schimar/ta_dna_snakemake_pbs),the [Snakemake Cluster Tutorial](https://github.com/SchlossLab/snakemake_cluster_tutorial.git) and the [Software Carpentry lesson repository](https://hpc-carpentry.github.io/hpc-python/17-cluster/). For more information on snakemake itself (https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).

The analysis is based on the [physalia-lcwgs: the Physalia course on Population genomic inference from low-coverage whole-genome sequencing data, Oct 19-22 2020](https://github.com/nt246/physalia-lcwgs).

======================================================

## conda and other [dependencies](https://github.com/schimar/ta_dna_snakemake_pbs/blob/main/envs/s21.yaml)   

create environment from yaml file (in envs/):
```
# run these two once, to create the environment:
conda init bash
conda env create -f envs/s21.yaml

# with this, you can activate the environment with all [dependencies](https://github.com/schimar/ta_dna_snakemake_pbs/blob/main/envs/s21.yaml):
conda activate eggs

# (also, when ssh'ing onto mach2, you can activate the env and then do a dry-run of your workflow) 
## how to submit the main snakemake job:
qsub code/clusterSnakemake.pbs

# if you've added new software to install to the conda environment, then you can update:
conda env update --name eggs --file envs/s21.yaml
```
## mamba

Mamba (https://github.com/mamba-org/mamba) is a reimplementation of the conda package manager in C++.

```
# First, remove conda environment
conda env remove -n eggs

# Load Anaconda on cluster (here mach2):
module load Anaconda3/2021.04/miniconda-base-2021.04

# To use conda commands in your current shell session, first do:
source $UIBK_CONDA_PROFILE

# Create environment from yaml file (in envs/):
conda init bash
mamba env create -f envs/s21.yaml

# Activate the environment
conda activate eggs

# if you've added new software to install to the conda environment, then you can update:
mamba env update --name eggs --file envs/s21.yaml


```

