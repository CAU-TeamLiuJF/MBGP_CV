mbBayesABLD (Bayesian multi-breed genomic prediction model that partitions each chromosome into non-overlapping blocks based on LD patterns)
This repository contains the data, scripts, and program files required to reproduce the paper *"Multi-trait Bayesian models enhance the accuracy of genomic prediction in multi-breed reference populations"*.

> Please cite the original paper when using the contents of this repository.
> Author contact information is given at the end.
> 

# User Manual

- [1. Introduction](#1-introduction)
- [2. Runtime environment and dependencies](#2-runtime-environment-and-dependencies)
- [3. Obtaining the code and data](#3-obtaining-the-code-and-data)
  - [3.1 Download compressed package (offline or no git access)](#31-download-compressed-package-offline-or-no-git-access)
  - [3.2 Clone from GitHub (online)](#32-clone-from-github-online)
- [4. Directory structure (detailed)](#4-directory-structure-detailed)
  - [4.1 Directory Structure](#41-directory-structure)
  - [4.2 Key files explained](#42-key-files-explained)
- [5. Initialization (`initialize.sh`) details](#5-initialization-initializesh-details)
- [6. `main.sh` workflow detailed explanation](#6-mainsh-workflow-detailed-explanation)
  - [6.1 Top-level loop — per data source (`Xie2021`, `Keller2022`)](#61-top-level-loop--per-data-source-xie2021-keller2022)
  - [6.2 Random seed handling](#62-random-seed-handling)
  - [6.3 Phenotype adjustment (only for Xie2021)](#63-phenotype-adjustment-only-for-xie2021)
  - [6.4 Within-breed GBLUP (within-breed prediction)](#64-within-breed-gblup-within-breed-prediction)
  - [6.5 Multi-breed single-trait and multi-trait predictions](#65-multi-breed-single-trait-and-multi-trait-predictions)
  - [6.6 Bayesian within-breed prediction](#66-bayesian-within-breed-prediction)
  - [6.7 Summary statistics (accur, var, time)](#67-summary-statistics-accur-var-time)
- [7. Important scripts and usage](#7-important-scripts-and-usage)
  - [7.1 `GP_cross_validation.sh` (`$GP_cross`)](#71-gp_cross_validationsh-gp_cross)
  - [7.2 `GP_single_breed.sh` (`$GP_single`)](#72-gp_single_breedsh-gp_single)
  - [7.3 `GP_multi_breed.sh` (`$GP_multi`)](#73-gp_multi_breedsh-gp_multi)
  - [7.4 R scripts](#74-r-scripts)
- [8. Input file formats and examples](#8-input-file-formats-and-examples)
  - [8.1 Genotype (plink bfile)](#81-genotype-plink-bfile)
  - [8.2 Phenotype file `phenotype.txt`](#82-phenotype-file-phenotypetxt)
  - [8.3 `breeds.txt`](#83-breedstxt)
  - [8.4 `traits_all.txt` and `traits_analysis.txt`](#84-traits_alltxt-and-traits_analysistxt)
  - [8.5 `breeds_combination.txt`](#85-breeds_combinationtxt)
- [9. Output files description](#9-output-files-description)
  - [9.1 `phe_adj_BLUP.*` (phenotype adjustment results)](#91-phe_adj_blup-phenotype-adjustment-results)
  - [9.2 `accur_*` (accuracy summaries)](#92-accur_-accuracy-summaries)
  - [9.3 `var_*` (variance component summaries)](#93-var_-variance-component-summaries)
  - [9.4 `time_*` (runtime statistics)](#94-time_-runtime-statistics)
  - [9.5 Simulation outputs (Simulation section)](#95-simulation-outputs-simulation-section)
- [10. Common run scenarios and example commands](#10-common-run-scenarios-and-example-commands)
  - [10.1 With Slurm (recommended)](#101-with-slurm-recommended)
  - [10.2 Without Slurm (local testing)](#102-without-slurm-local-testing)
  - [10.3 Run a single step (examples)](#103-run-a-single-step-examples)
- [11. Performance advice and parameter tuning](#11-performance-advice-and-parameter-tuning)
- [12. Troubleshooting](#12-troubleshooting)
  - [12.1 Common errors and how to locate them](#121-common-errors-and-how-to-locate-them)
  - [12.2 Log file locations](#122-log-file-locations)
- [13. FAQ (frequently asked questions)](#13-faq-frequently-asked-questions)
- [14. License and citation](#14-license-and-citation)
- [15. Contact and acknowledgements](#15-contact-and-acknowledgements)
- [Appendix A: Example full run (step-by-step)](#appendix-a-example-full-run-step-by-step)
- [Appendix B: Practical tips](#appendix-b-practical-tips)
- [Closing remarks](#closing-remarks)


# 1. Introduction

This project provides a complete analysis pipeline to reproduce the multi-breed prediction experiments and simulation results from the paper.
 The workflow includes real data processing (phenotype adjustment, single-/multi-breed models, comparison of methods), simulation data generation (QMSim), and result aggregation (accuracy, variance components, runtime).

The repository organizes programs into layers: a top-level control script (`main.sh`), core functionality script (`shell/GP_cross_validation.sh`, abbreviated `GP_cross`), R scripts, and compiled binaries or executables under `bin/`.

 [back to top](#user-manual)

# 2. Runtime environment and dependencies

**Operating system**

- Recommended: CentOS Linux release 7.x (tested: 7.6.1810)
- Other mainstream Linux distributions can be used, but paths and package management may need manual adjustment.

**Job scheduler**

- Recommended: Slurm workload manager (scripts submit jobs using `sbatch`).
- If Slurm is not available, see the instructions to comment out `sbatch` lines and run commands directly.

**Software dependencies (install beforehand)**

- R (recommended ≥ 4.1.0; tested with R 4.3.1)
- plink (for genotype processing)
- QMSim (for simulations)
- GNU bash (scripts are bash-based)
- make / gcc (if compiling tools in `bin/`)
- Standard Unix tools: awk, sed, grep, (python optional)

**R packages (initialization script will check)**

- Common: `Cairo`, `coda`, `data.table`, `dplyr`, `getopt`, `ggplot2`, `lattice`, `LAVA`, `MASS`, `Matrix`, `pedigree`, `Rcpp`, `RcppArmadillo`, `RcppEigen`, `reshape`, `reshape2`, `stringr`, `tidyr`.
- If initialization reports missing packages, run `install.packages("packageName")` in R or use your organization’s CRAN mirror.

**Hardware recommendations**

- High-performance nodes are recommended (≥ 50 CPU cores, ≥ 100 GB memory) for full runs.
- If resources are limited, follow the instructions to reduce thread counts and parallel job submissions.

 [back to top](#user-manual)

# 3. Obtaining the code and data

## 3.1 Download compressed package (offline or no git access)

Assuming you have the zip file on the server under `/public/home/liujf/liwn/download`, run:

```bash
cd /public/home/liujf/liwn/download
unzip -d /public/home/liujf/liwn/code/GitHub mbBayesABLD-main.zip
mv /public/home/liujf/liwn/code/GitHub/mbBayesABLD-main /public/home/liujf/liwn/code/GitHub/mbBayesABLD
```

## 3.2 Clone from GitHub (online)

If the server has Internet access, clone the repository:

```bash
cd /public/home/liujf/liwn/code/GitHub
git clone git@github.com:CAU-TeamLiuJF/mbBayesABLD.git
```

After cloning, verify that expected subdirectories (e.g., `data`, `bin`, `shell`) are present.

 [back to top](#user-manual)

# 4. Directory structure (detailed)

## 4.1 Directory Structure

Below is an example top-level layout and purpose of each folder/file:

```
mbBayesABLD/
├─ bin/                    # Executables and binary tools (ensure executable permission)
├─ code/                   # R language scripts for tasks such as phenotype simulation
├─ prm/                    # Template file for generating QMSim software parameter cards
├─ data/
│  ├─ Real/
│  │  ├─ Xie2021/          # Pig dataset (genotype, phenotype, configs)
│  │  └─ Keller2022/       # Bean dataset
│  └─ Simulation/          # Random number seeds used to reproduce simulated dataset
├─ shell/
│  ├─ GP_cross_validation.sh  # Core pipeline script (GP_cross)
│  └─ ...                     # Other helper shell scripts
├─ initialize.sh           # Initialization script: set PATH, check R packages
├─ main.sh                 # Top-level controller script orchestrating the workflow
└─ README.md               # This detailed README (plain-text)
```

## 4.2 Key files explained

- This project contains several important scripts and directories. Below are the key files and their purposes:
- **main.sh**
  The central pipeline script to run genomic prediction with cross-validation. It integrates data preprocessing, model execution, and result collection. Users typically start analysis by executing this script.
- **initialize.sh**
  A script to initialize the working environment, check required dependencies, and prepare necessary directories before running the main workflow.
- **QMSim.sh**
  A helper script to execute QMSim simulations. It automates simulation runs using parameter templates under the `prm/` directory.
- **R/accuracy_bias_calculation.R**
  Calculates prediction accuracy and potential bias in genomic prediction results. This is a key evaluation script after cross-validation.
- **R/multibreed_relationship_matrix.R**
  Builds relationship matrices across multiple breeds, an essential step in multi-breed genomic prediction scenarios.
- **R/pheno_simulation.R**
  Simulates phenotypic values under different genetic models and breeding structures. Often used in conjunction with QMSim for simulation studies.
- **shell/GP_cross_validation.sh**
  Implements genomic prediction with cross-validation in a shell environment. It is called by `main.sh` to conduct large-scale model training and validation.
- **shell/GP_multi_breed.sh**
  Focusing on genomic prediction with multi-breed datasets, we ran the m-STGBLUP, m-BayesR, m-MTGBLUP, mbBayesAB-fix, mbBayesAB-lava, and mbBayesABLD models.
- **shell/GP_single_breed.sh**
  Similar to the above but designed for single-breed prediction to run the w-STGBLUP and w-BayesABLD models.
- **shell/BayesR.sh**
  Executes BayesR analyses with prepared input files. Closely linked to R scripts for file formatting and post-analysis visualization.
- **bin/**
  Contains compiled executables such as `dmu4`, `QMSim`, and `mbBayesABLD`, which are the core computational engines for genetic parameter estimation, simulation, and Bayesian prediction.
- **data/**
  Includes both real datasets (e.g., Keller2022, Xie2021) and simulation datasets (TwoBreeds, ThreeBreeds). These serve as input examples for testing and reproduction of the pipeline.
- **prm/**
  Parameter template files for QMSim runs, including genome definition, population history, and subpopulation setups.

 [back to top](#user-manual)

# 5. Initialization (`initialize.sh`) details

The script `initialize.sh` is provided to set up the environment and paths for the project. It should be executed once after downloading or cloning the repository, or whenever the project folder is moved to a different location.

Run the script from the project root:

```
cd /public/home/liujf/liwn/code/GitHub/mbBayesABLD
./initialize.sh
```

This script performs the following tasks:

- **Locate project path**
   If the variable `code` is already defined, the script checks whether this folder exists. Otherwise, it automatically determines the path of the script using `readlink -f`.

- **Switch to main directory**
   Changes the current working directory to the main project folder where `main.sh` is located.

- **Update `main.sh`**
   Replaces the `code=` line inside `main.sh` so that it points to the correct project directory.
   Updates the reference to `initialize.sh` inside `main.sh` for clarity and reproducibility.

- **Load custom shell functions**
   Sources the script `shell/function.sh`.
   If this file is missing, the script will exit with an error.

- **Check R installation**
   Uses the function `check_command Rscript` to verify that `Rscript` can be found in the environment.
   If R is not available, it prints an error message and exits.
   On HPC systems, you may need to load R with:

  ```
  module load R
  ```

- **Determine R installation path**
   Saves the absolute path of `Rscript` (e.g. `/usr/bin/Rscript` or `/opt/R/4.3.1/bin/Rscript`).

- **Standardize R script headers**
   For every `*.R` script in the `R/` folder:

  - If the first line already begins with `#!...Rscript`, it will be replaced with the correct `Rscript` path.
  - If no shebang line exists, one will be inserted automatically at the top.
     This ensures that all R scripts in the project can be executed directly as standalone programs.

- **Assign executable permissions**
   Grants user execute permissions (`chmod u+x`) for the following scripts:

  - `main.sh`
  - All scripts in `bin/`, `R/`, and `shell/` directories

- **Check required R packages**
   Looks for the file `R/package_required.R`.
   If found, this script will be executed to verify whether all necessary R packages are installed.
   If missing, the initialization stops with an error.

- **Create log folder**
   Creates a folder named `log/` in the main directory to store log files generated during later runs.

- **Final message**
   If all steps succeed, the script prints:

  ```
  Initialization completed.
  ```

 [back to top](#user-manual)

# 6. `main.sh` workflow detailed explanation

`main.sh` is divided into two main parts: **real data processing** and **simulations (QMSim)**. The following explains the key logic and parameters.

## 6.1 Top-level loop — per data source (`Xie2021`, `Keller2022`)

Main variables:

- `code`: repository root path (e.g., `/public/home/liujf/liwn/code/GitHub/MBGP_CV`)
- `GP_cross`: path to `shell/GP_cross_validation.sh`
- `pro`: project output path (`$code/output/${source}`)
- `data_path`: source data path (`$code/data/Real/${source}`)
- `bfile`: plink bfile prefix (default `genotype`)
- `phef`: phenotype file (default `phenotype.txt`)
- `breeds`: read from `breeds.txt` (space-separated)
- `traits_all`: read from `traits_all.txt` (all traits for multi-trait models)
- `traits_cal`: read from `traits_analysis.txt` (traits to analyze)

## 6.2 Random seed handling

Script checks `data/Real/xxx/random.seed`. If present, reuse it; otherwise generate `$RANDOM` and write to this file to ensure reproducibility.

## 6.3 Phenotype adjustment (only for Xie2021)

For `Xie2021` the pipeline runs phenotype correction (`--type adj`). Jobs are submitted via `sbatch`. The script waits until `phe_adj_BLUP.SOL` exists and has a sufficient number of lines before continuing.

For the dataset in `Keller2022`, the authors have already corrected for non-genetic effects, so it is not necessary to recalculate the corrected phenotypes.

## 6.4 Within-breed GBLUP (within-breed prediction)

For each trait in `traits_cal` and each breed in `breeds`, submit `$GP_cross --type within` with parameters:

- `--thread` (default 50; adjust to available cores)  

  Number of tasks submitted during parallel execution of breeding value estimation program

- `--rep` (default 10)

  The number of repetitions during cross validation

- `--fold` (default 5)

  Number of subsets divided during cross validation

## 6.5 Multi-breed single-trait and multi-trait predictions

`breeds_combination.txt` contains breed combinations (each line a combo). For each combination and method (`GBLUP`, `Bayes`, etc.), the script submits jobs:

- `single` (single-trait methods: GBLUP, bayesR)
- `multi` (multi-trait methods: GBLUP, mbBayesAB)

For `mbBayesAB`, the script loops over `bin` choices (`fix`, `lava`, `cubic`), submitting corresponding jobs.

## 6.6 Bayesian within-breed prediction

When block files (e.g., `cubic_*.txt`) generated, for each breed it submits `mbBayesABLD` within-breed prediction jobs.

## 6.7 Summary statistics (accur, var, time)

Finally, the script calls `$GP_cross --type accur|var|time` to aggregate results and generate summary files (`accur_*`, `var_*`, `time_*`).

 [back to top](#user-manual)

# 7. Important scripts and usage

## 7.1 `GP_cross_validation.sh` (`$GP_cross`)

**Function:** Provides multiple working modes via the `--type` parameter. Common `--type` options:

- `adj`: phenotype adjustment (produces BLUP/SOL and `phe_adj_BLUP.*`)
- `within`: within-breed prediction (GBLUP or mbBayesAB)
- `single`: multi-breed single-trait prediction
- `multi`: multi-breed multi-trait prediction (GBLUP or mbBayesAB)
- `psim`: simulate phenotypes (QMSim or internal)
- `geno`: genotype filtering/processing (plink-based)
- `bin`: genome partitioning / block file generation
- `var`: summarize variance components (heritability, SD-A, SD-P, etc.)
- `accur`: compute accuracies (correlation, RMSE, bias, etc.)
- `time`: collect runtime statistics

**Common parameters (non-exhaustive):**

- `--proj`: output path
- `--breeds`: single or multiple breeds (space-separated)
- `--traits`: list of all traits
- `--trait`: target trait for the current job
- `--bfile`: plink bfile prefix
- `--phef`: phenotype file path
- `--seed`: random seed
- `--thread`: thread count
- `--rep`, `--fold`: cross-validation settings
- `--method`: GBLUP, bayesR, mbBayesAB, etc.
- `--bin`: block strategy (`fix`/`lava`/`cubic`)
- `--burnin`, `--iter`: Bayesian model parameters

**Example (direct execution):**

```
$GP_cross --type within --proj /path/to/output --breeds YY --traits "MS PFAI" --trait PFAI --bfile /path/to/genotype --phef /path/to/phenotype.txt --seed 12345 --thread 20 --rep 10 --fold 5
```

## 7.2 `GP_single_breed.sh` (`$GP_single`)

**Function:** Run within-breed genomic prediction and cross-validation workflows for a single breed. This script prepares model inputs (phenotype, pedigree if provided, genotype), calls the chosen prediction engine (GBLUP/mbBayesAB/bayesR), performs `rep × fold` cross-validation, and writes per-fold/per-rep accuracy outputs.

**Core behaviors**

- Accepts a phenotype file and genotype (plink) prefix, partitions data into training/validation folds according to `--rep` and `--fold`, and runs the specified prediction method across folds.
- Supports using adjusted phenotypes / true breeding values (`--tbvf`) when available (for simulation or when using `adj` outputs).
- Integrates with DMU-based tools (`run_dmu4` / `run_dmuai`) or other executables in `bin/` to estimate GEBVs.
- Optionally runs Bayesian methods (mbBayesABLD) when provided with block files or partitioning (`--binf`, `--bin`) and MCMC parameters (`--iter`, `--burnin`).
- Writes per-breed accuracy outputs (e.g., `accur_<method>.txt`) and stores logs in the working folder for the breed.

**Common parameters (non-exhaustive)**

- `--label` : breed label (used to name output files / folders).
- `--phef` : phenotype file for this breed (format: ID + phenotype columns).
- `--bfile` : plink prefix (bed/bim/fam) used for genotype input.
- `--method` : prediction method (`GBLUP`, `bayesR`, `mbBayesAB`, etc.).
- `--tbvf` : file with "true" breeding values or adjusted phenotypes (e.g., `phe_adj_BLUP.txt`) to compute accuracy against.
- `--phereal` : index of the phenotype column (1-based) to evaluate.
- `--all_eff` : fixed effects column indices (same semantics as `main.sh` usage).
- `--ran_eff` : random effect column index (if applicable).
- `--rep` : number of CV repeats.
- `--fold` : number of CV folds.
- `--thread` : threads available to inner programs (DMU/mbBayesABLD).
- `--iter`, `--burnin` : MCMC iterations and burn-in for Bayesian methods.
- `--seed` : random seed for partitioning and sampling.
- `--debug` : dry-run / skip heavy GEBV computation (useful for testing I/O and folder layout).
- `--binf` : block/interval file (when using Bayesian multi-component models).
- `--tbv_col` : column in phenotype file containing TBV when not using `tbvf`.
- `--out` : output filename (overrides default `accur_<method>.txt`).

**Inputs expected**

- Breed-specific working folder (script will create if missing).
- `bfile` plink files (or `${label}m` created by `GP_cross` when extracting breed genotypes).
- Phenotype file matching individuals in the `.fam` (IDs must be consistent).
- Optional pedigree and `phe_adj_BLUP.txt` for `--tbvf`.

**Outputs**

- Per-fold/per-rep predictions and accuracy metrics (correlation, bias) saved under the breed folder.
- Logs capturing stdout/stderr from the prediction engine and wrapper.
- Optionally posterior samples or DMU outputs (depending on method).

**Example (direct execution)**

```
$GP_single --label YY --phef /path/to/YY_pheno.txt --bfile /path/to/YYm \
  --method GBLUP --seed 12345 --all_eff "2 1" --ran_eff "1" \
  --rep 10 --fold 5 --phereal 2 --thread 20 --tbvf /path/to/phe_adj_BLUP.txt \
  --iter 30000 --burnin 20000 --out accur_GBLUP.txt
```

**Notes and caveats**

- Plink `.fam` family IDs must include the breed label or match the IDs in the phenotype file; otherwise extraction/merging will fail.
- For Bayesian methods, ensure `--binf` or block partitioning is available and that `mbBayesABLD` binaries are executable.
- Use `--debug` to validate file paths and folder creation without consuming heavy CPU/memory.
- When running on HPC under SLURM, `--thread` is typically matched to `SLURM_CPUS_ON_NODE` and memory requests should be tuned according to method (Bayesian methods require more memory).

## 7.3 `GP_multi_breed.sh` (`$GP_multi`)

**Function:** Perform multi-breed genomic prediction workflows — either single-trait (`--type single`) or multi-trait (`--type multi`) analyses — integrating genotype and phenotype information across multiple populations. Supports GBLUP and Bayesian multi-component models (e.g., `mbBayesAB`) with flexible genome partitioning, priors, and covariance structure options.

**Core behaviors**

- Combine per-breed inputs (genotypes, phenotypes or `tbvf`) into multi-breed analyses, optionally using union/blend/merged genotype sets or breed-specific partitions.
- For single-trait multi-breed prediction: assemble training sets across specified `--pops`, run chosen method, and produce prediction accuracies for target populations.
- For multi-trait multi-breed prediction: estimate (co)variance structures across traits and breeds, support constrained residual structures (GBLUP, `--noCov`).
- Support different block/bin partitioning strategies (`--bin`: `fix`, `lava`, `cubic`) and use block files (`--binf`) when running `mbBayesAB` to allocate effects per block.
- Pass MCMC control (`--iter`, `--burnin`) and other model flags to the backend Bayesian engine (`mbBayesABLD`).

**Common parameters (non-exhaustive)**

- `--pops` : space-delimited list of populations/breeds forming the combined analysis (e.g., `"YY LL"`).
- `--type` : `single` or `multi` (analysis scope).
- `--method` : `GBLUP`, `mbBayesAB`, `bayesR`, etc.
- `--tbvf` : file with true/adjusted phenotypes used as TBV reference (per-trait).
- `--phereal` : index of target trait column to evaluate.
- `--nsnp_win`, `--win` : window parameters for sliding-window partitioning and LD-based summaries.
- `--thread` : threads for inner programs.
- `--iter`, `--burnin` : MCMC length and burn-in for Bayesian approaches.
- `--ref` : SNP panel reference used when partitioning (e.g., `M` or other).
- `--seed` : random seed.
- `--bin` : partition method (`fix`, `lava`, `cubic`).
- `--binf` : block file path (precomputed) to feed to Bayesian models.
- `--tbv_col` : column index in phenotype file that stores TBV (if using a global `phef` format).
- `--noCov` : constrain residual covariance between traits to zero (flag passed to GBLUP model).
- `--debug` : skip computationally heavy final runs (useful to test I/O and parameter passing).

**Inputs expected**

- Per-trait working folder created by `GP_cross` (e.g., `${proj}/${trait}`) containing per-breed subfolders produced by `within` runs.
- `tbvf` (adjusted phenotypes) files under the trait folder or a global phenotype with TBV columns.
- Block/bin files if Bayesian partitioned models are used.

**Outputs**

- Multi-breed prediction results and accuracy summaries saved within the trait folder (possibly under subfolders named by `--dirPre`, `--bin`, method).
- Variance component estimates and (co)variance matrices when `--type multi` is used (these are later aggregated by `varcomp_summary.sh`).
- Logs for the multi-breed run capturing backend engine output and wrapper messages.

**Example (direct execution)**

```
$GP_multi --pops "YY LL" --type multi --method mbBayesAB --tbvf /path/to/phe_adj_BLUP.txt \
  --phereal 2 --nsnp_win 100 --win 50 --thread 50 --iter 30000 --burnin 20000 \
  --ref M --seed 12345 --bin cubic --binf /path/to/cubic_blocks.txt --dirPre iter30000 --suffix
```

**Notes and caveats**

- Ensure that all per-breed `within` analyses have completed and their folders exist (GP_cross checks these before calling `GP_multi`). Missing per-breed outputs will cause the multi-breed step to be skipped.
- Bayesian multi-trait models can be memory- and time-intensive; assign sufficient `--thread` and cluster memory. When in doubt, start with a small test (reduced iter/burnin and fewer bins).
- The `--bin` strategy selection affects model computational cost and interpretation: `fix` (equal SNP count) and `lava`/`cubic` (smooth partitioning based on local LD) are supported; precompute block files with `GP_cross --type bin` or provide `--binf`.

**Practical integration**

- `main.sh` orchestrates calls to `GP_multi` for many breed-combinations and parameter sweeps (different `--bin`, `--method`, etc.). Use `--dirPre` to keep outputs from different parameterizations separate.
- After `GP_multi` runs, call `GP_cross --type var` and `GP_cross --type accur` to aggregate variance components and accuracy statistics across runs.

## 7.4 R scripts

Core scripts are described individually; remaining helpers are summarized in one paragraph. These scripts are normally invoked by the shell drivers (`GP_cross_validation.sh`, `GP_single_breed.sh`, `GP_multi_breed.sh`) via `Rscript`.

**adj_pheno_cal.R**
 Purpose: post-process DMU / BLUP outputs to compute adjusted phenotypes (EBV + residual) and produce `phe_adj_BLUP.txt` (or equivalent TBV files).
 Input: DMU/BLUP output files (`.SOL`, `.lst`, etc.) and raw phenotype file.
 Output: adjusted phenotype table used as `tbvf` for accuracy evaluation and downstream analyses.
 When used: by `GP_cross --type adj` and before accuracy/validation steps.

**pheno_simulation.R**
 Purpose: simulate phenotypes on provided genotypes (or QMSim output) with user-specified architecture: heritability, QTL count, inter-breed/trait genetic correlations, QTL overlap rules, etc.
 Input: genotype files or QMSim genotype list and simulation parameters (`--h2`, `--nqtl`, `--rg`, etc.).
 Output: `pheno_sim.txt`, `qtl_info.txt`, and optional `qtl_snpid.txt`.
 When used: by `GP_cross --type psim` to create ground-truth phenotypes for method evaluation.

**bayesR_file.R**
 Purpose: prepare BayesR-compatible input files (marker matrix, phenotype file, parameter/mapping files).
 Input: PLINK genotype and phenotype files.
 Output: BayesR input files and mapping tables.
 When used: prior to running BayesR.

**accuracy_bias_calculation.R**
 Purpose: compute prediction accuracy and bias statistics (Pearson r, regression slope/intercept, RMSE, etc.).
 Input: predicted GEBVs/EBVs and TBV/observed phenotype file.
 Output: per-fold/per-rep accuracy summary tables used by accuracy aggregation scripts.
 When used: after per-fold prediction outputs (called by wrappers or post-processing steps).

**Block / partitioning scripts (cubic_smoothing_block.R, fix_frq_ld_bolck.R, block_LD_cor.R)**
 Purpose: generate genome partition/block definition files using different strategies (cubic smoothing, fixed SNP-count windows), and compute block LD/correlation statistics used to choose or validate blocks.
 Input: SNP map, LD statistics or genotype data.
 Output: block definition files (used as `--binf`) and block summary statistics.
 When used: before partitioned Bayesian models (e.g., `mbBayesAB`) or when `GP_cross --type bin` is invoked.

**Sample filtering & matching scripts (keep_pheno_geno_individuals.R, geno_individuals_select.R, pheno_miss_remove.R)**
 Purpose: harmonize genotype and phenotype samples, remove missing phenotypes, and select individuals by breed/generation/criteria.
 Input: PLINK files and phenotype tables.
 Output: cleaned phenotype files, ID lists, or PLINK subsets ready for model fitting.
 When used: early in the pipeline (geno/psim/adj steps) to ensure consistent sample sets.

**Other R helpers (brief summary)**
 The repository also contains multiple smaller/auxiliary R scripts that perform focused tasks: plotting and diagnostics (`LD_decay_plot.R`, `PCA_plot.R`), local block genetic-correlation estimation (`local_block_rg.R`), column/trait correlation checks (`columns_correlation.R`), distance/statistics utilities (`mean_distance.R`), QMSim marker selection (`QMSim_mrk_select.R`), validation fold generation (`validation_population_define.R`), variance-prior file generation (`variance_prior_setting.R`), array/panel combination helpers (`array_combination.R`), and a few format-conversion utilities (`bayesR_file.R` overlaps here). These scripts are typically invoked by the shell drivers as needed and produce intermediate tables, configuration files, or plots consumed by downstream shell wrappers and model executables.

 [back to top](#user-manual)

# 8. Input file formats and examples

This section describes common input files and their formats. Example paths use `data/Real/Xie2021`.

## 8.1 Genotype (plink bfile)

- Prefix: `genotype` (expected `genotype.bed`, `.bim`, `.fam`)
- Format: standard plink binary.
- If your prefix differs, change `bfile` in `main.sh` or pass `--bfile`.

## 8.2 Phenotype file `phenotype.txt`

- Typically tab- or space-separated. First column is ID, subsequent columns are trait values or covariates.
- Example:

```
ID  MS    PFAI  100SdW  Yield
ind1  2.34  0.45  12.4   3.21
ind2  2.10  0.40  11.8   3.05
...
```

Scripts may split or filter phenotypes by breed based on `breeds.txt`.

## 8.3 `breeds.txt`

- A space-separated list of breed names, e.g., `YY LL VEC VEF ADP`.
- Script treats this as array `breeds` and iterates.

## 8.4 `traits_all.txt` and `traits_analysis.txt`

- `traits_all.txt`: list of all traits (for multi-trait models).
- `traits_analysis.txt`: subset of traits to analyze (for iterations, phenotype adjustment, etc.).
- Format: space-separated or newline-separated trait names.

## 8.5 `breeds_combination.txt`

- Lists breed combinations, one combination per line, breeds separated by spaces. Example:

```
YY LL
YY VEC
LL VEF
VEC VEF ADP
```

 [back to top](#user-manual)

# 9. Output files description

After a successful run, many files appear under `output/`, organized by source/trait/breed. Below are common outputs and their meanings.

## 9.1 `phe_adj_BLUP.*` (phenotype adjustment results)

- Files: `phe_adj_BLUP.SOL`, `phe_adj_BLUP.lst`, `phe_adj_BLUP.txt`, etc.
- `phe_adj_BLUP.lst` often contains variance component and heritability estimates with lines such as:

```
          Intra Class
   Trait  correlation     V(t)         SE(t)         SD-A          SD-P

     1     0.22265       0.00554       0.07444       0.22882       0.48495
```

- `V(t)`: estimated heritability (or variance component depending on script semantics)
- `SE(t)`: its standard error

## 9.2 `accur_*` (accuracy summaries)

- Files like `accur_within_*.txt`, `accur_single_*.txt`, `accur_multi_*.txt`.
- Contain correlation, regression slope, RMSE, bias, and other accuracy metrics for each method/combination.

## 9.3 `var_*` (variance component summaries)

- Aggregated variance components and heritability estimates across methods/combos for comparison.

## 9.4 `time_*` (runtime statistics)

- Records runtime per step/job for performance analysis.

## 9.5 Simulation outputs (Simulation section)

- For each replicate: `pheno_sim.txt`, `merge` bfile, block files (e.g., `cubic_M_50.txt`), etc.
- Organized under `output/repX/dist/corY` etc.

 [back to top](#user-manual)

# 10. Common run scenarios and example commands

Below are examples for running pipeline with Slurm and without Slurm.

## 10.1 With Slurm (recommended)

To run the entire pipeline (assuming `initialize.sh` already ran):

```
bash main.sh
```

`main.sh` will submit many `sbatch` tasks; monitor logs under `output/<source>/log/`.

## 10.2 Without Slurm (local testing)

- Edit `main.sh` to comment out `sbatch` lines and call `$GP_cross` directly.
- Reduce `--thread` to the number of available cores (e.g., `--thread 8`).
- Run `main.sh` in fragments to monitor progress.

## 10.3 Run a single step (examples)

Run phenotype adjustment:

```bash
$GP_cross --proj /path/to/output/Xie2021 --breeds "YY LL" --traits "MS PFAI" --trait MS --bfile /path/to/genotype --phef /path/to/phenotype.txt --all_eff "3 2 1" --type adj
```

Run multi-trait `mbBayesAB` with cubic blocks and burnin 20000:

```bash
$GP_cross --type multi --method mbBayesAB --proj /path/to/output/Xie2021 --breeds "YY LL" --traits "MS PFAI 100SdW Yield" --trait MS --code /path/to/code --tbv_col same --seed 12345 --thread 50 --bin cubic --burnin 20000 --iter 30000 --suffix
```

 [back to top](#user-manual)

# 11. Performance advice and parameter tuning

- `--thread`: number of threads used by a job. If a node has 50 cores, set to 50. If 20 cores, set to 20.
- `--rep` and `--fold`: affect total job count (`rep × fold × breeds × traits × methods`). Reduce if resources limited (e.g., `rep=5`).
- `--burnin` and `--iter`: Bayesian convergence parameters. Increasing yields more stable estimates at the cost of runtime.
- `--mem` (Slurm): Bayesian methods can be memory-intensive. Defaults use 80–100G; adjust as needed.
- Parallel strategy: If cluster supports many small jobs, submit multiple smaller jobs rather than one massive job for better resource utilization.

 [back to top](#user-manual)

# 12. Troubleshooting

## 12.1 Common errors and how to locate them

- **R package not found**: Install the missing package(s) in R or set the library path.
- **`sbatch: command not found`**: Slurm not installed on the machine. Comment out `sbatch` lines and run commands directly or use another scheduler.
- **Jobs killed (OOM)**: Check Slurm logs for memory overuse. Increase `--mem` or decrease parallelism.
- **`phe_adj_BLUP.SOL` not generated or empty**: Check logs for the `--type adj` run; ensure phenotype and genotype paths are correct.
- **GBLUP / mbBayesAB slow or interrupted**: Inspect per-job logs under `output/<source>/log/`.

## 12.2 Log file locations

- Global logs: `output/<source>/log/*.log` (each `sbatch` writes stdout/stderr here).
- Task-level outputs: subdirectories like `output/<source>/<trait>/<breed>/` contain method-specific logs and outputs.

 [back to top](#user-manual)

# 13. FAQ (frequently asked questions)

Q1. My genotype prefix is not `genotype`. How to change it?
 A1. Edit `bfile=${data_path}/genotype` in `main.sh` to your prefix, or pass `--bfile /path/to/yourprefix` to `$GP_cross`.

Q2. I don’t have Slurm. Can I still run everything?
 A2. Yes. Comment out or replace `sbatch` lines in `main.sh` and manage the execution (serial or custom parallel) yourself. Be cautious about resource exhaustion.

Q3. How to run only the simulation section?
 A3. Remove or comment out the real-data top loop in `main.sh` and run the Simulation block, or extract the Simulation portion into a separate script.

Q4. How to ensure reproducibility?
 A4. Use the same `random.seed` file (the script reads/writes this file) and save all runtime parameters. Logging commands and variables helps reproduce runs exactly.

Q5. Runs take too long. How to find bottlenecks?
 A5. Check `time_*` outputs for timing statistics; run single tasks wrapped by `time` for profiling; look at per-job logs for long steps.

 [back to top](#user-manual)

# 14. License and citation

This project is licensed under the **GPL-3.0 License**.
 See the `LICENSE` file in the repository or https://www.gnu.org/licenses/gpl-3.0.en.html for details.

**Citation:** Please cite the paper and optionally this repository when using the data or methods. Format depends on journal requirements.

 [back to top](#user-manual)

# 15. Contact and acknowledgements

**Primary contact (corresponding author):** Li Weining (Li Weining)

- Email: li.wn@qq.com
- ORCID: 0000-0002-0578-3812
- Affiliation: College of Animal Science and Technology, China Agricultural University, Beijing, China

If you encounter problems (unclear data, script errors, large result discrepancies), please email with:

- Your runtime environment (OS, R version, plink version, QMSim version)
- Script name and line number if available
- Relevant log file snippets (include context around error messages)

 [back to top](#user-manual)

------

# Appendix A: Example full run (step-by-step)

A minimal end-to-end example (assuming adequate resources and Slurm availability):

1. Clone and enter repository:

```bash
cd /public/home/liujf/liwn/code/GitHub
git clone git@github.com:CAU-TeamLiuJF/mbBayesABLD.git
cd mbBayesABLD
```

2. Initialize environment:

```bash
./initialize.sh
# Install any missing R packages as indicated
```

3. Breed data folder:

```bash
ls data/Real/Xie2021
# Ensure genotype.*, phenotype.txt, breeds.txt, traits_all.txt, traits_analysis.txt, breeds_combination.txt exist
```

4. Execute the main script (run in parts or whole):

```bash
bash main.sh
```

Please note that this will execute the full reproduction of all results presented in the paper and may take a considerable amount of time. We recommend running the script line by line or selectively executing specific tasks within the project.

5. Monitor logs:

```bash
tail -f output/Xie2021/log/phe_adj.log
tail -f output/Xie2021/log/within_MS_YY.log
```

6. After all jobs finish, run summaries:

```bash
$GP_cross --type accur --proj /public/home/liujf/liwn/code/GitHub/MBGP_CV/output/Xie2021 --bin "fix lava cubic" --breeds "YY LL VEC VEF ADP" --traits "MS PFAI 100SdW Yield"
$GP_cross --type var --proj /public/home/liujf/liwn/code/GitHub/MBGP_CV/output/Xie2021 --bin "fix lava cubic" --breeds "YY LL VEC VEF ADP" --traits "MS PFAI 100SdW Yield"
```

 [back to top](#user-manual)

# Appendix B: Practical tips

- Back up scripts before modifying: `cp main.sh main.sh.bak`.
- If only some sub-jobs fail, re-run those sub-jobs instead of the whole pipeline.
- Clean or use a fresh `output/` directory for major re-runs to avoid mixing results.
- Use `git` to version-control any changes to `GP_cross` or `main.sh` for reproducibility and tracking.
- For debugging, run `$GP_cross` on a small toy dataset (few individuals / few SNPs) interactively to verify behavior.

------

# Closing remarks

This README aims to cover environment setup, pipeline execution, inputs/outputs, post-processing, and troubleshooting so that users unfamiliar with the original paper can run and reproduce results.
 If you want further expansion of any part (for example, a full parameter reference for every `GP_cross --type` mode, exhaustive output field definitions for `var`/`accur`, additional R post-processing scripts, or automated deployment scripts), tell me which section to expand (e.g., “expand GP_cross parameter manual” or “add detailed var/accur output description”), and I will continue to refine the documentation.

 [back to top](#user-manual)
