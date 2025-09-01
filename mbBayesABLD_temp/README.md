<!--
 * @Author: Li Weining
 * @LastEditTime: 2024-09-26 13:36
 * @FilePath: \mbBayesABLD\Readme.md
-->

# User Manual

- [1. Software Introduction](#1-software-introduction)
- [2. Installation Instructions](#2-installation-instructions)
- [3. Quick Start](#3-quick-start)
  - [3.1 Installation](#31-installation)
  - [3.2 Running the Software](#32-running-the-software)
  - [3.3 Parameters of the Software](#33-parameters-of-the-software)
- [4. File Description](#4-file-description)
  - [4.2 Input Files](#42-input-files)
  - [4.2.1 Phenotype File](#421-phenotype-file)
  - [4.2.2 Genotype File](#422-genotype-file)
  - [4.2.3 Marker Chromosome and Physical Location File](#423-marker-chromosome-and-physical-location-file)
  - [4.2.4 Genomic Block Partitioning File](#424-genomic-block-partitioning-file)
  - [4.2.5 Variance Component Prior File](#425-variance-component-prior-file)
- [5. Example Code](#5-example-code)
  - [5.1 Running with Phenotype and Genotype Files](#51-running-with-phenotype-and-genotype-files)
  - [5.2 Running with a Specified Block Partition File](#52-running-with-a-specified-block-partition-file)

# mbBayesABLD: Multi-Breed Genomic Prediction Bayesian Model

User Manual for Multi-Breed Joint Prediction Models that Accounting for Heterogeneous Genetic (Co)variance Software V1.0

## 1. Software Introduction

**mbBayesABLD** is a multi-breed joint genomic prediction tool designed for livestock and poultry populations, but it is also applicable to plant populations. **mbBayesABLD** allows for more accurate sharing of information between breeds, thereby improving the accuracy of breeding value estimation for small-scale breeds within joint reference populations. In this model, chromosomes are first partitioned into blocks based on linkage disequilibrium (LD) information, and a genetic correlation estimate between populations is obtained for each independent block. Information is precisely shared between breeds through the correlations between the genetic marker effects.

This tool enables the calculation of genomic estimated breeding values (GEBVs) for all individuals based on their phenotypic and genotypic data. It provides the ability to effectively leverage information from multiple breeds, significantly enhancing breeding selection efficiency and the overall productivity of the industry.

 [back to top](#user-manual)

## 2. Installation Instructions

The program was developed on the `CentOS Linux release 7.9.2003` platform and is based on Linux version `3.10.0-1127.19.1.el7.x86_64` kernel version. It should be able to run on any Unix-like system.

Considering that the program is mainly developed in C language, it should also be runnable on Windows or macOS systems after recompilation.

This software is designed for genomic prediction and requires input in the form of phenotypic records and genotypic information of the analyzed individuals (typically in the form of SNP arrays). Additionally, due to matrix storage and computational processes, extra storage demands may arise during analysis. We recommend that the computer running the software has at least 10GB of available memory. For larger populations, more system memory will be required depending on the population size.

This software relies on the `Intel® oneAPI Math Kernel Library` for scientific calculations (no version requirement). Users need to ensure the existence of this dependencies in their runtime environment. 

[back to top](#user-manual)

## 3. Quick Start

Please note that, unless otherwise specified, all of the following operations are carried out in Linux command-line mode.

### 3.1 Installation

The software provides pre-compiled binary files that can be directly downloaded and used:

```bash
wget https://github.com/CAU-TeamLiuJF/mbBayesABLD/archive/refs/heads/main.zip
unzip main.zip
cd mbBayesABLD-main
chmod 755 ./bin/mbBayesABLD
```

### 3.2 Running the Software

To run the software, simply enter the relevant commands in the command line:

```bash
cd example
../bin/mbBayesABLD --bfile /path/to/plink/binary/genotype/file/prefix --phef /path/to/phenotype/file [options] ...
```
**bfile** is the genotype file parameter, specifying the genotype information file required by the software. **phef** is the phenotype file parameter, specifying the phenotypic record information file required by the software. Both of these parameters are mandatory for the software to run, and can be specified using either relative or absolute paths.

The genotype file input format should follow the PLINK binary file format, consisting of three files with the extensions .fam, .bim, and .bed, respectively. For more details on the file formats, please refer to the [PLINK File Format Documentation](https://www.cog-genomics.org/plink/1.9/formats).

The phenotype file input format is a plain text file, using spaces or tabs as delimiters. Whether the file is in DOS or UNIX format is not critical, although it is recommended to use the UNIX format for more stable software operation. The file should not include a header. For more detailed information regarding the input file formats, please refer to the "File Formats" chapter.

### 3.3 Parameters of the Software

Other relevant parameters for the software can be queried using the following command:

```bash
mbBayesABLD --help
```
The specific details of each parameter are as follows:

- #### Required Parameters:
  
    - `-l, --bfile`: PLINK binary file prefix (string)
      - Specifies the prefix for the PLINK binary file, for example `/path/plink`. The required genotype file should contain `.bed`, `.bim`, and `.fam` files.
    - `-p, --phef`: Phenotype file path (string)
      - Specifies the path to the phenotype file. It can be either a relative or an absolute path.
    
    #### Optional Parameters:
    
    - `-g, --raw`: PLINK raw file path (string)
      - Specifies the path for the PLINK raw file, e.g., `/path/plink.raw`.
    - `-T, --thread`: Number of threads (integer [=1])
      - Specifies the number of threads for the OpenMP run-time library to use during the computation.
    - `-y, --phe`: Phenotype position (character string)
      - Specifies the positions of phenotypes in the phenotype file, e.g., `'4 5'` for selecting the 4th and 5th columns as phenotypes.
    - `-F, --fix`: Fixed effect position (character string)
      - Specifies the positions of fixed effects in the phenotype file, e.g., `'2 3'`. By default, this is set to `NULL`.
    - `-R, --ran`: Environmental random effect position (character string)
      - Specifies the positions of environmental random effects in the phenotype file. Default is set to `NULL`.
    - `-a, --id`: ID position (integer [=1])
      - Specifies the position of the individual ID column in the phenotype file. Default is `1`.
    - `-v, --varf`: SNP variance prior file path (string)
      - Specifies the file containing SNP variance priors, default value is `varp * 0.5`.
    - `-e, --region`: SNP grouping method (character string [=LD])
      - Specifies the method for dividing SNPs into groups. Available options: `LD` or `fix`.
    - `-b, --binf`: SNP interval file path (string)
      - Specifies the file containing the number of SNPs per interval. Default is `NULL`.
    - `-j, --bincol`: SNP count column position in bin file (integer [=1])
      - Specifies the position of the SNP count column in the interval file. Default is `1`.
    - `-f, --mapf`: PLINK map file path (string)
      - Specifies the path to the `.map` or `.bim` file to define bins. Default is `NULL`.
    - `-N, --nsnp`: Number of SNPs per interval (integer [=100])
      - Specifies the number of SNPs per interval.
    - `-E, --effOut`: Posterior mean of location effects file (string [=theta.txt])
      - Specifies the file to output the posterior mean of location effects.
    - `-A, --varOut`: Posterior mean of dispersion matrix file (string [=var.txt])
      - Specifies the file to output the posterior mean of the dispersion matrix.
    - `-V, --gebvOut`: Breeding value output file (string [=gebv.txt])
      - Specifies the file for genomic estimated breeding value (GEBV) output.
    - `-c, --varC`: First `c` predictions in MCMC chain (integer)
      - Specifies the number of initial predictions to output from the MCMC chain.
    - `-M, --mcmcOut`: MCMC chain predictions output file (string [=MCMC_var.txt])
      - Specifies the file to output predictions in the MCMC chain.
    - `-m, --miss`: Missing value in phenotype data (double [=-99])
      - Specifies the value used to represent missing data in the phenotype file.
    - `-d, --dfvara`: Degrees of freedom for variance (Va) (double [=4])
      - Specifies the degrees of freedom for additive genetic variance.
    - `-D, --dfvare`: Degrees of freedom for environmental variance (Ve) (double [=4])
      - Specifies the degrees of freedom for residual (environmental) variance.
    - `-G, --diagG`: Large integer value (integer [=100000])
      - Specifies a sufficiently large integer for diagonalization steps.
    - `-S, --iter`: Length of MCMC chain (integer [=100])
      - Specifies the number of iterations for the MCMC chain.
    - `-J, --maf`: Minimum allele frequency (integer [=-0.01])
      - Excludes SNPs with a minor allele frequency lower than this value.
    - `-K, --geno`: SNP missing data proportion (integer [=0.99])
      - Specifies the allowed proportion of missing data for SNPs.
    - `-P, --mind`: Individual genotype missing data proportion (integer [=0.99])
      - Specifies the allowed proportion of missing genotypes for individuals.
    - `-B, --burnin`: Burn-in steps for MCMC (integer [=10])
      - Specifies the number of burn-in steps for the MCMC chain.
    - `-t, --thin`: Thinning rate for MCMC (integer [=10])
      - Specifies the thinning rate for the MCMC chain.
    - `-s, --seed`: Random number seed (long)
      - Specifies the seed for random number generation. Default is `NULL`.
    - `-o, --outdir`: Output directory (character string [=.])
      - Specifies the output directory for all generated files.
    - `-r, --report`: MCMC sampling report interval (integer [=100])
      - Specifies the interval for reporting MCMC sampling progress.
    - `-L, --logf`: Log file path (string)
      - Specifies the file to redirect standard output (log) to.
    - `-i, --prioSa`: Additive effect prior scale distribution (character string)
      - Specifies the prior distribution for the scale of additive effects. Options: `NULL`, `wishart`, or `gamma`.
    - `-I, --prioSe`: Residual effect prior scale distribution (character string)
      - Specifies the prior distribution for the scale of residual effects. Options: `NULL`, `wishart`.
    
    #### Flag Parameters:
    
    - `-n, --nocov`: Set residual covariances to zero
      - Constrains the covariance between residuals to zero.
    - `-C, --dontcenter`: Disable centering of the gene content matrix
      - Disables centering of the gene content matrix.
    - `-H, --header`: Specify if input files have headers
      - If present, automatically detects headers in the phenotype and bins files.
    

[back to top](#user-manual)

## 4. File Description

The relevant files include the required pedigree files and genotype files. Files should be ASCII-formatted files, separated by one or multiple spaces (or newline characters). The variable names for categories can include English letters and numbers. The individual numbers should be consistent across different files.

### 4.2 Input Files

The input files include required phenotype files, genotype files, genomic partitioning files, and prior variance component files. The phenotype file should be in ASCII format, with data fields separated by one or more spaces (or newline characters). Category variable names can include alphanumeric characters, and individual IDs should be consistent across all files.

#### 4.2.1 Phenotype File

The phenotype file should be a plain text file format, primarily used to provide individual IDs and the phenotypic values of traits to be analyzed. It can also include the identifiers for fixed effects and random effect groupings. The individual IDs should be placed in the first column, followed by columns containing phenotypic values. If fixed effects and other non-genetic random effects need to be fitted into the model, additional columns should be included in the file. Here is an example of a phenotype file:

```text
10004 1 1 1.96
10011 1 2 1.39
10019 1 2 -99
10023 1 1 2.22
10030 1 2 2.91
...
```

The columns represent the following:

- **ID**: Individual identifier
- **mu**: Intercept term
- **Sex**: Gender
- **y**: Phenotypic value

Missing phenotypic values are represented as `-99` (can be adjusted using the appropriate parameter).

#### 4.2.2 Genotype File

The genotype file is used for sampling marker effects during the calculation and for quality control of markers and individuals. It is formatted in one of two common PLINK file types:

**Binary format files** (`*.fam`, `*.bim`, `*.bed`). 

You only need to specify the file prefix, which can be a relative (`./geno/prefix`) or absolute path (`/storage1/data/geno/prefix`).

**012-encoded raw file (*.raw)**

The *.raw file format can be generated by PLINK using the following command:

```bash
plink --bfile 1 –recode A --out 1
```

The genotype file must contain breed markers in the column corresponding to family IDs, which indicate the breed affiliation of each individual.

#### 4.2.3 Marker Chromosome and Physical Location File

If you are using a `*.raw` genotype file in 012 format, and you want the software to divide SNPs into blocks based on the number of markers or linkage disequilibrium (LD) information, you need to provide an additional map file that specifies the chromosome and physical location for each marker. This map file follows the standard PLINK format, as shown below:

```text
1	M52  	0.03678	36780  
1	M129	0.10465	104650  
1	M294	0.23881	238810  
1	M330	0.26796	267960  
1	M374	0.3079	307900  
... 
```

For more details on the format, refer to the PLINK documentation:
https://www.cog-genomics.org/plink/1.9/formats.

#### 4.2.4 Genomic Block Partitioning File

The model includes heterogeneous genetic covariance, requiring chromosome segmentation into blocks. When users provide a genomic partitioning file, it should contain the number of markers in each block. This file contains a single column, where each line indicates the number of markers in a block. Here's an example:

```text
148
98
87
80
109
...
```

Alternatively, the file can have a more detailed format, like the following:

```text
1 1 320342 9318788 154
2 1 9318788 13417613 86
3 1 13417613 17786299 89
4 1 17786299 23858687 78
5 1 23858687 29949287 114
...
```

In this case, additional columns are allowed, as long as you specify which column contains the number of markers per block (e.g., using `--bincol 5`). The columns represent block ID, chromosome, start position, stop position, and the number of SNPs per block. The software uses only the last column, so you must specify which column contains the marker count. The software assigns markers to specified blocks sequentially from the genotype file.

If users provide both a custom block file and set genotype quality control parameters, the number of markers after quality control must match the numbers in the partitioning file. Otherwise, the program will exit with an error.

The software folder (`src`) contains a script for dividing the genome into blocks based on LD information. You can use the following command to generate a block partitioning file for use with the mbBayesABLD model:

```bash
cd ./mbBayesABLD-main/example  
../src/lava_LD_block.sh --bfile genotype --type LD --out ld_blocks.txt  
```

The first 6 lines of the output file will look like this:

```text
LOC CHR START STOP nSNP  
1 1 320342 9318788 154
2 1 9318788 13417613 86
3 1 13417613 17786299 89
4 1 17786299 23858687 78
5 1 23858687 29949287 114
...
```

#### 4.2.5 Variance Component Prior File

The software allows users to provide prior information for the variance components, including the prior scale matrix for the variance-covariance matrix of marker effects in each block. Each line represents a variance-covariance matrix. For example, the prior variance-covariance matrix for Block 1 with two breeds is:
$$
\sigma(b_1, b_1) \ \ \sigma(b_1, b_2) \ \ \sigma(b_1, b_2) \ \ \sigma(b_2, b_2)
$$
For three blocks, the complete prior file would look like this:

```text
0.10 0.10 0.10 0.20  
0.23 0.23 0.23 0.33  
0.26 0.26 0.26 0.36  
0.30 0.30 0.30 0.40  
```

The first three lines provide the prior variance for each block, and the last line represents the prior residual variance.

### 4.3 Output Files

The output files mainly include the breeding value files, fixed effect estimates, variance component files, parameter iteration estimates, log files, and more.

#### 4.3.1 Breeding Values

The breeding values for each breed are output into separate files. For joint prediction of three breeds, you can specify the following parameter for the output files:

```bash
--gebvOut EBV
```

The prefix "EBV" will be used for the breeding value output files, and the resulting files will be named `EBV_y1.txt`, `EBV_y2.txt`, and `EBV_y3.txt`, containing the breeding values for each breed. Here’s an example format for the file `EBV_y1.txt`:

```text
ID        EBV  
710387    2.406609E-01  
710392    4.330189E-01  
710396    3.861383E-01  
...
```

The file contains two columns with headers. The first column provides the individual IDs, and the second column contains the corresponding breeding values. All breeding value output files follow the same format.

#### 4.3.2 Posterior Means of Marker Effects

After the program completes successfully, it outputs the posterior means of the marker effects. The file format looks as follows:

```text
-1.509555E-05  -8.801209E-05  
-2.197057E-04   6.185192E-06  
2.724469E-04   -1.391984E-04  
...
```

The number of rows equals the number of quality-controlled markers, and the number of columns corresponds to the number of breeds.

#### 4.3.3 Fixed Effect Estimates (`sol.txt`)

If fixed effects (including intercept) were specified in the parameters, the posterior means of the fixed effects will be output to a designated file. Here is an example format of the fixed effect estimates file:

```
1.469106E+00   2.337024E+00  
-3.382830E+00  -3.398553E+00  
3.350822E+00   2.968550E+00  
3.702356E+00   -2.376497E+00  
-2.281983E+00  2.093072E-01  
```

This file contains three fixed effects: population mean, gender, and breed. The levels of each fixed effect are 1, 2, and 2, respectively, resulting in a total of 5 rows. Each row represents the estimated value for one level of a fixed effect across two “traits” (which in practice are different breeds).

#### 4.3.4 Variance Component File

The format of the output posterior means of variance components follows the same structure as the input prior variance component file. For more details, refer to section 5.2.4.

#### 4.3.5 Parameter Iteration Estimates

If the user specifies the output of all iterations of selected parameters using the `--mcmcOut` option, the values for each iteration are written to the specified file. Each row represents the estimated values for unknown parameters during a single iteration, and the number of columns corresponds to the number of unknown parameters specified. The number of rows equals `(iter-burnin)/thin`.

#### 4.3.6 Log Files

If the user does not want the log information to be printed directly to the console during runtime, they can redirect the log information to a file by specifying a path for the log file. The format of the log file is as follows:

```text
Program mbBayesABLD
Run started at Sun Dec  3 14:33:51 2023
Loading genotypes file...
Filtering SNPs or samples...

0 variants removed due to minor allele threshold (--maf).
0 variants removed due to missing genotype data (--geno).
0 samples removed due to missing genotype data (--mind).
40679 variants and 1800 people pass filters and QC.

Input genotype file  : merge
Input phenotype file : pheno.txt
output directory     : .
Bins file            : bins.txt
No. of individuals   : 1440
No. of loci          : 40679
No. of regions       : 185
Missing phenotype    : -99.000000
Prior Vara           : 0.0001    0.0000    0.0000    0.0000    0.0001    0.0000    0.0000    0.0000    0.0001    
Prior Vare           : 1.6202    0.0000    0.0000    0.0000    1.5124    0.0000    0.0000    0.0000    1.5863    
No. of cycles        : 100
Burnin               : 10
Thinning rate        : 10
Seed                 : 8123

### Begin Gibbs sampling...
Iter 20: 3 s    elasp: 0 min    eta: 0 min    First 2 fixed and 2 random: 1.783 1.138 0.009 0.012
Iter 40: 2 s    elasp: 0 min    eta: 0 min    First 2 fixed and 2 random: 1.664 1.059 0.007 0.010
Iter 60: 2 s    elasp: 0 min    eta: 0 min    First 2 fixed and 2 random: 1.680 1.090 0.005 0.004
Iter 80: 2 s    elasp: 0 min    eta: 0 min    First 2 fixed and 2 random: 1.692 1.063 0.004 -0.006
Iter 100: 2 s   elasp: 0 min    eta: 0 min    First 2 fixed and 2 random: 1.586 1.090 -0.003 0.002

MCMC finished.

Run ended at Sun Dec  3 14:34:09 2023
```

[back to top](#user-manual)

## 5. Example Code

Sample datasets can be found in the `example` folder.

### 5.1 Running with Phenotype and Genotype Files

```bash
mbBayesABLD --bfile genotype --phef phenotype.txt
```

```text
Program mbBayesABLD
Run started at Sun Sep 29 11:49:23 2024
Loading genotypes file...
Filtering SNPs or samples...

0 variants removed due to minor allele threshold (--maf).
0 variants removed due to missing genotype data (--geno).
0 samples removed due to missing genotype data (--mind).
37304 variants and 200 people pass filters and QC.

Input genotype file  : genotype
Input phenotype file : phenotype.txt
output directory     : .
No. of individuals   : 196
No. of loci          : 37304
No. of regions       : 372
Missing phenotype    : -99.000000
Prior Vara           : 0.0001	0.0000	0.0000	0.0001	
Prior Vare           : 0.8265	0.0000	0.0000	0.8501	
No. of cycles        : 100
Burnin               : 10
Thinning rate        : 10
Seed                 : 1727581763

### Begin Gibbs sampling...
Iter 100: 4 s	elasp: 0 min	eta: 0 min	First 2 fixed and 2 random: 2.144 1.711 -0.003 0.006

MCMC finished.

Run ended at Sun Sep 29 11:49:28 2024
```

The output files are the same as in the previous example.

### 5.3 Running with a Specified Block Partition File

```bash
mbBayesABLD --bfile genotype --phef phenotype.txt --binf ld_blocks.txt --bincol 5
```

```text
Program mbBayesABLD
Run started at Sun Sep 29 11:54:35 2024
Loading genotypes file...
Filtering SNPs or samples...

0 variants removed due to minor allele threshold (--maf).
0 variants removed due to missing genotype data (--geno).
0 samples removed due to missing genotype data (--mind).
37304 variants and 200 people pass filters and QC.

Input genotype file  : genotype
Input phenotype file : phenotype.txt
output directory     : .
No. of individuals   : 196
No. of loci          : 37304
No. of regions       : 442
Missing phenotype    : -99.000000
Prior Vara           : 0.0000	0.0000	0.0000	0.0000	
Prior Vare           : 0.8265	0.0000	0.0000	0.8501	
No. of cycles        : 100
Burnin               : 10
Thinning rate        : 10
Seed                 : 1727582075

### Begin Gibbs sampling...
Iter 100: 4 s	elasp: 0 min	eta: 0 min	First 2 fixed and 2 random: 2.288 1.841 -0.007 0.003

MCMC finished.

Run ended at Sun Sep 29 11:54:40 2024
```

[back to top](#user-manual)

**Note:**
This software is made available solely for the purpose of reproducing the results presented in the scientific article titled _Enhancing Multi-Breed Genomic Prediction for Small-Scale Breeds by Modeling Heterogeneous Genetic (Co)Variance Blockwise Accounting for Linkage Disequilibrium_. Any other use of this software, including commercial or non-research applications, is strictly prohibited without prior permission. For inquiries or permission requests, please contact the author at: liujf@cau.edu.cn.
