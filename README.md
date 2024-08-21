# mbBayesABLD: Multi-breed joint prediction model fitting heterogeneous genetic (co)variance
Data and program needed to reproduce the research results of the article "xxx"

cite: 

# Running environment

All scripts and program files need to run on the Linux system, and the testing environment is:

**CentOS Linux release 7.6.1810 (Core)**

configured [Slurm](https://slurm.schedmd.com/documentation.html) Job Management System

# run steps

## 1.Download/Clone GitHub repository

### 1.1 Download compressed package

Download the compressed package containing data, scripts and program files from the link [https://github.com/CAU-TeamLiuJF/mbBayesABLD/archive/refs/heads/main.zip](https://github.com/CAU-TeamLiuJF/mbBayesABLD/archive/refs/heads/main.zip), and then decompress it to the specified path. Assuming that the compressed package downloaded from the local computer is copied to the remote server /public/home/liujf/liwn/download path, you can then use the following command to decompress the folder to the /public/home/liujf/liwn/code/GitHub path:

```bash
cd /public/home/liujf/liwn/download
unzip -d /public/home/liujf/liwn/code/GitHub mbBayesABLD-main.zip
mv /public/home/liujf/liwn/code/GitHub/mbBayesABLD-main /public/home/liujf/liwn/code/GitHub/mbBayesABLD ## Change folder name
```
### 1.2 Clone 

Assuming that the Linux system is connected to the Internet, you can clone the Github repository with the following command:

```bash
cd /public/home/liujf/liwn/code/GitHub
git clone git@github.com:CAU-TeamLiuJF/mbBayesABLD.git
```

## 2.Script initialization

### 2.1 Run the initialization script initialize.sh

In order to modify the working directory in the script and call the R language script directly in the command line terminal, it is necessary to update some paths in the bash script and the path of the R language script interpreter according to the current computer environment. Note that at this time, the command terminal is required to be able to call R language directly through R, that is, the R language installation path has been loaded into the environment variable:

```bash
R
# R version 4.3.1 (2023-06-16) -- "Beagle Scouts"
# Copyright (C) 2023 The R Foundation for Statistical Computing
# Platform: x86_64-pc-linux-gnu (64-bit)
# 
# ...
# 
# > 
```

At the same time, this step will also check whether the required R packages have been installed in the current R language version.

First, switch the working path to the main directory of the project folder, that is, the directory where the script initialize.sh is located, such as:

```bash
cd /public/home/liujf/liwn/code/GitHub/mbBayesABLD  ## Need to be modified
./initialize.sh
```

### 2.2 Modify/comment some commands according to the running environment

The main script of the project is main.sh. The development environment of this project is CentOS 7, and the Slurm workload manager is installed. The number of cores of the computing nodes is more than 50. If you want to run commands in a Linux system without the Slurm job management system, you need to comment out the line where the sbatch command is located in the script main.sh, such as line 72:

```bash
...
    ## Parameters need to be ...
    ## Note: If the Slurm workload manager is not installed ...
#    sbatch -c2 --mem=4G \
    $GP_cross \
      --proj ${pro} \
...
```

The prediction accuracy is calculated by 10 repetitions of the 5-fold cross validation step, so one case needs to run 10x5 subprocesses. If the number of CPUs in the operating environment computer is less than 50, please modify the --thread parameter. For example, if the number of cores is 20, change it to --thread 20

## 3.Run the program

Run the commands in the main script **main.sh** line by line in the Linux command line terminal.

# Contact

Liweining li.wn@qq.com
