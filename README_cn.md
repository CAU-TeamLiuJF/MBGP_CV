# mbBayesAB：拟合异质遗传（协）方差的多品种联合评估模型

用于复现文章《Multi-trait Bayesian models enhance the accuracy of genomic prediction in multi-breed reference populations》的数据和脚本、程序文件等

引用：Li, W.; Zhang, M.; Du, H.; Wu, J.; Zhou, L.; Liu, J. Multi-Trait Bayesian Models Enhance the Accuracy of Genomic Prediction in Multi-Breed Reference Populations. Agriculture **2024**, 14, 626. https://doi.org/10.3390/agriculture14040626

## ⚠️ 声明

> 本项目相关技术已申请专利。**除个人阅读外**，任何使用、复制、修改或商业应用前，**必须**联系通讯作者获取授权。  
> 📧 联系邮箱：liujf@cau.edu.cn

# 运行环境

所有的脚本和程序文件都需要运行在Linux系统中，测试环境为：

**CentOS Linux release 7.6.1810 (Core)**

同时配置了[Slurm](https://slurm.schedmd.com/documentation.html)作业管理系统

# 运行步骤

## 1.下载/克隆GitHub仓库

### 1.1 下载压缩包 

通过链接[https://github.com/CAU-TeamLiuJF/mbBayesABLD/archive/refs/heads/main.zip](https://github.com/CAU-TeamLiuJF/mbBayesABLD/archive/refs/heads/main.zip)下载包含数据、脚本和程序文件的压缩包，然后解压到指定路径中。假设将本地电脑下载的压缩包复制到远程服务器/public/home/liujf/liwn/download路径中，然后可用以下命令将文件夹解压至/public/home/liujf/liwn/code/GitHub路径：

```bash
cd /public/home/liujf/liwn/download
unzip -d /public/home/liujf/liwn/code/GitHub mbBayesABLD-main.zip
mv /public/home/liujf/liwn/code/GitHub/mbBayesABLD-main /public/home/liujf/liwn/code/GitHub/mbBayesABLD ## 修改文件夹名称
```

### 1.2 克隆

假设Linux系统已连接互联网，可以通过以下命令克隆Github仓库：

```bash
cd /public/home/liujf/liwn/code/GitHub
git clone git@github.com:CAU-TeamLiuJF/mbBayesABLD.git
```

## 2.脚本初始化

### 2.1 运行初始化脚本initialize.sh

为了修改脚本中的工作目录，以及直接在命令行终端调用R语言脚本，需要根据当前计算机环境更新bash脚本中的部分路径及R语言脚本解释器的路径。注意，此时要求命令终端可以直接通过R调用R语言，即R语言安装路径已加载到环境变量中：

```bash
R
# R version 4.3.1 (2023-06-16) -- "Beagle Scouts"
# Copyright (C) 2023 The R Foundation for Statistical Computing
# Platform: x86_64-pc-linux-gnu (64-bit)
# 
# R是自由软件，不带任何担保。
# 在某些条件下你可以将其自由散布。
# 用'license()'或'licence()'来看散布的详细条件。
# 
# R是个合作计划，有许多人为之做出了贡献.
# 用'contributors()'来看合作者的详细情况
# 用'citation()'会告诉你如何在出版物中正确地引用R或R程序包。
# 
# 用'demo()'来看一些示范程序，用'help()'来阅读在线帮助文件，或
# 用'help.start()'通过HTML浏览器来看帮助文件。
# 用'q()'退出R.
# 
# [原来保存的工作空间已还原]
# 
# > 
```

同时该步骤还会检查当前R语言版本中是否已安装所需的R包。

首先将工作路径切换到项目文件夹的主目录，即脚本initialize.sh所在的目录，如：

```bash
cd /public/home/liujf/liwn/code/GitHub/mbBayesABLD  ## 需根据实际修改
./initialize.sh
```

### 2.2 根据运行环境修改/注释某些命令

项目的主脚本为main.sh，本项目的开发环境为CentOS 7，同时配置了Slurm作业管理系统，计算节点的核心数在50以上。若是想要在没有配置Slurm作业管理系统的Linux系统中运行命令，则需要注释脚本main.sh中sbatch命令所在的行，如72行：

```bash
...
    ## Parameters need to be ...
    ## Note: If the Slurm workload manager is not installed ...
#    sbatch -c2 --mem=4G \
    $GP_cross \
      --proj ${pro} \
...
```

预测准确性由10次重复的五折交叉验证步骤计算得到，因此一种情形需要运行10x5个子进程，若是运行环境计算机中CPU数目小于50，请修改--thread参数，如核心数为20，则修改为--thread 20

## 3.运行程序

在Linux命令行终端中逐行运行主脚本main.sh中的命令即可。

# 联系方式

李伟宁 li.wn@qq.com
