#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript
## phenotype simulation  ##
## liwn 2021-09-22 ##

## 获取命令行参数
spec <- matrix(
  c(
    "gt",       "A", 1, "character", "[Required] qtl plink file prefix of A",
    "h2",       "1", 1, "character", "[Optional] heritability of populations [0.5]",
    "rg",       "r", 1, "character", "[Optional] genetic correlation between breeds [0.2]",
    "win",      "w", 1, "integer",   "[Optional] number of snp in each bins [100]",
    "min",      "M", 1, "integer",   "[Optional] minimum SNPs in the bins where rg exist [nsnp_cor]",
    "bin",      "b", 1, "character", "[Optional] win / chr [win]",
    "binf",     "B", 1, "character", "[Optional] ld block file path [NULL]",
    "binc",     "C", 1, "integer",   "[Optional] column number of nSNP in block file [last col]",
    "miss",     "m", 1, "character", "[Optional] miss value in phenotype [-99]",
    "dist_cor", "c", 1, "character", "[Optional] distribution of genetic correlation of qtl [identical/uniform]",
    "nbin_cor", "q", 1, "integer",   "[Optional] [0.05]",
    "nsnp_cor", "Q", 1, "integer",   "[Optional] [10]",
    "nqtl",     "p", 1, "integer",   "[Optional] Proportion of QTL to total SNPs [0.01]",
    "seed",     "s", 1, "integer",   "[Optional] seed of random effect sampling",
    "mean",     "e", 1, "character", "[Optional] population mean of poopulations [1.0]",
    "out",      "o", 1, "character", "[Optional] output phenotype file name",
    "qtlf",     "t", 1, "character", "[Optional] output QTL effects file name",
    "overlap",  "v", 0, "logical",   "[Optional] allow QTLs and markers to overlap",
    "evenly",   "E", 0, "logical",   "[Optional] Distribute QTLs evenly within the selected bins",
    "fid",      "f", 0, "logical",   "[Optional] whether output family id",
    "help",     "h", 0, "logical",   "This is Help!"
  ),
  byrow = TRUE, ncol = 5
)
opt <- getopt::getopt(spec = spec)

## 检查参数
n_need <- length(c(opt$gt))
if (!is.null(opt$help) || n_need < 1) {
  cat(paste(getopt::getopt(spec = spec, usage = TRUERUE), "\n"))
  quit()
}

## 自定义函数
argp_parse <- function(char, num, sep = " ", label = "parameters") {
  argv <- as.numeric(unlist(strsplit(char, sep)))

  if (length(argv) == 1) {
    return(rep(argv, num))
  } else if (length(argv) != num) {
    cat("The number of", label, "should be:", num, "!\n")
    quit()
  } else {
    return(argv)
  }
}

find_local_extreme <- function(x, diff = 0.2, type = "min") { # nolint
  ## 寻找极值点
  n <- length(x)
  peaks <- which(diff(sign(diff(x))) < 0) + 1
  valleys <- which(diff(sign(diff(x))) > 0) + 1

  ## 加上端点
  if (x[1] > x[2]) {
    peaks <- c(1, peaks)
  } else if (x[1] > x[2]) {
    valleys <- c(1, valleys)
  }
  if (x[n] < x[n - 1]) {
    valleys <- c(valleys, n)
  } else if (x[n] > x[n - 1]) {
    peaks <- c(peaks, n)
  }

  final <- c()
  for (i in valleys) {
    left_peak <- TRUE
    if (any(peaks < i)) {
      left_peak <- max(peaks[peaks < i])
      diff_rate <- abs(x[i] - x[left_peak]) / x[i]
      if (diff_rate < diff) left_peak <- FALSE
    }

    right_peak <- TRUE
    if (any(peaks > i)) {
      right_peak <- min(peaks[peaks > i])
      diff_rate <- abs(x[i] - x[right_peak]) / x[i]
      if (diff_rate < diff) right_peak <- FALSE
    }

    if (type == "min" && left_peak && right_peak) {
      final <- c(final, i)
    } else if (type == "max" && (left_peak || right_peak)) {
      final <- c(final, i)
    }
  }

  ## 端点处不要
  final <- final[final != 1 & final != n]

  return(final)
}

## 加载需要的程序包（未安装的话请提前安装）
cat("Loading required packages... \n")
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("MASS"))
suppressPackageStartupMessages(library("Matrix"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyr"))

## 默认参数
if (is.null(opt$seed)) opt$seed <- Sys.time()
if (is.null(opt$miss)) opt$miss <- -99
if (is.null(opt$h2)) opt$h2 <- "0.5"
if (is.null(opt$rg)) opt$rg <- "0.2"
if (is.null(opt$nbin_cor)) opt$nbin_cor <- 0.05 ## 大于1为具体存在相关的区间数，小于1为相关的区间占总区间的比例
if (is.null(opt$nqtl)) opt$nqtl <- 0.005 ## QTL数占总SNP数的比例
if (is.null(opt$nsnp_cor)) opt$nsnp_cor <- 10
if (is.null(opt$min)) opt$min <- opt$nsnp_cor
if (is.null(opt$dist_cor)) opt$dist_cor <- "identical"
if (is.null(opt$mean)) opt$mean <- "1.0"
if (is.null(opt$bin)) opt$bin <- "win"
if (is.null(opt$win)) opt$win <- 100
if (!is.null(opt$binf)) opt$bin <- opt$binf

## 随机数种子
set.seed(as.integer(opt$seed))
cat("random number seed: ", opt$seed, "\n")

## 检查基因型文件是否存在
bfiles <- trimws(opt$gt, which = c("both"), whitespace = "[ \t\r\n]")  ## 去除首尾空格
bfiles <- unlist(strsplit(bfiles, " "))
for (prefix in bfiles) {
  for (suffix in c(".map", ".ped")) {
    f <- paste0(prefix, suffix)
    if (!file.exists(f)) {
      cat("file needed not found:", f, "\n")
      quit()
    }
  }
}

## 品种数
np <- length(bfiles)

## 相关系数矩阵上三角元素个数
nrg <- choose(np, 2)

## 参数解析
rgs <- argp_parse(opt$rg, nrg, label = "rg")
h2s <- argp_parse(opt$h2, np, label = "h2")
means <- argp_parse(opt$mean, np, label = "mean")

## 报告遗传力等参数
cat("\ncorrelations of QTL effect: ", rgs, "\n")
cat("heritabilities: ", h2s, "\n")
cat("population means: ", means, "\n")
cat("bins definition: ", opt$bin, "\n")
cat("bins with genetic correlation: ", opt$nbin_cor, "\n")
cat("number of SNPs with genetic correlation in each bins: ", opt$nsnp_cor, "\n")
if (opt$bin == "win") cat("nsnp in bins: ", opt$win, "\n")
cat("distribution of genetic correlation : ", opt$dist_cor, "\n\n")

## 读取plink格式基因型
cat("Loading plink genotype files... \n")
maps <- list()
peds <- list()
nsnps <- c()
for (i in 1:np) {
  maps[[i]] <- fread(paste0(bfiles[i], ".map"))
  peds[[i]] <- fread(paste0(bfiles[i], ".ped"))
  names(maps[[i]]) <- c("CHR", "SNP", "cM", "POS")
  nsnps <- c(nsnps, nrow(maps[[i]]))
}

## 共有标记位点
if (length(unique(nsnps)) > 1) {
  cat("The number of markers in the genotype file is inconsistent with that in the ldblock file!\n")
  quit()
}

## 生成snp所属区间
map <- maps[[1]]
chrs <- unique(map$CHR)
if (opt$bin == "win") {
  ## 根据提供win进行区间划分，即每个区间的SNP数为win
  bin <- nsnp_bin <- NULL
  for (i in seq_along(chrs)) {
    nsnp_i <- sum(map$CHR == chrs[i])  ## 染色体i上的SNP数
    win_num <- floor(nsnp_i / opt$win) ## 染色体i可划分的区间数
    if (win_num < 1) win_num <- 1
    nsnp_bini <- rep(opt$win, win_num)
    nsnp_bini[win_num] <- nsnp_bini[win_num] + nsnp_i - sum(nsnp_bini)
    nsnp_bin <- c(nsnp_bin, nsnp_bini)
  }
  bin <- data.frame(bin = seq_len(length(nsnp_bin)), nsnp = nsnp_bin)
} else if (opt$bin == "chr") {
  ## 每条染色体为一个区间，每条染色体上的QTL数相等
  bin <- as.data.frame(table(as.factor(map$CHR)))
  opt$nbin_cor <- nrow(bin)
  names(bin) <- c("bin", "nsnp")
} else if (file.exists(opt$bin)) {
  ## 根据提供的区间划分文件，在每个区间中选取QTL
  bin <- fread(paste0(opt$bin))
  plink_ld_name <- c("CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2")
  if (all(names(bin) %in% plink_ld_name)) {
    ## 根据局部LD情况确定SNP所属区间
    cat("calculating mean LD in windows...\n")
    map$index <- FALSE
    for (chr in chrs) {
      snps <- map$SNP[map$CHR == chr]
      ldi <- bin[bin$CHR_A == chr]
      mapi <- map[map$CHR == chr, ]
      for (j in seq_along(snps)) {
        index <- (max(j - opt$win, 1)):(min(j + opt$win, length(snps)))
        mapi$R2[j] <- mean(unlist(subset(ldi, SNP_A %in% snps[index] & SNP_B %in% snps[index], select = "R2")))
      }
      max <- find_local_extreme(mapi$R2, type = "max")
      map$index[map$CHR == chr & map$SNP %in% mapi$SNP[max]] <- TRUE
    }

    ## 计算每个SNP与上下游SNP的LD均值
    # ld_mean <- bin %>%
    #   group_by(SNP_A) %>%
    #   summarize(CHR = unique(CHR_A), R2 = mean(R2, na.rm = TRUE))
    # ld_mean2 <- bin %>%
    #   group_by(SNP_B) %>%
    #   summarize(R2 = mean(R2, na.rm = TRUE))
    # ld_mean_2 <- left_join(ld_mean, ld_mean2, by = c("SNP_A" = "SNP_B"))
    # ld_mean_m <- ld_mean_2 %>%
    #   rowwise() %>%
    #   mutate(mean_R2 = mean(c(R2.x, R2.y), na.rm = TRUE)) %>%
    #   select(CHR, SNP_A, mean_R2)
    # names(ld_mean_m) <- c("CHR", "SNP", "R2")
    # map <- left_join(map, ld_mean_m, by = c("CHR", "SNP"))
    # ## 填充缺失值
    # map <- map %>%
    #   fill(R2, .direction = "down")
    # map$index <- FALSE
    # for (i in unique(map$CHR)) {
    #   mapi <- map[map$CHR == i, ]
    #   mapi$R2s <- smooth.spline(mapi$R2, spar = 0.20)$y
    #   max <- find_local_extreme(mapi$R2s, type = "max")
    #   map$index[map$CHR == i & map$SNP %in% mapi$SNP[max]] <- TRUE
    #   # max <- find_local_extreme(mapi$R2s, type = "min")
    #   # map$index[map$CHR == i & map$SNP %in% mapi$SNP[max]] <- 2
    # }
    index <- unlist(sapply(which(map$index), function(x) (x - opt$win + 1):(x + opt$win + 1)))
    index <- index[index > 0 & index <= nrow(map)]
    map$index <- FALSE
    map$index[index] <- TRUE

    map$index <- c(FALSE, map$index[-length(map)])
    map$bin <- cumsum(!map$index)
  } else {
    if (is.null(opt$binc)) opt$binc <- ncol(bin)
    names(bin)[opt$binc] <- "nsnp"
    bin_nsnp_above <- which(bin$nsnp > opt$min)
  }
} else {
  cat("wrong parameter for bin:", opt$bin, "\n")
  quit()
}

## 指明每个SNP所属的区间
if ("bin" %in% names(map)) {
  cat("number of regions: ", max(map$bin), "\n")
} else {
  cat("number of regions: ", nrow(bin), "\n")
  map$bin <- rep(seq_len(nrow(bin)), times = bin$nsnp)
}

### 挑选QTL ###
## 每个区间挑选一个标记位点
qmap <- map %>%
  subset(!is.na(bin)) %>%
  group_by(bin) %>%
  slice_sample()

## 随机挑选nbin_cor个区间
bins_sel <- c()
if (opt$nbin_cor > 0) {
  if (opt$nbin_cor <= 1) {
    ## 提供的参数为百分比，根据百分比计算要挑选的区间数
    opt$nbin_cor <- floor(opt$nbin_cor * nrow(bin))
    cat("The number of bins with non-zero correaltion is:", opt$nbin_cor, "\n")
  }

  ## snp数符合要求的区间
  if (!exists("bin_nsnp_above")) {
    bin_nsnp_above <- map$bin[duplicated(map$bin)]
  }

  ## 可供挑选存在遗传相关的QTL的区间
  bin_snps <- unique(map$bin[map$bin %in% bin_nsnp_above])

  ## 检测区间数是否够抽样所需
  if (length(bin_snps) < opt$nbin_cor) {
    cat("The number of candidate intervals cannot meet the sampling requirements\n")
    cat(length(bin_snps), "/", opt$nbin_cor, "\n")
    quit()
  }

  ## 检测区间数是否够抽样所需
  if (length(bin_snps) < opt$nbin_cor) {
    cat("The number of candidate intervals cannot meet the sampling requirements\n")
    cat(length(bin_snps), "/", opt$nbin_cor, "\n")
    quit()
  }

  ## 抽取nbin_cor个区间，这些区间内的QTL在品种间存在关联
  bins_sel <- sample(qmap$bin[qmap$bin %in% bin_snps], opt$nbin_cor)

  ## 为保证区间内相关，再挑选nsnp_cor-1个标记位点
  if (opt$nsnp_cor > 1) {
    if (is.null(opt$evenly)) {
      qmap2 <- subset(map, bin %in% bins_sel & !SNP %in% qmap$SNP) %>%
        subset(!is.na(bin)) %>%
        group_by(bin) %>%
        slice_sample(n = opt$nsnp_cor - 1)
      qmap <- bind_rows(qmap, qmap2)
    } else {
      for (i in seq_along(bins_sel)) {
        bini <- bins_sel[i]
        snp_count <- length(map$bin[map$bin == bini])

        # 计算均匀分布的间隔
        interval <- floor(snp_count / (opt$nsnp_cor))

        # 按间隔选择QTL
        qtl_indices <- seq(from = 1, to = snp_count, by = interval)
        qtl_indices <- qtl_indices[1:(opt$nsnp_cor)]
        if (length(qtl_indices) != opt$nsnp_cor) break

        # 根据选定的QTL索引从区间内选择QTL
        qmap2 <- map %>%
          filter(bin == bini) %>%
          slice(qtl_indices)

        ## 删除原来挑选的单个SNP
        qmap <- subset(qmap, bin != bini)

        # 将选定的QTL添加到结果数据框
        qmap <- bind_rows(qmap, qmap2)
      }
    }
  }
}

## qtl数
if (opt$nqtl < 1) {
  nqtl <- round(opt$nqtl * nrow(map))
} else {
  nqtl <- opt$nqtl
}

## 保证qtl数符合要求
over <- nrow(qmap) - nqtl
if (over > 0) {
  keep_single <- nqtl - sum(qmap$bin %in% bins_sel)
  keep_qtl <- sample(qmap$SNP[!qmap$bin %in% bins_sel], keep_single)
  qmap <- subset(qmap, SNP %in% keep_qtl | bin %in% bins_sel)
} else if (over < -0.5 * nrow(bin)) {
  over <- abs(over)
  add_each_bin <- ceiling(over / nrow(bin))
  qmap3 <- subset(map, !SNP %in% qmap$SNP) %>%
    subset(!is.na(bin)) %>%
    group_by(bin) %>%
    slice_sample(n = add_each_bin)
  qmap <- bind_rows(qmap, qmap3)
}

## 排序
qmap <- qmap[order(qmap$bin, qmap$POS), ]

## 报告QTL数目
cat("total number of qtl: ", nrow(qmap), "\n")

## 提取QTL基因型信息
qtl <- ebvs <- list()
fids <- c()
for (i in 1:np) {
  qtl_index <- which(maps[[i]]$SNP %in% qmap$SNP)
  col_index <- c(1:6, 2 * qtl_index + 6, 2 * qtl_index - 1 + 6)

  ## QTL基因含量矩阵
  qtl[[i]] <- subset(peds[[i]], select = col_index)
  ebvs[[i]] <- qtl[[i]][, 2:1]
  names(ebvs[[i]]) <- c("iid", "fid")

  ## 品种ID
  fids <- c(fids, unique(qtl[[i]]$V1))
}

## 判断fid合法性
if (length(fids) != np) {
  cat("The 'fid' in the ped file is not a population identifier!\n")
  quit()
}

## 加性效应协方差阵
cormat <- diag(np)
index <- 1
for (i in 1:(np - 1)) {
  for (j in (i + 1):np) {
    cormat[i, j] <- cormat[j, i] <- rgs[index]
    index <- index + 1
  }
}

## 抽样QTL效应
if (opt$nbin_cor > 0) { # nolint
  ## 遗传相关的列
  rg_cols <- c()
  for (i in 1:(np - 1)) {
    for (j in (i + 1):np) {
      rg_cols <- c(rg_cols, paste(fids[i], fids[j], sep = "_"))
    }
  }

  ## 存在遗传相关的标记效应
  if (opt$dist_cor == "identical") {
    qtl_cor <- which(qmap$bin %in% bins_sel)
    nsample <- length(qtl_cor)
    qmap[qtl_cor, fids] <- mvrnorm(nsample, rep(0, np), cormat)
    qmap[qtl_cor, rg_cols] <- matrix(rep(rgs, times = nsample), nsample, nrg, byrow = TRUE)
  } else {
    ## 各区间遗传相关
    if (opt$dist_cor == "normal") {
      rg_bins <- rnorm(opt$nbin_cor, mean = 0, sd = 0.5)
    } else if (opt$dist_cor == "uniform") {
      rg_bins <- runif(opt$nbin_cor, -1, 1)
    } else {
      cat("dist_cor can noly be identical, normal or uniform.\n")
      quit()
    }

    ## 抽样
    qmap[, fids] <- 0.00
    rgsi <- rep(0, nrg)
    for (i in seq_along(bins_sel)) {
      ## 标记效应协方差阵
      cormati <- cormat
      ## 赋值
      index <- 1
      for (j in 1:(np - 1)) {
        for (k in (j + 1):np) {
          if (cormat[j, k] == 0) next
          rgij <- rg_bins[i] - mean(rg_bins) + cormat[j, k]

          ## 保证合法性
          if (rgij > 1) rgij <- 1.0
          if (rgij < -1) rgij <- -1.0

          rgsi[index] <- cormati[k, j] <- cormati[j, k] <- rgij
          index <- index + 1
        }
      }

      ## 保证矩阵正定性
      cormat_pd <- nearPD(cormati, keepDiag = TRUE) # default
      if (!cormat_pd$converged) {
        stop("warning! The positive definiteness of genetic covariance matrix cannot be guaranteed.\n")
      } else if (cormat_pd$iterations > 1) {
        cat("add a small value to genetic effect covariance matrix\n")
        cormati <- cormat_pd$mat
      }

      ## 标记效应抽样
      qtl_cor <- which(qmap$bin == bins_sel[i])
      nsample <- length(qtl_cor)
      qmap[qtl_cor, fids] <- matrix(mvrnorm(nsample, rep(0, np), cormati), nsample, np)

      ## 保存相关系数
      qmap[qtl_cor, rg_cols] <- matrix(rep(rgsi, times = nsample), nsample, nrg, byrow = TRUE)
    }
  }

  ## 不存在遗传相关的标记效应
  qtl_no_cor <- which(qmap[[fids[1]]] == 0.0 | is.na(qmap[[fids[1]]]))
  if (length(qtl_no_cor) > 0) {
    qmap[qtl_no_cor, fids] <- mvrnorm(length(qtl_no_cor), rep(0, np), diag(np))
    qmap[qtl_no_cor, rg_cols] <- 0
  }
}

cat("global genetical correlation is:\n")
cor(subset(qmap, select = fids))

## 参考碱基(群体间一致)
gene <- qtl[[1]][, 7:ncol(qtl[[1]])]
ref_allele <- unlist(subset(gene, c(TRUE, rep(FALSE, nrow(gene) - 1)),
  select = seq(1, ncol(gene), 2)
))

## 生成表型信息
cat("Creating phenotypes files... \n")
pheno <- data.frame()
for (i in seq_len(np)) {
  gene <- qtl[[i]][, 7:ncol(qtl[[i]])]

  ## 如果基因型为字母表示，转换为1、2表示
  if (gene[1, 1] %in% LETTERS) {
    cat("Converting genotype data to 012 format...\n")

    ## 根据参考碱基转换为1、2标识符形式
    gene12 <- as.matrix(gene)
    for (j in 1:(ncol(gene) / 2)) {
      gene12[, 2 * j - 1] <- as.integer(subset(gene12, select = 2 * j - 1) == ref_allele[j]) + 1
      gene12[, 2 * j] <- as.integer(subset(gene12, select = 2 * j) == ref_allele[j]) + 1
    }

    ## 字符转数值
    gene <- apply(gene12, 2, as.numeric)
  }

  ## 基因含量转成012格式
  # cat("Converting genotype data to 012 format...\n")
  gene <- sapply(
    seq(1, ncol(gene) - 1, by = 2),
    function(i) rowSums(gene[, i:(i + 1)]) - 2
  )

  ## 计算加性效应
  ebvs[[i]]$tbv_raw <- apply(gene, 1, function(x) sum(x * qmap[[fids[i]]]))

  ## 校正为标准正态
  ebv_mean <- mean(ebvs[[i]]$tbv_raw)
  ebv_sd <- sd(ebvs[[i]]$tbv_raw)
  ebvs[[i]]$tbv_nor <- (ebvs[[i]]$tbv_raw - ebv_mean) / ebv_sd

  ## 残差效应(不相关)
  ve <- ((1 / h2s[i]) - 1) * var(ebvs[[i]]$tbv_nor)
  ebvs[[i]]$envir <- rnorm(nrow(ebvs[[i]]), 0, sqrt(ve))

  ## 加性 + 残差 + 群体均值 = 表型
  ebvs[[i]]$phe <- ebvs[[i]]$tbv_nor + ebvs[[i]]$envir + means[i]

  ## 合并各个品种的表型
  pheno <- rbind(pheno, ebvs[[i]])
}

## 增加品种信息和排序
if (is.null(opt$fid)) {
  pheno <- pheno[, c("iid", "breed", "phe", "tbv_raw", "tbv_nor", "envir")]
} else {
  ## 群(品种)效应标识符
  pheno$breed <- as.integer(as.factor(pheno$fid))
  pheno <- pheno[, c("iid", "fid", "breed", "phe", "tbv_raw", "tbv_nor", "envir")]
}

## 输出表型文件
if (is.null(opt$out)) opt$out <- paste(c(fids, "pheno.txt"), collapse = "_")
write.table(pheno, opt$out, row.names = FALSE, quote = FALSE)
cat("phenotypes output to file:", opt$out, "\n")

## 加性效应真值
if (is.null(opt$qtlf)) opt$qtlf <- paste(c(fids, "qtl.txt"), collapse = "_")
qmap <- qmap[order(qmap$bin), ]
write.table(qmap, opt$qtlf, row.names = FALSE, quote = FALSE, na = "0.0")
cat("true qtl effects output to file:", opt$qtlf, "\n")

## 标记与QTL是否重叠
if (is.null(opt$overlap) && !is.null(opt$bin)) {
  ## 更新每个区间内的标记数
  qtl_num <- as.data.frame(table(as.factor(qmap$bin)), stringsAsFactors = FALSE)
  qtl_num$Var1 <- as.integer(qtl_num$Var1)
  if (nrow(qtl_num) < nrow(bin)) {
    bins <- seq_len(nrow(bin))
    add_qtl_num <- data.frame(Var1 = bins[!bins %in% qmap$bin], Freq = 0)
    qtl_num <- bind_rows(qtl_num, add_qtl_num)
  }
  qtl_num <- qtl_num[order(qtl_num$Var1), ]
  bin$nsnp <- bin$nsnp - qtl_num$Freq
  newbinf <- paste0(basename(opt$bin), ".noqtl")
  write.table(bin$nsnp, newbinf, row.names = FALSE, quote = FALSE, col.names = FALSE)
  cat("new bins file output to:", newbinf, "\n")
}
