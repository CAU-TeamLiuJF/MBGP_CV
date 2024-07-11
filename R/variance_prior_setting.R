#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript

## liwn@cau.edu.cn 2023-07-05
## 根据表型方差估计近似遗传和残差方差，用作初始值

# 加载需要的程序包
suppressPackageStartupMessages(library("getopt"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("Matrix"))

## 命令行参数
spec <- matrix(
  c("filea",    "p", 1, "character", "[Required] phenotype file 1",
    "fileb",    "P", 1, "character", "[Optional] phenotype file 2",
    "pcol",     "c", 1, "double",    "[Optional] phenotype column in phenotype file [first real type column]",
    "type",     "t", 1, "character", "[Optional] Format of the output variance file. blend/union/multi [union]",
    "rep",      "r", 1, "integer",   "[Optional] repeat times [1]",
    "fold",     "f", 1, "integer",   "[Optional] cross validation fold [1]",
    "var",      "v", 1, "character", "[Optional] pheno/predict [pheno]",
    "h2",       "a", 1, "double",    "[Optional] heritability [0.5]",
    "h2B",      "A", 1, "double",    "[Optional] heritability of trait 2 [h2]",
    "rg",       "g", 1, "double",    "[Optional] genetic correlation [0.001]",
    "re",       "e", 1, "double",    "[Optional] residual correlation [0.001]",
    "add_rf",   "d", 1, "integer",   "[Optional] additive effect random group in dmu [1]",
    "rg_local", "l", 1, "character", "[Optional] file contains correlations of each bins [1]",
    "miss",     "n", 1, "double",    "[Optional] missing value [-99]",
    "out",      "o", 1, "character", "[Optional] output file name",
    "overwri",  "O", 2, "character", "[Optional] overwrie the existing result file",
    "norec",    "R", 0, "logical",   "[Optional] non-existing covariances between residuals",
    "help",     "h", 0, "logical",   "This is Help!"),
  byrow = TRUE, ncol = 5)
opt <- getopt(spec = spec)

if (!is.null(opt$help) || is.null(opt$filea)) {
  cat(paste(getopt(spec = spec, usage = TRUE), "\n"))
  quit()
}

## 默认参数
if (is.null(opt$h2)) opt$h2 <- 0.5
if (is.null(opt$h2B)) opt$h2B <- opt$h2
if (is.null(opt$rg)) opt$rg <- 0.001
if (is.null(opt$re)) opt$re <- 0.001
if (is.null(opt$miss)) opt$miss <- -99
# if (is.null(opt$out)) opt$out = "var_prior.txt"
if (is.null(opt$type)) opt$type <- "union"
if (is.null(opt$add_rf)) opt$add_rf <- 1
if (is.null(opt$var)) opt$var <- "pheno"
if (is.null(opt$fileb)) {
  opt$fileb <- opt$filea
  cat("warn: fileb not set! set to filea\n")
}

output <- FALSE
for (r in 1:opt$rep) { # nolint
  for (f in 1:opt$fold) {
    ## 更换路径中的占位符
    filea <- gsub("#val#", f, gsub("#rep#", r, opt$filea))
    fileb <- gsub("#val#", f, gsub("#rep#", r, opt$fileb))
    out <- gsub("#val#", f, gsub("#rep#", r, opt$out))

    ## 检查文件夹中是否已有表型文件
    if (file.exists(out) && opt$overwri != "true") next

    output <- TRUE

    ## 检查文件是否存在
    if (!file.exists(filea) || !file.exists(fileb)) {
      cat("warn: file a or b not found, prior information of rep", r, "fold", f, "will not be generated!\n")
      next
    }

    ## 读取文件
    data1 <- fread(filea)
    data2 <- fread(fileb)

    ## 获取方差组分
    if (opt$var == "pheno") {
      ## 命名
      if (is.null(opt$pcol)) {
        int_index <- apply(data1, 2, function(x) all(floor(x) == x))
        opt$pcol <- which(!int_index)[1]
      } else if (opt$pcol > ncol(data1)) {
        stop("pcol (", opt$pcol, ") cannot be bigger than file columns (", ncol(data1), ")")
      }
      names(data1)[opt$pcol] <- names(data2)[opt$pcol] <- "phe"
      ## 缺失值处理
      data1$phe[data1$phe == opt$miss] <- NA
      data2$phe[data2$phe == opt$miss] <- NA
      ## 计算表型方差
      pvar1 <- var(data1$phe, na.rm = TRUE)
      pvar2 <- var(data2$phe, na.rm = TRUE)
      ## 获取方差组分
      vara1 <- pvar1 * opt$h2
      vara2 <- pvar2 * opt$h2B
      vare1 <- pvar1 * (1 - opt$h2)
      vare2 <- pvar2 * (1 - opt$h2B)
    } else if (opt$var == "predict") {
      if (opt$add_rf >= max(data1$V1)) stop("add_rf cannot be larger than the number of random groups in dmu!")
      if (any(data1$V2 != data2$V2)) stop("The random effects groups are different in the two files!")
      ## 获取方差组分
      vara1 <- data1$V4[opt$add_rf]
      vara2 <- data2$V4[opt$add_rf]
      vare1 <- tail(data1$V4, 1)
      vare2 <- tail(data2$V4, 1)
    }

    ## 计算协方差
    cova <- opt$rg * sqrt(vara1) * sqrt(vara2)
    cove <- opt$re * sqrt(vare1) * sqrt(vare2)

    ## 保证协方差矩阵正定性
    ## 遗传方差
    vara <- matrix(c(vara1, cova, cova, vara2), 2, 2)
    vara_pd <- nearPD(vara, keepDiag = TRUE) # default
    if (!vara_pd$converged) {
      stop("warning! The positive definiteness of genetic covariance matrix cannot be guaranteed.\n")
    } else if (vara_pd$iterations > 1) {
      cat("add a small value to genetic effect covariance matrix\n")
      cova <- vara_pd$mat[2, 2]
    }
    ## 残差方差
    vare <- Matrix(c(vare1, cova, cova, vare2), 2, 2)
    vare_pd <- nearPD(vare, keepDiag = TRUE) # default
    if (!vare_pd$converged) {
      stop("warning! The positive definiteness of genetic covariance matrix cannot be guaranteed.\n")
    } else if (vare_pd$iterations > 1) {
      cat("add a small value to genetic effect covariance matrix\n")
      cova <- vare_pd$mat[2, 2]
    }

    ## 方差组分形式
    if (opt$type == "union") {
      prior <- data.frame(
        group = c(1, 1, 1, 2, 2, 2),
        rindx = c(1, 2, 2, 1, 2, 2),
        cindx = c(1, 1, 2, 1, 1, 2),
        var = c(vara1, cova, vara2, vare1, cove, vare2)
      )
      ## 约束残差协方差
      if (!is.null(opt$norec)) {
        prior <- prior[-5, ]
      }
    } else if (opt$type == "multi") {
      if (!is.null(opt$rg_local)) {
        ## 检查文件是否存在
        if (!file.exists(opt$rg_local)) {
          cat(opt$rg_local, "not found!\n")
          quit(status = 1)
        }

        ## 获取相关系数
        rg_local <- fread(opt$rg_local)
        names(rg_local)[ncol(rg_local)] <- "cor"
        names(rg_local)[3] <- "nsnp"
        nsnp_total <- sum(rg_local$nsnp)
        ## 计算每个区间的方差组分
        prior <- data.frame(
          var1 = vara1 * rg_local$nsnp / nsnp_total,
          var2 = vara2 * rg_local$nsnp / nsnp_total
        )
        prior$cov1 <- prior$cov2 <- rg_local$cor * sqrt(prior$var1) * sqrt(prior$var2)
        prior <- prior[, c("var1", "cov1", "cov2", "var2")]
        ## 加上残差方差组分
        prior[nrow(prior) + 1, ] <- c(vare1, cove, cove, vare2)
        ## 保证矩阵的正定性
        for (i in seq_len(nrow(prior))) {
          vari <- matrix(unlist(prior[i, ]), 2, 2)
          vari_pd <- nearPD(vari, keepDiag = TRUE) # default
          if (!vari_pd$converged) {
            cat("warning! The positive definiteness of bin ", i, " matrix cannot be guaranteed.\n")
          } else if (vari_pd$iterations > 1) {
            cat("The covariance matrix of bin ", i, " is set to: \n", vari_pd$mat, "\n")
            prior[i, c("var1", "cov1", "cov2", "var2")] <- as.vector(vari_pd$mat)
          }
        }
      } else {
        prior <- data.frame(matrix(c(
          vara1, cova, cova, vara2,
          vare1, cove, cove, vare2
        ), 2, 4, byrow = TRUE))
      }
    } else if (opt$type == "blend") {
      if (opt$var == "predict") {
        prior <- data1
        prior$V4 <- (prior$V4 + data2$V4) / 2
      } else {
        prior <- data.frame(
          group = c(1, 2),
          rindx = c(1, 1),
          cindx = c(1, 1),
          var = c(mean(vara1, vara2), mean(vare1, vare2))
        )
      }
    } else {
      stop("Unknown type!")
    }

    if (!is.null(opt$out)) {
      ## 输出文件
      write.table(prior, out, row.names = FALSE, col.names = FALSE, quote = FALSE)
      # cat("prior file output to:", out, "\n")
    } else {
      prior
    }
  }
}

if (output) cat("variance component file generated.\n")
