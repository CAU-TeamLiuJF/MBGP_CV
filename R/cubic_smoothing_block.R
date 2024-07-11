#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript
## 根据自编的C语言程序计算的mean R2值进行区域划分
# debug
# opt <- list(files="breedAsq breedBsq", mapf = "breedAsq.map")

## 加载需要的程序包
cat("Loading required packages... \n\n")
suppressPackageStartupMessages(library("getopt"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("Cairo"))

## 参数列表
spec <- matrix(
  c("r2",   "I", 1, "character", "[Required] R2 mean\n",
    "spar", "S", 1, "double",    "[Optional] Smoothness (0-1) when fitting curves [0.2]\n",
    "diff", "D", 1, "double",    "[Optional] [0.05]\n",
    "bim",  "B", 1, "character", "[Optional] bim file name [NULL]\n",
    "out",  "O", 1, "character", "[Optional] Output file name [cubic_LD_block.txt]\n",
    "plot", "P", 0, "logical",   "[Optional] Output plot [FALSE]\n",
    "help", "h", 0, "logical",  "This is Help!"),
  byrow = TRUE, ncol = 5
)
opt <- getopt(spec = spec)

## 检查参数
if (!is.null(opt$help) || is.null(opt$r2)) {
  cat(paste(getopt(spec = spec, usage = TRUERUE), "\n"))
  quit()
}

## 默认参数
if (is.null(opt$out)) opt$out <- "cubic_LD_block.txt"
if (is.null(opt$spar)) opt$spar <- 0.2
if (is.null(opt$diff)) opt$diff <- 0.05

## 函数定义
find_local_min <- function(x, diff = 0.2) { # nolint
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

    if (left_peak && right_peak) final <- c(final, i)
  }

  ## 端点处不要
  final <- final[final != 1 & final != n]

  return(final)
}

## 读取文件
r2 <- fread(opt$r2)

## 拟合曲线
fit <- smooth.spline(r2$V1, spar = 0.2)

## 搜索符合要求的极低值点
point <- find_local_min(fit$y, 0.05)

## 计算每个区间的SNP数
point2 <- c(0, point, length(fit$y))
n <- length(point2)
nsnp <- point2[2:n] - point2[1:(n - 1)]
nsnp[length(nsnp)] <- nsnp[length(nsnp)] + 1

## SNP数目检查
if (sum(nsnp) != (nrow(r2) + 1)) {
  cat("The number of SNPs in the result is not equal to that in the input file!\n")
  quit()
}

## 输出中包含其他信息
if (!is.null(opt$bim)) {
  if (file.exists(opt$bim)) {
    bim <- fread(opt$bim)
    if (nrow(bim) == sum(nsnp)) {
      # CHR START STOP nSNP
      point2[1] <- 1
      nsnp <- data.frame(CHR = bim$V1[1], START = bim$V4[point2[1:(n - 1)]], STOP = bim$V4[point2[2:n]], nSNP = nsnp)
    } else {
      cat("The number of SNPs in the result is not equal to that in the bim file!\n")
      quit()
    }
  } else {
    cat(opt$bim, "not found!\n")
    quit()
  }
}

## 输出SNP数目文件
write.table(nsnp, opt$out, row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("bins definition file has been output to:", opt$out, "\n")

## 作图
if (!is.null(opt$plot)) {
  CairoPNG(paste0(opt$out, ".png"), width = 1600, height = 800)
  plot(r2$V1, col = "grey", pch = 20,
       ylab = "local mean r2",
       xlab = "SNP rank")
  lines(fit, col = "red", lwd = 2)
  for (h in point) {
    abline(v = h, col = "blue", lwd = 1)
  }
  dev.off()
  cat("plot has been output to:", paste0(opt$out, ".png"), "\n")
}
