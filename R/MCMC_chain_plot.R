#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript
## MCMC链 作图
# debug
# setwd('/BIGDATA2/cau_jfliu_2/liwn/mbGS/Real/YCJY/rmodel/AGE/multi/val1/rep1/YY_fix_0.1_100')
# opt <- list(files="MCMC_samples_heritability.txt", start=50010, end=70000, thin=10)

## 命令行参数
spec <- matrix(
  c("files",   "F", 2, "character", "[Required] file(s) contain a vector, or a matrix with one column per variable",
    "start",   "S", 2, "integer",   "[Required] the iteration number of the first observation",
    "end",     "E", 2, "integer",   "[Required] the iteration number of the last observation",
    "thin",    "T", 2, "integer",   "[Required] the thinning interval between consecutive observations",
    "names",   "N", 2, "character", "[Optional] Name(s) of each variable marked in output figure/var*",
    "out",     "O", 2, "character", "[Optional] Output file name prefix/MCMC_*_chain_xy(/density).png",
    "help",    "h", 0, "logical",   "This is Help!"),
  byrow = TRUE, ncol = 5)
opt <- getopt::getopt(spec = spec)

## 检查参数
if (!is.null(opt$help) || is.null(opt$files)) {
  cat(paste(getopt::getopt(spec = spec, usage = TRUE), "\n"))
  quit()
}

# 加载需要的程序包
cat("Loading required packages... \n\n")
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("lattice"))
suppressPackageStartupMessages(library("coda"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("Cairo"))


## 默认参数
if (is.null(opt$start)) opt$win <- 50000
if (is.null(opt$end))   opt$slid <- 70000
if (is.null(opt$thin))  opt$seg <- 10

## 文件
files <- unlist(strsplit(opt$files, " "))

for (i in seq_along(files)) {
  mc_i <- fread(files[i])

  ## 定义变量的名称
  if (is.null(opt$names)) {
    opt$names <- paste0("var", seq_len(ncol(mc_i)))
  } else {
    opt$names <- unlist(strsplit(opt$names, " "))
    if (length(opt$names) != ncol(mc_i)) {
      cat("The supplied variable names does not match the number of columns in the file.\n")
      quit()
    }
  }

  ## 命名
  names(mc_i) <- opt$names

  ## 转化成mcmc对象
  mc_i <- mcmc(mc_i, start = opt$start, end = opt$end, thin = opt$thin)

  ## 合并不同链
  if (i == 1) {
    mc_list <- mcmc.list(mc_i)
  } else {
    mc_list[i] <- mcmc.list(mc_i)
  }
}

## 链数量
num_chain <- length(mc_list)

## 输出文件名前缀
if (is.null(opt$out)) opt$out <- paste0("MCMC_", num_chain, "_chain")

## 折线图 xyplot
filename <- paste0(opt$out, "_xy.png")
CairoPNG(filename, width = 1200, height = 900, bg = "white")
# png(filename, width=1200, height=900, bg = "white")
xyplot(mc_list,
  scales = list(tck = c(1, 0), x = list(cex = 2.5), y = list(cex = 1.4)),
  xlab = list(label = "Iteration number", cex = 2.5),
  lwd = 2,
  strip = strip.custom(
    bg = "lightgrey",
    par.strip.text = list(
      color = "black",
      cex = 2,
      font = 3
    )
  ),
  par.settings = list(layout.heights = list(strip = 2))
)
hide_message <- dev.off()

## 密度图 densityplot
filename <- paste0(opt$out, "_density.png")
CairoPNG(filename, width = 1200, height = 900, bg = "white")
# png(filename, width=1200, height=900, bg = "white")
densityplot(mc_list,
  scales = list(tck = c(1, 0), x = list(cex = 2.5), y = list(cex = 1.4)),
  xlab = list(label = "Iteration number", cex = 2.5),
  lwd = 2,
  strip = strip.custom(
    bg = "lightgrey",
    par.strip.text = list(
      color = "black",
      cex = 2,
      font = 3
    )
  ),
  par.settings = list(layout.heights = list(strip = 2))
)
hide_message <- dev.off()

cat("MCMC plot finished.\n")
