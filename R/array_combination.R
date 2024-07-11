#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript
## 字符串自由组合 ###
## 如输入：A B C
## 在不加其他限制的情况下会输出3种组合：'A B' 'A C' 'B C' 'A B C'

## 加载需要的程序包（未安装的话请提前安装）
#cat('Loading required packages... \n\n')
suppressPackageStartupMessages(library("getopt"))

## 获取命令行参数
spec <- matrix(
  c("array",  "a", 2, "character", "[Required] character array separated by space",
    "label",  "l", 2, "character", "[Optional] combination must include/exclude this string",
    "type",   "t", 2, "character", "[Optional] include or exclude [include]",
    "min",    "m", 2, "integer",   "[Optional] minimum number of elements contained in a combination",
    "max",    "M", 2, "integer",   "[Optional] maximum number of elements contained in a combination",
    "out",    "o", 2, "character", "[Optional] output file name",
    "append", "A", 2, "logical",   "[Optional] append to output file",
    "help",   "h", 0, "logical",   "This is Help!"),
  byrow = TRUE, ncol = 5)
opt <- getopt::getopt(spec = spec)

## 检查参数
if (!is.null(opt$help) || is.null(opt$array)) {
  cat(paste(getopt::getopt(spec = spec, usage = TRUE), "\n"))
  quit()
}

## 默认参数
if (is.null(opt$min)) opt$min <- 2
# if (is.null(opt$out)) opt$out <- "combination.txt"
if (is.null(opt$type)) opt$type <- "include"
if (!is.null(opt$append)) append <- TRUE else append <- FALSE

## 分割字符串
array <- strsplit(opt$array, " ")[[1]]
np <- length(array)
if (is.null(opt$max)) opt$max <- np

## 列出指定元素数量所有可能的组合的函数
com <- function(x, m) {
  if (m < 2) {
    return(x)
  }

  ## 自由组合
  df <- combn(x, m)

  ## 将每行变量粘贴成字符串
  result <- apply(df, 2, paste, collapse = " ")

  return(result)
}

## 组合
all_com <- c()
for (i in opt$min:opt$max) {
  com_i <- com(array, i)
  all_com <- c(com_i, all_com)
}

## 只保留有目标标签的组合
if (!is.null(opt$label)) {
  ## 判断字符串是否在向量中
  if (!opt$label %in% array) {
    cat("label not in array!\n")
    quit(1)
  }
  
  ## 找出含有某个标签的组合
  idx <- grepl(opt$label, all_com)

  ## 提取
  if (opt$type == "include") {
    all_com <- all_com[idx]
  } else {
    all_com <- all_com[!idx]
  }
}

## 输出文件
n_com <- length(all_com)
cat("The number of all possible combinations is:", n_com, "\n")
if (!is.null(opt$out)) {
  write.table(all_com, opt$out, quote = FALSE, col.names = FALSE, row.names = FALSE, append = append)
} else {
  print(all_com)
}
