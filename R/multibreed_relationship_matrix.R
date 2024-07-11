#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript
# 李伟宁 2021-12-08
## 计算多群体关系矩阵 ###
## 输入两个群体合并基因型文件(plink的012格式，*.raw) ###
## 输出不同方法计算的多群体关系矩阵[(ida+idb) x (ida+idb)] ###
## 输出格式为三列，id1 id2 value，tab分隔符
## 使用前需要安装相应的包，以及调用c++时需要的组件
## debug
# opt = list(rawf='popA_popB.raw', idaf='popA.ids', idbf='popB.ids')

# Load packages
cat("Loading packages needed...\n")
suppressPackageStartupMessages(library("getopt"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("Rcpp"))
suppressPackageStartupMessages(library("RcppArmadillo"))
suppressPackageStartupMessages(library("RcppEigen"))

## command parameters
spec <- matrix(
  c("rawf",    "R", 2, "character", "[Required] Raw file after merging two populations",
    "idaf",    "A", 2, "character", "[Required] ID file of population A with only 1 column",
    "idbf",    "B", 2, "character", "[Required] ID file of population B with only 1 column",
    "phef",    "P", 2, "character", "[Optional] If columns 1 and 2 in G matrix are row/column index,
                                                the ID in the phenotype file is renumbered",
    "method",  "M", 2, "character", "[Optional] Method of calculating G matrix/[Yvonne(default)/No/Chen/Mean]",
    "workdir", "D", 2, "character", "[Optional] working directory",
    "out",     "O", 2, "character", "[Optional] output file name(Can include path)",
    "tol",     "T", 0, "double",    "[Optional] tol of make matrix positive define",
    "inv",     "V", 0, "logical",   "[Optional] Output inverse matrix",
    "logdet",  "L", 0, "logical",   "[Optional] Add a row of '0 0 log(det)' (for DMU)",
    "header",  "H", 0, "logical",   "[Optional] include header in output files/NULL",
    "index",   "I", 0, "logical",   "[Optional] The first two columns of the output file are line index\n",
    "help",    "h", 0, "logical",   "This is Help!"),
  byrow = TRUE, ncol = 5
)
opt <- getopt::getopt(spec = spec)

## 检查必要参数
if (!is.null(opt$help) || is.null(opt$rawf) || is.null(opt$idaf) || is.null(opt$idbf)) {
  ## print help message
  cat(paste(getopt::getopt(spec = spec, usage = TRUERUE), "\n"))
  quit(status = -1)
}

## 默认参数设置
if (is.null(opt$header)) header <- FALSE else header <- TRUE
if (is.null(opt$tol)) opt$tol <- 1e-6
if (is.null(opt$out)) opt$out <- "Yvonne2017"
if (is.null(opt$method)) opt$method <- "Yvonne"

## 检查method参数
if (!opt$method %in% c("Yvonne", "Chen", "Mean", "No")) {
  cat("Method can only have one of these options: Yvonne/No/Chen/Mean\n")
  quit(status = -1)
}

## 更换工作路径
if (!is.null(opt$workdir)) setwd(opt$workdir)

## 准备加速矩阵乘法的c++脚本
cpp <- "// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]"
cpp <- c(cpp, "#include <RcppArmadillo.h>")
cpp <- c(cpp, "#include <RcppEigen.h>")
cpp <- c(cpp, "// [[Rcpp::export]]")
cpp <- c(cpp, "SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){")
cpp <- c(cpp, "    Eigen::MatrixXd C = A * B;")
cpp <- c(cpp, "    return Rcpp::wrap(C);")
cpp <- c(cpp, "}")
cpp <- c(cpp, "// [[Rcpp::export]]")
cpp <- c(cpp, "arma::mat getInv(arma::mat M) {")
cpp <- c(cpp, "    return arma::inv(M);")
cpp <- c(cpp, "}")

write.table(cpp, "Matrix_multiplication.cpp",
  sep = "\n", col.names = FALSE,
  row.names = FALSE, quote = FALSE
)

## 加载c++函数
cat("Loading C++ functions...\n")
sourceCpp("Matrix_multiplication.cpp")

## 读取基因型文件
cat("Loading genotype...\n")
raw <- fread(opt$rawf)

## 获取不同群体id
ida <- read.table(opt$idaf)
idb <- read.table(opt$idbf)

## 群体个体数
num_a <- nrow(ida)
num_b <- nrow(idb)
num_all <- num_a + num_b

## 不同群体基因含量矩阵
geno_a <- as.matrix(raw[raw$IID %in% ida$V1, -c(1:6)])
geno_b <- as.matrix(raw[raw$IID %in% idb$V1, -c(1:6)])

## 标记数
nsnp_a <- ncol(geno_a)
nsnp_b <- ncol(geno_b)
nsnp_all <- ncol(geno_a)

## 基因频率
cat("Calculating gene frequency...\n")
af_a <- apply(geno_a, 2, sum) / num_a / 2
sum2pq_a <- sum(2 * af_a * (1 - af_a))

af_b <- apply(geno_b, 2, sum) / num_b / 2
sum2pq_b <- sum(2 * af_b * (1 - af_b))

sum2pq_ab <- 2 * sum(sqrt(af_a * (1 - af_a) * af_b * (1 - af_b)))

af_all <- apply(raw[, -c(1:6)], 2, sum) / num_all / 2
sum2pq_all <- sum(2 * af_all * (1 - af_all))

## Z-matrices
pa <- matrix(2 * af_a, byrow = TRUE, nrow = num_a, ncol = nsnp_a)
za <- geno_a - pa

pb <- matrix(2 * af_b, byrow = TRUE, nrow = num_b, ncol = nsnp_b)
zb <- geno_b - pb

pall <- matrix(2 * af_all, byrow = TRUE, nrow = num_all, ncol = nsnp_all)
zall <- as.matrix(raw[, -c(1:6)]) - pall

## 删除变量，释放内存
rm(geno_a)
rm(geno_b)
rm(raw)
gclog <- gc(verbose = FALSE)

## 关系矩阵容器
gmat <- matrix(0, nrow = num_all, ncol = num_all)

## 群体在Gmat中的位置索引
ida_index <- c(rep(TRUE, num_a), rep(FaLSE, num_b))

## 基因组关系矩阵校正因子
scalea <- sum2pq_a
scaleb <- sum2pq_b
if (opt$method == "Yvonne") {
  scaleab <- sqrt(sum2pq_a) * sqrt(sum2pq_b)
} else if (opt$method == "Chen") {
  scaleab <- sum2pq_ab
} else if (opt$method == "Mean") {
  scalea <- scaleb <- scaleab <- sum2pq_all
} else if (opt$method == "No") {
  scalea <- scaleb <- scaleab <- 1
}

## 生成G矩阵
cat("Generating G matrix...\n")
if (opt$method != "Mean") {
  gmat[ida_index, ida_index] <- eigenMapMatMult(za, t(za)) / scalea
  gmat[!ida_index, !ida_index] <- eigenMapMatMult(zb, t(zb)) / scaleb
  gmat[!ida_index, ida_index] <- eigenMapMatMult(zb, t(za)) / scaleab
  gmat[ida_index, !ida_index] <- eigenMapMatMult(za, t(zb)) / scaleab
} else {
  gmat[, ] <- eigenMapMatMult(zall, t(zall)) / sum2pq_all
}

## 保证矩阵正定
make_positive_definite <- function(mat, tol = 1e-6) {
  eig <- eigen(mat, symmetric = TRUE)
  rtol <- tol * eig$values[1]
  if (min(eig$values) < rtol) {
    cat("Make the matrix positive definite...\n")
    vals <- eig$values
    vals[vals < rtol] <- rtol
    srev <- eigenMapMatMult(eig$vectors, vals * t(eig$vectors))
    dimnames(srev) <- dimnames(mat)
    return(srev)
  } else {
    return(mat)
  }
}

## 保证矩阵正定
gmat_pd <- make_positive_definite(gmat, tol = opt$tol)

## 是否输出逆矩阵
if (!is.null(opt$inv)) {
  cat("Inverting matrix...\n")
  opt$out <- paste0(opt$out, ".grm.inv")
  gmat_out <- getInv(gmat_pd)
} else {
  opt$out <- paste0(opt$out, ".grm")
  gmat_out <- gmat_pd
}

## 是否以id命名
if (is.null(opt$index)) {
  ## 输出文件中1、2列为个体基因型文件中id
  rownames(gmat_out) <- colnames(gmat_out) <- c(ida$V1, idb$V1)
} else {
  if (!is.null(opt$phef)) {
    ## 输出的基因型矩阵中1、2列为行/列号，对表型文件中id进行重编
    phe <- fread(opt$phef) ## 默认第一列为id
    names(phe)[1] <- "org"
    mtab <- data.frame(org = c(ida$V1, idb$V1), new = 1:num_all)
    pcol <- ncol(phe)
    phe2 <- merge(phe, mtab, by = "org", sort = FALSE)
    phe <- subset(phe2, select = c(pcol + 1, 2:pcol, 1))
    fwrite(phe, paste0(basename(opt$phef), ".new"), col.names = FALSE, sep = " ")
    fwrite(mtab, "grm_id_match.txt", sep = " ")
    cat("The ID in phenotype file has been replaced\n")
  }
}

## 将矩阵转换成三列形式(下三角)
gmat_out[upper.tri(gmat_out)] <- NA
gmat3col <- melt(gmat_out, na.rm = TRUE)

## 是否输出矩阵的log(determinant)
if (!is.null(opt$logdet)) {
  cat("Calculating determinant of matrix...\n")
  logdet <- determinant(gmat_pd)$modulus
  row1add <- matrix(c(0, 0, logdet), nrow = 1, ncol = 3)
  if (is.na(row1add[, 3])) {
    cat("The determinant of matrix is NA, and log(det) is not output\n")
  } else {
    colnames(row1add) <- colnames(gmat3col)
    gmat3col <- rbind(row1add, gmat3col)
  }
}

## 写出矩阵
cat("Writing out matrix...\n")
fwrite(gmat3col, opt$out, sep = " ", col.names = header)

cat("G matrix has been exported to: ", opt$out, "\n")
