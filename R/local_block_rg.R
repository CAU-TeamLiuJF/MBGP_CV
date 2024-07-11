#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript
## 根据提供的plink二进制文件估计LD独立区域 ##

## 命令行参数
spec <- matrix(
  c(
    "summ1",  "1", 1, "character", "[Required] GWAS descriptive statistics for trait A",
    "summ2",  "2", 1, "character", "[Required] GWAS descriptive statistics for trait B",
    "bfile",  "b", 1, "character", "[Required] plink binary file prefix",
    "block",  "g", 1, "character", "[Required] block define file",
    "nA",     "n", 1, "integer",   "[Required] individuals of breed A [traitA]",
    "nB",     "N", 1, "integer",   "[Required] individuals of breed B [traitA]",
    "traitA", "A", 1, "character", "[Optional] label of breed A [traitA]",
    "traitB", "B", 1, "character", "[Optional] label of breed B [traitB]",
    "univ",   "u", 1, "double",    "[Optional] P-value threshold for the univariate test [0.05]",
    "keep",   "k", 2, "character", "[Optional] keep intermediate files",
    "out",    "O", 1, "character", "[Optional] output file name prefix/NULL",
    "help",   "h", 0, "logical",   "This is Help!"
  ),
  byrow = TRUE, ncol = 5
)
opt <- getopt::getopt(spec = spec)

## 检查参数
needed <- c(opt$summ1, opt$summ2, opt$block, opt$nA, opt$nB, opt$bfile)
if (!is.null(opt$help) || length(needed) < 6) {
  cat(paste(getopt::getopt(spec = spec, usage = TRUERUE), "\n"))
  quit()
}

## 加载需要的程序包
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("LAVA"))
suppressPackageStartupMessages(library("dplyr"))

## 默认参数
if (is.null(opt$traitA)) opt$traitA <- "traitA"
if (is.null(opt$traitB)) opt$traitB <- "traitB"
if (is.null(opt$univ)) opt$univ <- 0.05
if (is.null(opt$out)) opt$out <- "ld_block_rg.txt"

## 检查文件是否存在
if (!all(file.exists(opt$summ1, opt$summ2, opt$block, paste0(opt$bfile, ".bim")))) {
  cat("One of the following files does not exist:\n")
  cat(opt$summ1, opt$summ2, opt$block, paste0(opt$bfile, ".bim"), sep = "\n")
  quit()
}

## 对提供的文件进行列重命名
summ1 <- fread(opt$summ1)
summ2 <- fread(opt$summ2)
block <- fread(opt$block)

names(block)[1:4] <- c("LOC", "CHR", "START", "STOP")

names_dict <- data.frame(
  new = c("CHR", "SNPID", "BP", "N", "A1", "A2", "MAF", "BETA", "SE", "P"),
  row.names = c("chr", "rs", "ps", "n_miss", "allele1", "allele0", "af", "beta", "se", "p_wald")
)

rename_col <- which(names(summ1) %in% rownames(names_dict))
names(summ1)[rename_col] <- names_dict[names(summ1)[rename_col], ]

rename_col <- which(names(summ2) %in% rownames(names_dict))
names(summ2)[rename_col] <- names_dict[names(summ2)[rename_col], ]

## 将文件中"N"列的缺失基因型个体改为非缺失基因型个体
summ1$N <- opt$nA - summ1$N
summ2$N <- opt$nB - summ2$N

## 写出文件
tmp_dir <- paste0(getwd(), "/rg_tmp", as.integer(runif(1) * 10000))
dir.create(tmp_dir)
write.table(summ1, file = paste0(tmp_dir, "/", opt$traitA, ".sumstats.txt"), row.names = FALSE, quote = FALSE)
write.table(summ2, file = paste0(tmp_dir, "/", opt$traitB, ".sumstats.txt"), row.names = FALSE, quote = FALSE)
write.table(block, file = paste0(tmp_dir, "/", opt$traitB, ".block"), row.names = FALSE, quote = FALSE)
info <- c(
  "phenotype cases controls filename",
  paste0(opt$traitA, " NA NA ", paste0(tmp_dir, "/", opt$traitA, ".sumstats.txt")),
  paste0(opt$traitB, " NA NA ", paste0(tmp_dir, "/", opt$traitB, ".sumstats.txt"))
)
writeLines(info, paste0(tmp_dir, "/input.info.txt"))

## 输入文件信息
input <- process.input(
  input.info.file = paste0(tmp_dir, "/input.info.txt"),
  sample.overlap.file = NULL,
  ref.prefix = opt$bfile,
  phenos = c(opt$traitA, opt$traitB)
)

### Read in locus info file
loci <- read.loci(opt$block)
nbin <- nrow(loci)

#### ANALYSE ####
print(paste("Starting LAVA analysis for", nbin, "loci..."))
results <- NULL
## 计算rg值
for (i in 1:nbin) {
  # process locus
  locus <- process.locus(loci[i, ], input)

  ## in some cases the locus object cannot be created due to
  ## e.g too few SNPs or negative variances in all analysed phenotypes, hence this check
  if (!is.null(locus)) {
    # extract general locus info for output
    loc_info <- data.frame(
      locus = locus$id, chr = locus$chr, start = locus$start,
      stop = locus$stop, n.snps = locus$n.snps, n.pcs = locus$K
    )

    # run univ results
    univ <- run.univ(locus)

    if (nrow(univ) > 1) {
      # run bivar analysis functions & store results
      bivar <- run.bivar(locus)
    } else {
      negative_trait <- setdiff(c(opt$traitA, opt$traitB), unlist(univ$phen))
      univ <- bind_rows(univ, data.frame(phen = negative_trait))
      bivar <- data.frame(phen1 = opt$traitA, phen2 = opt$traitB)
    }

    univ2 <- data.frame(
      h2.obs.phen1 = univ$h2.obs[1], h2.p.phen1 = univ$p[1],
      h2.obs.phen2 = univ$h2.obs[2], h2.p.phen2 = univ$p[2]
    )
    bivar2 <- bind_cols(bivar, univ2)

    results <- bind_rows(results, cbind(loc_info, bivar2))
  }
}

## 将NA设为0
results$rho[is.na(results$rho)] <- 0

## 将两个性状都在该区间显示出显著遗传力的区间打上标记
results$h2.p.phen1[is.na(results$h2.p.phen1)] <- 1.0
results$h2.p.phen2[is.na(results$h2.p.phen2)] <- 1.0
phen1_sig <- results$h2.p.phen1 < opt$univ
phen2_sig <- results$h2.p.phen2 < opt$univ
results$h2_sig <- "FALSE"
results$h2_sig[phen1_sig & phen2_sig] <- "TRUE"
cat("A total of", sum(phen1_sig), "blocks with significant heritability were identified in", opt$traitA, "\n")
cat("A total of", sum(phen2_sig), "blocks with significant heritability were identified in", opt$traitB, "\n")

## 鉴定到显著rg的区间打上标签
results$p[is.na(results$p)] <- 1.0
rg_sig <- results$p < opt$univ
results$hg_sig <- "FALSE"
results$hg_sig[rg_sig] <- "TRUE"

## 输出文件
cat("A total of", sum(rg_sig), "blocks with significant genetic correlation were identified.\n")
write.table(results, opt$out, row.names = FALSE, quote = FALSE, sep = "\t")

## 删除中间文件
if (is.null(opt$keep)) {
  cat(tmp_dir, "has been deleted.\n")
  unlink(tmp_dir, recursive = TRUE)
}
