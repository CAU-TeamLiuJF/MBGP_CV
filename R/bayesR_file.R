#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript

########################################################################################################################
## Version: 1.3.0
## Author:    Liweining liwn@cau.edu.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-08-20
##
## Function：
## Generate plink genotype (fam) files and fixed effects structure matrix files required for BayesR
##
##
## Usage: ./bayesR_file.R --bfile "/plink/binary/file/prefix" ...(Please refer to --help for detailed parameters)
##
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


# load packages
suppressPackageStartupMessages(library("getopt"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dplyr"))

## Command Line Parameters
spec <- matrix(
  c("bfile",   "b", 1, "character", "[Required] PLINK binary file prefix",
    "phe",     "p", 1, "character", "[Required] phenotype file name",
    "fix",     "f", 1, "character", "[Optional] fix effect(s) column(s), e.g. '2,3' [NULL]",
    "phe_col", "P", 1, "integer",   "[Optional] phenotype column location [last one]",
    "miss",    "m", 1, "integer",   "[Optional] missing phenotype value [-99]",
    "gt_out",  "o", 1, "character", "[Optional] output PLINK binary file prefix [overwrite]",
    "fix_out", "F", 1, "character", "[Optional] output file name [fix_BayesR.txt]",
    "out_dir", "D", 1, "character", "[Optional] output file directory [NULL]",
    "help",    "h", 0, "logical",   "This is Help!"),
  byrow = TRUE, ncol = 5)
opt <- getopt(spec = spec)

## check parameters
if (is.null(opt$bfile) || is.null(opt$phe)) {
  cat(paste(getopt(spec = spec, usage = TRUE), "\n"))
}

## default parameters
if (is.null(opt$miss)) opt$miss <- -99

## load phenotype and genotype
phe <- fread(opt$phe)
fam <- fread(paste0(opt$bfile, ".fam"))

## output directory
if (!is.null(opt$out_dir)) {
  setwd(opt$out_dir)
}

## order of individuals
if (!all(fam$V2 == phe$V1)) {
  if (!all(nrow(fam) %in% nrow(phe))) {
    cat("Error: The number of individuals in genotype and phenotype files is not equal! \n")
    quit(status = -1)
  } else {
    phe2 <- subset(phe, select = c(1, opt$phe_col))
    names(phe2) <- c("V2", "phe")
    fam2 <- left_join(fam, phe2, by = "V2")
    fam <- fam2[, -6]
  }
} else {
  fam$V6 <- subset(phe, select = opt$phe_col)
}

## prepare phenotype (fam) file
names(fam)[6] <- "phe"
fam$phe[fam$phe <= opt$miss] <- "NA"
if (!is.null(opt$gt_out)) {
  write.table(fam, paste0(opt$gt_out, ".fam"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  file.copy(paste0(opt$bfile, ".bed"), paste0(opt$gt_out, ".bed"), overwrite = TRUE)
  file.copy(paste0(opt$bfile, ".bim"), paste0(opt$gt_out, ".bim"), overwrite = TRUE)
  cat("The genotype files has been output to:", paste0(opt$gt_out, ".*"), "\n")
} else {
  ## overwrite
  write.table(fam, paste0(opt$bfile, ".fam"), row.names = FALSE, col.names = FALSE, quote = FALSE)
}

## prepare fix effect file
if (!is.null(opt$out_dir)) {
  cols <- as.numeric(unlist(strsplit(paste0(opt$fix, " "), "[, ]+")))
  i <- 1
  for (j in seq_along(cols)) {
    if (length(unique(unlist(subset(phe, select = cols[j])))) == 1) next  ## ignore intercept
    ci <- model.matrix(~as.factor(unlist(subset(phe, select = cols[j]))))

    if (i == 1) {
      fix <- data.frame(V1 = subset(ci, select = -1))
    } else {
      fix <- cbind(fix, subset(ci, select = -1))
    }
    i <- i + 1
  }

  if (i > 1) {
    write.table(fix, opt$fix_out, row.names = FALSE, col.names = FALSE, quote = FALSE)
    cat("The fixed effects design matrix has been output to:", opt$fix_out, "\n")
  } else {
    cat("There are no other fixed effects except for the intercept term.\n")
  }
}
