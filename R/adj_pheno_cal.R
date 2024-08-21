#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript

########################################################################################################################
## Version: 1.3.0
## Author:    Liweining liwn@cau.edu.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-08-20
##
## Function：
## Calculate the corrected phenotype (corrected for fixed effects and non additive effects) based on the
## DMU result file and phenotype file
##
##
## Usage: ./adj_pheno_cal.R --DIR "/path/to/DIR/file/prefix" ...(Please refer to --help for detailed parameters)
##
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


## Load packages
suppressPackageStartupMessages(library("getopt"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dplyr"))

## Command-line parameters
spec <- matrix(
  c("DIR",     "d", 1, "character", "[Required] Full dataset parameter card prefix",
    "phe",     "p", 1, "character", "[Required] Full dataset phenotype file name",
    "idc",     "c", 1, "integer",   "[Optional] id columns in phenotype file [1]",
    "traiti",  "i", 1, "integer",   "[Optional] trait rank [1]",
    "nTrait",  "n", 1, "integer",   "[Optional] number of traits [1]",
    "add_sol", "a", 1, "integer",   "[Optional] number of traits [1]",
    "append",  "A", 1, "character", "[Optional] append to the result file\n",
    "out",     "o", 1, "character", "[Optional] output file name/phe_adj.txt\n",
    "help",    "h", 0, "logical",   "This is Help!"),
  byrow = TRUE, ncol = 5)
opt <- getopt(spec = spec)

## Check parameters
if (!is.null(opt$help) || is.null(opt$phe) || is.null(opt$DIR)) {
  cat(paste(getopt(spec = spec, usage = TRUE), "\n"))
}

## Check if the required files exist
ex <- TRUE
if (!file.exists(opt$phe)) {
  file_name <- opt$phe
} else if (!file.exists(paste0(opt$DIR, ".SOL"))) {
  file_name <- paste0(opt$DIR, ".SOL")
} else if (!file.exists(paste0(opt$DIR, ".RESIDUAL"))) {
  file_name <- paste0(opt$DIR, ".RESIDUAL")
} else {
  ex <- FALSE
}
if (ex) {
  cat(file_name, "not found!\n")
  quit(status = 1)
}

## Default parameters
if (is.null(opt$traiti)) opt$traiti <- 1
if (is.null(opt$nTrait)) opt$nTrait <- 1
if (is.null(opt$idc)) opt$idc <- 1
if (is.null(opt$out)) opt$out <- "phe_adj.txt"

## Phenotype file, used to match residual IDs
phe <- fread(opt$phe)
names(phe)[opt$idc] <- "id"

## ebv
sol <- fread(paste0(opt$DIR, ".SOL"))
if (is.null(opt$add_sol)) opt$add_sol <- max(sol$V1)
if (!opt$add_sol %in% unique(sol$V1)) {
  cat("add_sol cant be ", opt$add_sol, "!\n")
  quit(status = 1)
}
ebv_all <- subset(sol, V1 == opt$add_sol & V2 == opt$traiti, c(5, 8))
names(ebv_all) <- c("id", "ebv")

## Extract residuals (currently only applicable to single-trait models)
res <- fread(paste0(opt$DIR, ".RESIDUAL"))
re <- subset(res, select = c(1, 4))
names(re) <- c("rows", "re")
if (max(re$rows) > nrow(phe)) {
  cat("sol file does not match RESIDUAL file!\n")
  quit(status = 1)
}
re$id <- phe$id[re$rows]

## Calculate the corrected phenotype
y_adj <- left_join(ebv_all, re, by = "id")
y_adj$Yhat <- y_adj$ebv + y_adj$re

## Remove missing values
adj <- subset(y_adj, !is.na(Yhat), select = c("id", "Yhat"))
nadj <- nrow(adj)

## Whether to append the results to an existing file
if (opt$append == "true") {
  append <- TRUE
} else {
  append <- FALSE
}

if (nadj > 0) {
  ## Output correction phenotype
  fwrite(adj, opt$out, col.names = FALSE, append = append, sep = " ")
  cat("phenotypes (", nrow(adj), ") corrected for all other effects output to:", opt$out, "\n")
} else {
  cat("error! no records in results file.\n")
}
