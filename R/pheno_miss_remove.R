#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript

########################################################################################################################
## Version: 1.3.0
## Author:    Liweining liwn@jaas.ac.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-08-20
##
## Function：
## Delete individuals with missing phenotypes
##
##
## Usage: ./pheno_miss_remove.R --file "/path/to/file" ...(Please refer to --help for detailed parameters)
##
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


## Load package required
suppressPackageStartupMessages(library("getopt"))

## Command-line parameters
spec <- matrix(
  c("file",     "I",  1, "character", "[Required] input file name",
    "col",      "C",  1, "integer",   "[Required] column to find missing value",
    "idC",      "i",  1, "integer",   "[Required] id column/1",
    "miss",     "M",  1, "double",    "[Optional] value to be treat as missing value/-99",
    "backName", "B",  1, "character", "[Optional] new file names to save origin file/*.bc",
    "map",      "m",  1, "character", "[Optional] plink map/fam file",
    "missid",   "s",  1, "character", "[Optional] file contain missing phenotype id/No_phe.id",
    "out",      "o",  1, "character", "[Optional] file contain missing phenotype id/No_phe.id",
    "NotRm",    "n",  0, "logical",   "[Optional] do not remove id with missing phenotype",
    "help",     "h",  0, "logical",   "This is Help!"),
  byrow = TRUE, ncol = 5)
opt <- getopt(spec = spec)

## Check parameters
if (!is.null(opt$help) || is.null(opt$file) || is.null(opt$col)) {
  cat(paste(getopt(spec = spec, usage = TRUE), "\n"))
  quit()
}

## Default parameters
if (is.null(opt$idC)) opt$idC <- 1
if (is.null(opt$miss)) opt$miss <- -99
if (is.null(opt$backName)) opt$backName <- paste0(opt$file, ".bc")
if (is.null(opt$missid)) opt$missid <- "No_phe.id"

## Read phenotype
phe <- read.table(opt$file)

## Rename
if (opt$col > ncol(phe) || opt$idC > ncol(phe)) {
  cat("parameter --col (", opt$col, ") or --idC (", opt$idC, ") bigger than phenotype file columns!\n")
  quit()
} else {
  names(phe)[opt$col] <- "value"
  names(phe)[opt$idC] <- "id"
}

## Search for missing phenotypic individuals
missid <- phe$id[phe$value <= opt$miss]
if (length(missid) > 0) {
  if (!is.null(opt$map)) {
    ## Output genotype individuals with missing phenotypes
    map <- read.table(opt$map)
    fid_iid <- map[map$V2 %in% missid, 1:2]
    write.table(fid_iid, opt$missid, row.names = FALSE, col.names = FALSE, quote = FALSE)
  } else {
    write.table(missid, opt$missid, row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  cat("individuals ID with the missing phenotype output to:", opt$missid, "\n")

  ## Delete individuals with missing phenotypes in phenotype files
  if (is.null(opt$NotRm)) {
    phe_new <- phe[!phe$id %in% missid, ]
    write.table(phe_new, opt$out, row.names = FALSE, col.names = FALSE, quote = FALSE)
    # write.table(phe, opt$backName, row.names = FALSE, col.names = FALSE, quote = FALSE)
    cat(
      "phenotype (", nrow(phe_new), ") that no missing values in column", opt$col,
      "has been output to", opt$out, "\n"
    )
    cat("remove", length(missid), "individuals with the missing value in input files\n")
  }
} else {
  write.table(phe, opt$out, row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat("No individuals with the missing phenotype found, output file copy.\n")
}
