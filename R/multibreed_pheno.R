#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript


########################################################################################################################
## Version: 1.3.0
## Author:    Liweining liwn@cau.edu.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-08-20
##
## Function：
## Merge the DMU phenotype files of multiple breeds for joint prediction
##
##
## Usage: ./multibreed_pheno.R --phef "/path/to/phenotype/file" ...(Please refer to --help for detailed parameters)
##
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################

## command parameters
spec <- matrix(
  c("phef",    "P", 1, "character", "[Required] phenotype file name",
    "pops",    "N", 1, "character", "[Required] population names",
    "nInt",    "n", 1, "integer",   "[Optional] number of integer columns",
    "rep",     "r", 1, "integer",   "[Optional] repeat times [1]",
    "fold",    "f", 1, "integer",   "[Optional] cross validation fold [1]",
    "type",    "m", 1, "character", "[Optional] single/multi(two trait model) [single]",
    "method",  "e", 1, "character", "[Optional] Genomic prediction model, BLUP/GBLUP/ssGBLUP/bayesR/mbBayesAB [GBLUP]",
    "pheCol",  "p", 1, "integer",   "[Optional] Number of columns of phenotypic in real columns [1]",
    "miss",    "M", 1, "double",    "[Optional] Identifier of the missing phenotype [-99]",
    "fixPop",  "F", 2, "character", "[Optional] Add fixed breed effects to the model [NULL]",
    "addpop",  "a", 1, "character", "[Optional] Add breed labels on some effects, like:'2:3', e.g. 2-breed, 3-breed",
    "header",  "H", 1, "logical",   "[Optional] whether include header in output phenotype files/NULL",
    "out",     "o", 1, "character", "[Optional] Output file prefix\n",
    "overwri", "O", 2, "character", "[Optional] whether overwrite file if pheno.txt already exist\n",
    "help",    "h", 0, "logical",   "This is Help!"),
  byrow = TRUE, ncol = 5)
opt <- getopt::getopt(spec = spec)

## Parameters check
if (!is.null(opt$help) || is.null(opt$phef) || is.null(opt$pops)) {
  ## print help message
  cat(paste(getopt::getopt(spec = spec, usage = TRUE), "\n"))
  quit(status = -1)
}

# Load packages
# cat('Loading packages needed...\n')
suppressPackageStartupMessages(library("data.table"))

## Default parameters
if (is.null(opt$rep)) opt$rep <- 1
if (is.null(opt$fold)) opt$fold <- 1
if (is.null(opt$out)) opt$out <- "pheno.txt"
if (is.null(opt$miss)) opt$miss <- -99
if (is.null(opt$pheCol)) opt$pheCol <- 1
if (is.null(opt$type)) opt$type <- "single"
if (is.null(opt$method)) opt$method <- "GBLUP"
if (is.null(opt$header)) header <- FALSE else header <- TRUE

## Parse parameters
pops <- strsplit(opt$pops, " ")[[1]]
np <- length(pops)

## Merge the phenotypes
for (r in 1:opt$rep) { # nolint
  for (f in 1:opt$fold) {
    ## Output filename
    out <- gsub("#rep#", r, opt$out)
    out <- gsub("#val#", f, out)

    ## Ensure that the target folder exists
    outdir <- dirname(out)
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

    ## Check if the output file exists
    if (file.exists(out) && opt$overwri != "true") next

    ## Merge phenotype files of various breeds
    phe <- data.table()
    for (i in seq_len(np)) {
      nint <- opt$nInt
      ## File name preparation
      phefi <- gsub("#rep#", r, opt$phef)
      phefi <- gsub("#val#", f, phefi)
      phefi <- gsub("#breed#", pops[i], phefi)

      ## Check if phenotype files of within-breed exist
      if (!file.exists(phefi)) {
        cat(paste("Error: ", phefi, " not exist!\n"))
        quit(status = -1)
      }

      ## load phenotype
      phei <- fread(phefi)

      ## Count the integer variable column
      if (is.null(nint)) nint <- sum(sapply(phei, class) == "integer")

      ## Extract the required columns
      phei2 <- subset(phei, select = c(1:nint, nint + opt$pheCol))
      names(phei2) <- c(paste0("I", 1:nint), "P")

      # ## The missing phenotype is represented by NA in bayesR software
      # if (opt$method == "bayesR") {
      #   phei2$P[phei2$P <= opt$miss] <- "NA"
      # }

      ## add breed fixed effect
      if (!is.null(opt$fixPop)) {
        phei2$pop <- i
        phei2 <- subset(phei2, select = c(1:nint, ncol(phei2), ncol(phei2) - 1))
        nint <- nint + 1
      }

      ## merge
      if (opt$type == "multi") {
        ## multi-trait model
        phei2[, paste0("P", seq_len(np))] <- opt$miss
        phei2[[paste0("P", i)]] <- phei2[[nint + 1]]
        phei2 <- subset(phei2, select = names(phei2) != "P") # drop the original phenotype
      } else if (opt$type != "single") {
        cat("The type parameter is incorrect. It can only be single/multi\n")
        quit(status = -1)
      }
      phe <- rbind(phe, phei2)
    }

    if (!is.null(opt$addpop)) {
      ## Add breed identifiers to the designated fixed effects group for differentiation
      cols <- as.numeric(unlist(strsplit(opt$addpop, ":")))
      for (c in cols) {
        if (c > opt$nInt) {
          cat("Non-integer columns cannot add breed identifiers.\n")
          quit(status = -1)
        }
        phe[[paste0("V", c)]] <- paste0(phe$pop, phe[[paste0("V", c)]])
      }
    }

    ## output to file
    fwrite(phe, out, col.names = header, sep = " ")

    ## report
    if (r == 1 && opt$fold > 1) {
      ## Number of reference
      if (opt$type == "single") {
        nind_ref <- sum(subset(phe, select = ncol(phe)) != opt$miss)
      } else {
        nind_val <- sum(apply(phe[, (ncol(phe) - np + 1):ncol(phe)], 1, function(x) all(x == opt$miss)))
        nind_ref <- nrow(phe) - nind_val
      }

      cat("records in fold ", f, " reference populations: ", nind_ref, "\n", sep = "")
    }
  }
}
