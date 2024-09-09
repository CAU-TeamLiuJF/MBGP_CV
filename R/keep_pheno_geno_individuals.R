#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript

########################################################################################################################
## Version: 1.3.0
## Author:    Liweining liwn@jaas.ac.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-08-20
##
## Function：
## Used to screen individuals with genotype information that simultaneously possess phenotype information
##
##
## Usage: ./keep_pheno_geno_individuals.R --famf "/path/to/plink/famf/file" ...(Please refer to --help for detailed parameters)
##
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


## Load the required packages
suppressPackageStartupMessages(library("getopt"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dplyr"))

## Command-line arguments
spec <- matrix(
  c("famf",        "G", 1, "character", "[Required] plink fam file",
    "phef",        "P", 1, "character", "[Required] phenotype file",
    "phec",        "p", 1, "integer",   "[Optional] phenotype column in phenotype file",
    "phe_idc",     "I", 1, "integer",   "[Optional] id column in phenotype file",
    "miss",        "m", 1, "double",    "[Optional] miss value for real value",
    "num_int",     "i", 1, "integer",   "[Optional] number of integer variables",
    "mean",        "M", 1, "integer",   "[Optional] population mean columns",
    "breed",       "b", 1, "integer",   "[Optional] breed effect column",
    "out",         "O", 1, "character", "[Optional] output file name prefix/phef\n",
    "rm",          "r", 1, "character", "[Optional] remove all/single trait missing recodes/NULL\n",
    "rmOut",       "R", 1, "character", "[Optional] remove ids\n",
    "header",      "H", 0, "logical",   "[Optional] add header in output phenotype file",
    "sort",        "d", 0, "logical",   "[Optional] sort individuals in phenotype file (consistent with famf)",
    "help",        "h", 0, "logical",   "This is Help!"),
  byrow = TRUE, ncol = 5)
opt <- getopt(spec = spec)

## Check parameters
if (!is.null(opt$help) || is.null(opt$phef) || is.null(opt$famf)) {
  cat(paste(getopt(spec = spec, usage = TRUE), "\n"))
  quit()
}

## Default parameters
if (is.null(opt$header)) header <- FALSE else header <- TRUE
if (is.null(opt$miss)) opt$miss <- -99
if (is.null(opt$phe_idc)) opt$phe_idc <- 1
if (is.null(opt$rmOut)) opt$rmOut <- "miss_phe.id"
if (is.null(opt$out)) opt$out <- opt$phef

## Phenotype file
if (file.exists(opt$phef)) {
  phe <- fread(opt$phef)
} else {
  cat(opt$phef, "not found!\n")
  quit()
}

## Number of integer columns
if (is.null(opt$num_int)) {
  index <- apply(phe[1, ], 1, function(x) as.integer(x) == as.numeric(x))
  opt$num_int <- which(!index)[1] - 1
}

## Add a helper mean column
if (ncol(phe) < 2) phe$mean <- 1

## Naming
names(phe)[opt$phe_idc] <- "id"
if (!is.null(opt$mean)) names(phe)[opt$mean] <- "mean"
if (!is.null(opt$breed)) names(phe)[opt$breed] <- "breed"
names(phe)[(opt$num_int + 1):ncol(phe)] <- paste0("R", 1:(ncol(phe) - opt$num_int))
phe_names <- names(phe)

## Genotyped individual IDs
fam <- fread(opt$famf)
names(fam) <- c("fid", "id", "sire", "dam", "sex", "phe")
cat("genotyped individuals in plink:", nrow(fam), "\n")

if (all(fam$id %in% phe$id) && is.null(opt$rm)) {
  ## All genotyped individuals are present in the phenotype file
  phe <- phe[phe$id %in% fam$id, ]

  if (!is.null(opt$out)) {
    cat("Phenotypes of all genotyped individuals (", nrow(fam), ") output to:", opt$out, "\n")
    write.table(phe, opt$out, row.names = FALSE, col.names = header, quote = FALSE)
  }
  quit(status = 1)
} else if (all(!fam$id %in% phe$id)) {
  ## None of the genotyped individuals are present in the phenotype file
  cat("ID in genotype file and phenotype file are all inconsistent.\n")
  quit(status = 1)
}

## Merge files
phe_g <- left_join(fam, phe, by = "id")

## Fill in population mean
if (!is.null(opt$mean)) {
  phe_g$mean <- 1
}

## Breed fixed effect variable
if (!is.null(opt$breed)) {
  phe_breed <- unique(phe_g$breed)
  if (any(is.na(phe_breed))) {
    if (length(unique(phe_g$fid)) > 10) {
      cat("Please provide the breed label in the first column of the fam file!\n")
      print(head(fam, n = 3))
      quit(status = 1)
    } else {
      phe_g$breed <- as.numeric(as.factor(phe_g$fid))
    }
  }
}

## Select original phenotype columns
phe_out <- subset(phe_g, select = phe_names)

## Fill missing values in other columns
for (i in seq_len(ncol(phe_out))) {
  ## Integer missing values are set to 0, real numbers are set to the specified value
  if (i <= opt$num_int) {
    miss <- 0
  } else {
    miss <- opt$miss
  }

  ## Fill missing values
  na_index <- is.na(phe_out[[i]])
  phe_out[[i]][na_index] <- miss
}

## Remove records of individuals missing all phenotypes
if (!is.null(opt$rm)) {
  if (opt$rm == "all") {
    all_miss_id <- apply(phe_out[, (opt$num_int + 1):ncol(phe_out)], 1, function(x) all(x == opt$miss))
    all_miss_id <- phe_out$id[all_miss_id]
    phe_out <- subset(phe_out, !id %in% all_miss_id)
    mv_id <- fam[fam$id %in% all_miss_id, 1:2]
    if (nrow(mv_id) > 0) {
      cat(nrow(mv_id), "individual IDs missing all phenotypes output to:", opt$rmOut, "\n")
      write.table(mv_id, opt$rmOut, row.names = FALSE, col.names = header, quote = FALSE)
    }
  } else if (opt$rm == "single") {
    if (is.null(opt$phec)) opt$phec <- opt$num_int + 1
    names(phe_out)[opt$phec] <- "phe"
    miss_index <- phe_out$phe == opt$miss

    if (sum(miss_index) > 0) {
      ## Individuals missing specified phenotypes
      mv_id <- fam[fam$id %in% phe_out$id[miss_index], c("fid", "id")]

      ## Filter out individuals with phenotypes
      phe_out <- phe_out[!miss_index, ]

      cat(nrow(mv_id), "individuals missing specified phenotypes output to:", opt$rmOut, "\n")
      write.table(mv_id, opt$rmOut, row.names = FALSE, col.names = header, quote = FALSE)
    }
  }
}

## Ensure individuals in the phenotype file are ordered consistently with the genotype file famf
if (!is.null(opt$sort)) {
  phe_out <- phe_out[phe_out$id %in% fam$id, ]
}

## Output file
if (!is.null(opt$out)) {
  write.table(phe_out, opt$out, row.names = FALSE, col.names = header, quote = FALSE)
  cat("Phenotypes of all genotyped individuals (", nrow(phe_out), ") output to: ", opt$out, "\n", sep = "")
}
