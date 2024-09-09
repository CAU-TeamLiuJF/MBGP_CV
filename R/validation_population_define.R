#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript

########################################################################################################################
## Version: 1.3.0
## Author:    Liweining liwn@jaas.ac.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-08-20
##
## Function：
## According to the requirements, divide the reference and validation used for calculating the prediction accuracy
##
##
## Usage: ./validation_population_define.R --phef "/path/to/phenotype/file" ...(Please refer to --help for detailed parameters)
##
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


# Load packages
cat("Loading packages needed...\n")
suppressPackageStartupMessages(library("getopt"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("reshape"))
suppressPackageStartupMessages(library("pedigree"))

## command parameters
spec <- matrix(
  c(
    "phef",      "I", 1, "character", "[Required] phenotype file name path",
    "pheCol",    "P", 1, "character", "[Required] phenotype columns",
    "seed",      "S", 1, "integer",   "[Optional] random number seed",
    "fold",      "F", 1, "integer",   "[Optional] fold",
    "rep",       "r", 1, "integer",   "[Optional] repeat times",
    "iyse",      "C", 1, "character", "[Optional] id[:year:year_star:year_end] [1]",
    "nonmiss",   "m", 1, "character", "[Optional] There must be no missing values in these columns, e.g. '1 2 3'",
    "fam",       "A", 1, "character", "[Optional] plink fam file name",
    "label",     "L", 1, "character", "[Optional] Label used to extract a specific rows",
    "labCol",    "B", 1, "character", "[Optional] Label columns",
    "gen",       "G", 1, "double",    "[Optional] last G geneorations as validations",
    "pedf",      "E", 1, "character", "[Optional] pedigree use for determine generations",
    "year",      "Y", 1, "double",    "[Optional] Minimum year of birth of candidate group individuals",
    "miss",      "M", 1, "double",    "[Optional] Identifier of the missing phenotype",
    "header",    "H", 1, "logical",   "[Optional] whether include header in output phenotype files/NULL",
    "outphe",    "O", 1, "character", "[Optional] Output phenotype file names",
    "outvid",    "o", 1, "character", "[Optional] Output validation id names",
    "missinval", "v", 1, "character", "[Optional] cols need to be set as missing in the validations, e.g. '2 3'",
    "keepid",    "k", 1, "character", "[Optional] fid iid file containing reference and vaildation",
    "outdir",    "D", 1, "character", "[Optional] Directory of output files",
    "rminvail",  "R", 0, "logical",   "[Optional] Remove individuals in phef that cannot be used as reference",
    "help",      "h", 0, "logical",   "[Optional] This is Help!"
  ),
  byrow = TRUE, ncol = 5
)
opt <- getopt::getopt(spec = spec)


## check parameters
if (!is.null(opt$help) || is.null(opt$phef) || is.null(opt$pheCol)) {
  ## print help message
  cat(paste(getopt::getopt(spec = spec, usage = TRUE), "\n"))
  quit(status = -1)
}

## Random number seed
if (is.null(opt$seed)) {
  opt$seed <- Sys.time()
  write.table(opt$seed, paste0(opt$outdir, "/pheno_group.seed"), col.names = FALSE, row.names = FALSE, quote = FALSE)
} else {
  set.seed(opt$seed)
}

## default parameters
if (is.null(opt$outphe)) opt$outphe <- "pheno.txt"
if (is.null(opt$miss)) opt$miss <- -99
if (is.null(opt$outvid)) opt$outvid <- "val.id"
if (is.null(opt$fold)) opt$fold <- 1
if (is.null(opt$rep)) opt$rep <- 1
if (is.null(opt$outdir)) opt$outdir <- getwd()
if (is.null(opt$iyse)) opt$iyse <- "1"
if (is.null(opt$header)) header <- FALSE else header <- TRUE

##  parse parameter
cols <- unlist(strsplit(opt$iyse, ":"))
id_col <- as.numeric(cols[1])
year_col <- as.numeric(cols[2])
ystar_col <- as.numeric(cols[3])
yend_col <- as.numeric(cols[4])

## Columns in the phenotype file
phe_cols <- as.numeric(unlist(strsplit(opt$pheCol, ",")))

## Phenotype file
pheno <- fread(opt$phef)
names(pheno)[id_col] <- "numid"
phe_ncol <- ncol(pheno)

## Genotype file
if (!is.null(opt$fam)) {
  fam <- read.table(opt$fam)
  geno_num <- fam[, 2]
  if (all(!geno_num %in% pheno$numid)) {
    cat("ID in the fam file is inconsistent with that in the phenotype file!\n")
    quit(status = -1)
  }
} else {
  geno_num <- unlist(subset(pheno, select = id_col))
}

## Year
if (!is.null(opt$year)) {
  year <- as.numeric(substr(pheno[, year_col], ystar_col, yend_col))
  year_index <- year >= opt$year
  cat("There are ", sum(year_index), " individuals with year >= ", opt$year, "\n")
} else {
  year_index <- rep(TRUE, nrow(pheno))
}

## Columns for fixed/random effects (integer columns) with missing values cannot be used as candidates,
## as residuals cannot be obtained
if (!is.null(opt$nonmiss)) {
  nonmiss_cols <- as.integer(unlist(strsplit(opt$nonmiss, " ")))
  effect_miss <- apply(subset(pheno, select = nonmiss_cols), 1, function(x) any(x == 0))
  if (any(effect_miss)) {
    cat("There are", sum(effect_miss), "individuals with missing values in these columns: ", opt$nonmiss, "\n")
  }
} else {
  effect_miss <- rep(FALSE, nrow(pheno))
}

## Individuals with missing phenotype values cannot be used as candidates, as residuals cannot be obtained
if (!is.null(opt$nonmiss)) {
  phe_miss <- apply(subset(pheno, select = phe_cols), 1, function(x) any(x == opt$miss))
  if (any(phe_miss)) {
    cat("There are", sum(phe_miss), "individuals with missing values in these columns: ", opt$pheCol, "\n")
  }
} else {
  phe_miss <- rep(FALSE, nrow(pheno))
}

## Generation
if (!is.null(opt$gen)) {
  if (is.null(opt$pedf)) {
    cat("Please provide a pedigree file to obtain individual generations.\n")
    quit(status = -1)
  } else {
    ped <- fread(opt$pedf)
    names(ped)[1:3] <- c("numid", "SIRE", "DAM")

    ## Pedigree sorting
    peds <- ped[order(orderPed(ped[, 1:3])), ]

    ## Calculate generations
    peds$Gen <- countGen(peds)

    ## Generations of individuals with phenotype data
    phe_gen <- pheno[, "numid"]
    phe_gen <- merge(phe_gen, peds, by = "numid", all.x = TRUE, sort = FALSE)
    phe_gen[is.na(phe_gen)] <- 0
    gens_phe <- unique(phe_gen$Gen)

    ## Report generations
    cat("Generation of phenotypic individuals:\n")
    print(summary(as.factor(phe_gen$Gen)))

    ## Selected generations
    gens_phe <- gens_phe[order(gens_phe)]
    gen_sel <- gens_phe[(length(gens_phe) - opt$gen + 1):length(gens_phe)]
    gen_index <- phe_gen$Gen %in% gen_sel
  }
} else {
  gen_index <- rep(TRUE, nrow(pheno))
}

## Sampling candidate population IDs
allid_index <- pheno$numid %in% geno_num & year_index & gen_index & !effect_miss & !phe_miss
allid <- unique(pheno$numid[allid_index])
nind <- length(allid)
nind_val <- nind %/% opt$fold

## phenotype containing qualified reference and candidate, and output fid-id for extracting genotype of qualified individuals
if (!is.null(opt$rminvail)) {
  ## Remove individuals who cannot be used as reference due to unmet data requirements
  ind_remove <- effect_miss | phe_miss | (!pheno$numid %in% geno_num)
  if (any(ind_remove)) {
    pheno <- subset(pheno, allid_index)
    write.table(pheno, opt$phef, row.names = FALSE, col.names = FALSE, quote = FALSE)
    cat("A total of", sum(!allid_index), "individuals that cannot be used as references are deleted.\n")
    cat("Output phenotype file containing reference and validation to:", opt$phef, "\n")

    if (!is.null(opt$fam) && !is.null(opt$keepid)) {
      keep_index <- fam$V2 %in% pheno$numid
      if (any(!keep_index)) {
        write.table(subset(fam, keep_index, select = 1:2), opt$keepid,
          row.names = FALSE, col.names = FALSE, quote = FALSE
        )
        cat("fid iid file containing reference and validation to:", opt$keepid, "\n")
      }
    }
  }
}

## Start sampling
for (r in 1:opt$rep) { # nolint
  ## IDs of all populations
  leftid <- allid

  for (f in 1:opt$fold) {
    ## Replace rep and fold in path parameters
    outdir <- gsub("#rep#", r, opt$outdir)
    outdir <- gsub("#val#", f, outdir)
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

    ## Replace rep and fold in filenames
    outphe <- gsub("#rep#", r, opt$outphe)
    outphe <- gsub("#val#", f, outphe)
    outvid <- gsub("#rep#", r, opt$outvid)
    outvid <- gsub("#val#", f, outvid)

    ## Randomly divide the candidate population into subsets
    if (f == opt$fold) {
      ref_id <- leftid
    } else {
      ref_id <- sample(leftid, nind_val, replace = FALSE)
      leftid <- leftid[!(leftid %in% ref_id)]
    }

    ## Generate phenotype file for reference population
    if (!is.null(opt$miss)) {
      ## Set the phenotype of the validation population to missing
      ref_pop <- pheno
      ref_index <- ref_pop$numid %in% ref_id
      ref_pop[ref_index, phe_cols] <- opt$miss

      ## Some effects/covariates need to be set to missing, as levels are not available when obtaining genetic information
      if (!is.null(opt$missinval)) {
        miss_val <- as.integer(unlist(strsplit(opt$missinval, " ")))
        int_index <- apply(ref_pop, 2, function(x) all(floor(x) == x))
        num_int <- which(!int_index)[1] - 1
        for (i in miss_val) {
          if (i > num_int) {
            ref_pop[ref_index, i] <- opt$miss
          } else {
            ref_pop[ref_index, i] <- 0 ## Missing values for integer variables in DMU are indicated by 0
          }
        }
      }
    } else {
      ## Remove records of the validation population from the phenotype file
      ref_pop <- pheno[!pheno$numid %in% ref_id, ]
    }

    ## Check for errors
    if (nrow(ref_pop) < 1 || length(ref_id) < 1) {
      cat("error! rep", r, "fold", f, ":", length(ref_id), "/", nrow(ref_pop), "\n")
      quit(status = 1)
    }

    ## Log
    if (r == 1) {
      cat("rep", r, "fold", f, ":", length(ref_id), "/", nrow(ref_pop), "\n")
    }

    ## Write out the reference population phenotype with validation population individuals removed (or set to missing)
    write.table(ref_pop, paste0(outdir, "/", outphe), row.names = FALSE, col.names = header, quote = FALSE)

    ## Write out validation population IDs
    if (!is.null(opt$fam)) {
      val_fam <- fam[fam$V2 %in% ref_id, 1:2]
    } else {
      val_fam <- data.frame(fid = ref_id, iid = ref_id)
    }
    write.table(val_fam, paste0(outdir, "/", outvid), row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

cat("phenotype grouping completed.\n")
