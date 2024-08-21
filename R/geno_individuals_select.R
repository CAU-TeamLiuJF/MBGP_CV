#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript

########################################################################################################################
## Version: 1.3.0
## Author:    Liweining liwn@cau.edu.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-08-20
##
## Function：
## Select genotype individuals based on the phenotype file_data_ * * *. txt output by QMSim
##
##
## Usage: ./geno_individuals_select.R --dir_all "/path/to/phenotype/file" ...(Please refer to --help for detailed parameters)
##
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


## Select genotyped individuals based on the phenotype file _data_***.txt output from QMSim

# Load packages
suppressPackageStartupMessages(library("getopt"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dplyr"))

## Command parameters
spec <- matrix(
  c("dataf",    "I", 1, "character", "[Required] phenotype file name/path",
    "nlitter",  "l", 1, "double",    "[Optional] numbers individuals select/litter [2]",
    "nsel",     "N", 1, "double",    "[Optional] total selected individuals number [2]",
    "gen_all",  "g", 1, "character", "[Optional] generations in genotype file. 8-10 [last one]",
    "gen_sel",  "G", 1, "character", "[Optional] generation to select eg. 10 [last one]",
    "fid",      "F", 1, "character", "[Optional] Output ID file with FID column",
    "out",      "O", 1, "character", "[Optional] Output file name prefix [keep_geno_ids.txt]\n",
    "seed",     "s", 2, "integer",   "[Optional] [NULL]",
    "outIndex", "i", 0, "logical",   "[Optional] selected Ind index in dataf \n",    
    "help",     "h", 0, "logical",   "This is Help!"),
  byrow = TRUE, ncol = 5)
opt <- getopt::getopt(spec = spec)

## Check parameters
if (!is.null(opt$help) || is.null(opt$dataf)) {
  ## Print help message
  cat(paste(getopt::getopt(spec = spec, usage = TRUE), "\n"))
  quit()
}

## Default parameters
if (is.null(opt$out)) opt$out <- "keep_geno_ids.txt"
if (is.null(opt$nlitter)) opt$nlitter <- 2

## Random seed
if (is.null(opt$seed)) {
  opt$seed <- Sys.time()
  write.table(opt$seed, paste0(out, ".seed"))
}
set.seed(opt$seed)

## Read phenotype file
if (file.exists(opt$dataf)) {
  data <- fread(opt$dataf)
} else {
  cat(opt$dataf, "not found.\n")
  quit()
}

## Generations of genotyped individuals output from QMSim
if (!is.null(opt$gen_all)) {
  gen <- as.numeric(unlist(strsplit(opt$gen_all, "-")))
  if (length(gen) > 1) gen <- gen[1]:gen[2]
  data <- subset(data, G %in% gen)
}

## Selected generation
if (!is.null(opt$gen_sel)) {
  gen <- as.numeric(unlist(strsplit(opt$gen_sel, "-")))
  if (length(gen) > 1) gen <- gen[1]:gen[2]
  data_gen <- subset(data, G %in% gen)
  cat("Selected generation:", gen, "\n")
} else {
  gen <- unique(data$G)
  gen <- gen[length(gen)]
  data_gen <- subset(data, G %in% gen)
}
cat("Number of candidates:", nrow(data_gen), "\n")
cat("Number of litters:", length(unique(data_gen$Dam)), "\n")

## Select nlitter individuals per litter
sample_litter <- data_gen %>% group_by(G, Dam) %>% slice_sample(n = opt$nlitter)

## Select only a portion of litters
if (!is.null(opt$nsel)) {
  ## Choose the appropriate number of litters
  nlitter_per_gen <- opt$nsel %/% opt$nlitter %/% length(gen)
  dam_sel <- sample_litter %>% select(G, Dam) %>% unique %>% group_by(G) %>% slice_sample(n = nlitter_per_gen)

  ## Filter individuals from selected litters
  sample_litter <- subset(sample_litter, Dam %in% dam_sel$Dam)
}

if (is.null(opt$outIndex)) {
  ## Whether the output file contains FID (for plink to extract individual genotypes)
  if (!is.null(opt$fid)) {
    keep <- data.frame(iid = opt$fid, iid = sample_litter$Progeny)
  } else {
    keep <- sample_litter$Progeny
  }
  cat("Output selected individuals ID\n")
  nsel_final <- nrow(keep)
} else {
  ## Output the position of selected individuals in the original file (indicated by a 0/1 variable column)
  keep <- as.integer(data$Progeny %in% sample_litter$Progeny)
  nsel_final <- sum(keep)
  cat("Output index of selected individuals in input file\n")
}

## Output file
write.table(keep, opt$out, row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("Total selected:", nsel_final, "\n")
cat("IDs selected output to:", opt$out, "\n")
