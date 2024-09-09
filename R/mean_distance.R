#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript

########################################################################################################################
## Version: 1.3.0
## Author:    Liweining liwn@jaas.ac.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-08-20
##
## Function：
## Calculate the mean distance between breeds
##
##
## Usage: ./mean_distance.R --prefix "/plink/binary/file/prefix" ...(Please refer to --help for detailed parameters)
##
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


### Load packages required
cat("Loading required packages... \n\n")
suppressPackageStartupMessages(library("getopt"))
suppressPackageStartupMessages(library("data.table"))

## parameters list
spec <- matrix(
  c(
    "prefix", "f", 1, "character", "[Required] mdist file path\n",
    "indexf", "i", 1, "character", "[Optional] family id\n",
    "out",    "o", 1, "character", "[Optional] output file name prefix [mdist_summary]\n",
    "help",   "h", 0, "logical", "This is Help!"
  ),
  byrow = TRUE, ncol = 5
)
opt <- getopt(spec = spec)

## Check parameters
if (!is.null(opt$help) || is.null(opt$prefix)) {
  cat(paste(getopt(spec = spec, usage = TRUE), "\n"))
  quit()
}

## Default parameters
if (is.null(opt$out)) opt$out <- "mdist_summary"

## File name
dist_file <- paste0(opt$prefix, ".mdist")
if (is.null(opt$indexf)) {
  id_file <- paste0(opt$prefix, ".mdist.id")
} else {
  id_file <- opt$indexf
}

## Distance matrix
if (file.exists(dist_file)) {
  mdist <- fread(dist_file)
  mdist <- as.matrix(mdist)
} else {
  cat(dist_file, "not found!\n")
  quit()
}

## family ID
if (file.exists(id_file)) {
  id <- fread(id_file)
  names(id)[1] <- "fid"
} else {
  cat(id_file, "not found!\n")
  quit()
}

## Check the rationality of the breed ID
fids <- unique(id$fid)
if (nrow(id) != nrow(mdist)) {
  cat("The number of id is not equal to the number of rows in mdist file!\n")
  quit()
} else if (length(fids) == nrow(mdist)) {
  cat("The number of family IDS is equal to the individual IDS, please check\n")
  quit()
} else if (length(fids) == 1) {
  cat("Only one family id:", fids, "\n")
  quit()
}

## Calculate the mean
mean <- matrix(NA, nrow = length(fids), ncol = length(fids))
for (i in seq_len(length(fids))) {
  for (j in 1:i) {
    if (i == j) {
      mean[i, j] <- NA
    } else {
      fid_i <- id$fid == fids[i]
      fid_j <- id$fid == fids[j]
      mean[i, j] <- mean(mdist[fid_i, fid_j])
    }
  }
}

## rename
colnames(mean) <- fids
rownames(mean) <- fids

## Output
output_name <- paste0(opt$out, ".csv")
write.csv(mean, output_name, row.names = TRUE, na = "")

cat("average Genetic distances output to:", output_name, "\n")
