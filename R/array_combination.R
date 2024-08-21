#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript

########################################################################################################################
## Version: 1.3.0
## Author:    Liweining liwn@cau.edu.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-08-20
##
## Function：
## Output provides a random combination of strings
##
##
## Usage: ./array_combination.R --array "breedA breedB breedC" ...(Please refer to --help for detailed parameters)
##
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


## Load packages
suppressPackageStartupMessages(library("getopt"))

## Command-line parameters
spec <- matrix(
  c("array",  "a", 2, "character", "[Required] character array separated by space",
    "label",  "l", 2, "character", "[Optional] combination must include/exclude this string",
    "type",   "t", 2, "character", "[Optional] include or exclude [include]",
    "min",    "m", 2, "integer",   "[Optional] minimum number of elements contained in a combination",
    "max",    "M", 2, "integer",   "[Optional] maximum number of elements contained in a combination",
    "out",    "o", 2, "character", "[Optional] output file name",
    "append", "A", 2, "logical",   "[Optional] append to output file",
    "help",   "h", 0, "logical",   "This is Help!"),
  byrow = TRUE, ncol = 5)
opt <- getopt::getopt(spec = spec)

## Check parameters
if (!is.null(opt$help) || is.null(opt$array)) {
  cat(paste(getopt::getopt(spec = spec, usage = TRUE), "\n"))
  quit()
}

## Default parameters
if (is.null(opt$min)) opt$min <- 2
# if (is.null(opt$out)) opt$out <- "combination.txt"
if (is.null(opt$type)) opt$type <- "include"
if (!is.null(opt$append)) append <- TRUE else append <- FALSE

## Splitting a string
array <- strsplit(opt$array, " ")[[1]]
np <- length(array)
if (is.null(opt$max)) opt$max <- np

## Function that lists all possible combinations of a specified number of elements
com <- function(x, m) {
  if (m < 2) {
    return(x)
  }

  ## Free combination
  df <- combn(x, m)

  ## Paste each line of variables into a string
  result <- apply(df, 2, paste, collapse = " ")

  return(result)
}

## combination
all_com <- c()
for (i in opt$min:opt$max) {
  com_i <- com(array, i)
  all_com <- c(com_i, all_com)
}

## Only keep combinations with target labels
if (!is.null(opt$label)) {
  ## Determine whether the string is in the vector
  if (!opt$label %in% array) {
    cat("label not in array!\n")
    quit(1)
  }

  ## Find combinations containing a certain tag
  idx <- grepl(opt$label, all_com)

  ## Extract
  if (opt$type == "include") {
    all_com <- all_com[idx]
  } else {
    all_com <- all_com[!idx]
  }
}

## Output the result to file
n_com <- length(all_com)
cat("The number of all possible combinations is:", n_com, "\n")
if (!is.null(opt$out)) {
  write.table(all_com, opt$out, quote = FALSE, col.names = FALSE, row.names = FALSE, append = append)
} else {
  print(all_com)
}
