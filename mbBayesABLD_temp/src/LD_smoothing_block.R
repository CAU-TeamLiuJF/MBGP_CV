#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript

########################################################################################################################
## Version: 1.3.0
## Author:    Liweining liwn@jaas.ac.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-08-20
##
## Function：
## Genomic block partitioning based on the mean R2 value calculated by a self written C language program
##
##
## Usage: ./cubic_smoothing_block.R --r2 "/path/to/r2/file" ...(Please refer to --help for detailed parameters)
##
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


## Loading required packages
cat("Loading required packages... \n\n")
suppressPackageStartupMessages(library("getopt"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("Cairo"))

## Parameter list
spec <- matrix(
  c("r2",   "I", 1, "character", "[Required] R2 mean\n",
    "spar", "S", 1, "double",    "[Optional] Smoothness (0-1) when fitting curves [0.2]\n",
    "diff", "D", 1, "double",    "[Optional] [0.05]\n",
    "bim",  "B", 1, "character", "[Optional] bim file name [NULL]\n",
    "out",  "O", 1, "character", "[Optional] Output file name [cubic_LD_block.txt]\n",
    "plot", "P", 0, "logical",   "[Optional] Output plot [FALSE]\n",
    "help", "h", 0, "logical",  "This is Help!"),
  byrow = TRUE, ncol = 5
)
opt <- getopt(spec = spec)

## Check parameters
if (!is.null(opt$help) || is.null(opt$r2)) {
  cat(paste(getopt(spec = spec, usage = TRUE), "\n"))
  quit()
}

## Default parameters
if (is.null(opt$out)) opt$out <- "cubic_LD_block.txt"
if (is.null(opt$spar)) opt$spar <- 0.2
if (is.null(opt$diff)) opt$diff <- 0.05

## Function definition
find_local_min <- function(x, diff = 0.2) { # nolint
  n <- length(x)
  peaks <- which(diff(sign(diff(x))) < 0) + 1
  valleys <- which(diff(sign(diff(x))) > 0) + 1

  ## Add endpoints
  if (x[1] > x[2]) {
    peaks <- c(1, peaks)
  } else if (x[1] > x[2]) {
    valleys <- c(1, valleys)
  }
  if (x[n] < x[n - 1]) {
    valleys <- c(valleys, n)
  } else if (x[n] > x[n - 1]) {
    peaks <- c(peaks, n)
  }

  final <- c()
  for (i in valleys) {
    left_peak <- TRUE
    if (any(peaks < i)) {
      left_peak <- max(peaks[peaks < i])
      diff_rate <- abs(x[i] - x[left_peak]) / x[i]
      if (diff_rate < diff) left_peak <- FALSE
    }

    right_peak <- TRUE
    if (any(peaks > i)) {
      right_peak <- min(peaks[peaks > i])
      diff_rate <- abs(x[i] - x[right_peak]) / x[i]
      if (diff_rate < diff) right_peak <- FALSE
    }

    if (left_peak && right_peak) final <- c(final, i)
  }

  ## Exclude endpoints
  final <- final[final != 1 & final != n]

  return(final)
}

## Read file
r2 <- fread(opt$r2)

## Fit curve
fit <- smooth.spline(r2$V1, spar = 0.2)

## Search for local minimum points that meet the criteria
point <- find_local_min(fit$y, 0.05)

## Calculate the number of SNPs in each region
point2 <- c(0, point, length(fit$y))
n <- length(point2)
nsnp <- point2[2:n] - point2[1:(n - 1)]
nsnp[length(nsnp)] <- nsnp[length(nsnp)] + 1

## SNP count check
if (sum(nsnp) != (nrow(r2) + 1)) {
  cat("The number of SNPs in the result is not equal to that in the input file!\n")
  quit()
}

## Output SNP count file with additional information
if (!is.null(opt$bim)) {
  if (file.exists(opt$bim)) {
    bim <- fread(opt$bim)
    if (nrow(bim) == sum(nsnp)) {
      # CHR START STOP nSNP
      point2[1] <- 1
      nsnp <- data.frame(CHR = bim$V1[1], START = bim$V4[point2[1:(n - 1)]], STOP = bim$V4[point2[2:n]], nSNP = nsnp)
    } else {
      cat("The number of SNPs in the result is not equal to that in the bim file!\n")
      quit()
    }
  } else {
    cat(opt$bim, "not found!\n")
    quit()
  }
}

## Output SNP count file
write.table(nsnp, opt$out, row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("Bins definition file has been output to:", opt$out, "\n")

## Plotting
if (!is.null(opt$plot)) {
  CairoPNG(paste0(opt$out, ".png"), width = 1600, height = 800)
  plot(r2$V1, col = "grey", pch = 20,
    ylab = "local mean r2",
    xlab = "SNP rank")
  lines(fit, col = "red", lwd = 2)
  for (h in point) {
    abline(v = h, col = "blue", lwd = 1)
  }
  dev.off()
  cat("Plot has been output to:", paste0(opt$out, ".png"), "\n")
}
