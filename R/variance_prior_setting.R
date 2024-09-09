#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript


########################################################################################################################
## Version: 1.3.0
## Author:    Liweining liwn@jaas.ac.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-08-20
##
## Function：
## Merge the DMU phenotype files of multiple breeds for joint prediction
##
##
## Usage: ./variance_prior_setting.R --filea "/path/to/file" ...(Please refer to --help for detailed parameters)
##
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################

## Load required packages (install them in advance if not already installed)
suppressPackageStartupMessages(library("getopt"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("Matrix"))

## Get command-line parameters
spec <- matrix(
  c("filea",    "p", 1, "character", "[Required] phenotype file 1",
    "out",      "o", 1, "character", "[Required] output file name",
    "fileb",    "P", 1, "character", "[Optional] phenotype file 2",
    "pcol",     "c", 1, "double",    "[Optional] phenotype column in phenotype file [first real type column]",
    "type",     "t", 1, "character", "[Optional] Format of the output variance file. single/multi [multi]",
    "method",   "m", 1, "character", "[Optional] genomic prediction models. GBLUP/mbBayesAB [GBLUP]",
    "rep",      "r", 1, "integer",   "[Optional] repeat times [1]",
    "fold",     "f", 1, "integer",   "[Optional] cross validation fold [1]",
    "var",      "v", 1, "character", "[Optional] pheno/predict [pheno]",
    "h2",       "a", 1, "double",    "[Optional] heritability [0.5]",
    "h2B",      "A", 1, "double",    "[Optional] heritability of trait 2 [h2]",
    "rg",       "g", 1, "double",    "[Optional] genetic correlation [0.001]",
    "re",       "e", 1, "double",    "[Optional] residual correlation [0.001]",
    "add_rf",   "d", 1, "integer",   "[Optional] additive effect random group in dmu [1]",
    "rg_local", "l", 1, "character", "[Optional] file contains correlations of each bins [1]",
    "miss",     "n", 1, "double",    "[Optional] missing value [-99]",
    "overwri",  "O", 2, "character", "[Optional] overwrie the existing result file",
    "norec",    "R", 0, "logical",   "[Optional] non-existing covariances between residuals",
    "help",     "h", 0, "logical",   "This is Help!"),
  byrow = TRUE, ncol = 5)
opt <- getopt(spec = spec)

if (!is.null(opt$help) || is.null(opt$filea) || is.null(opt$out)) {
  cat(paste(getopt(spec = spec, usage = TRUE), "\n"))
  quit()
}

## Default parameters
if (is.null(opt$h2)) opt$h2 <- 0.5
if (is.null(opt$h2B)) opt$h2B <- opt$h2
if (is.null(opt$rg)) opt$rg <- 0.001
if (is.null(opt$re)) opt$re <- 0.001
if (is.null(opt$miss)) opt$miss <- -99
if (is.null(opt$type)) opt$type <- "multi"
if (is.null(opt$method)) opt$method <- "GBLUP"
if (is.null(opt$add_rf)) opt$add_rf <- 1
if (is.null(opt$var)) opt$var <- "pheno"
if (is.null(opt$fileb)) {
  opt$fileb <- opt$filea
  cat("warn: fileb not set! set to filea\n")
}

output <- FALSE
for (r in 1:opt$rep) { # nolint
  for (f in 1:opt$fold) {
    ## Replace placeholders in the path
    filea <- gsub("#val#", f, gsub("#rep#", r, opt$filea))
    fileb <- gsub("#val#", f, gsub("#rep#", r, opt$fileb))
    out <- gsub("#val#", f, gsub("#rep#", r, opt$out))

    ## Check if there are any phenotype files in the folder
    if (file.exists(out) && opt$overwri != "true") next

    output <- TRUE

    ## Check if the file exists
    if (!file.exists(filea) || !file.exists(fileb)) {
      cat("warn: file a or b not found, prior information of rep", r, "fold", f, "will not be generated!\n")
      next
    }

    ## Read the files
    data1 <- fread(filea)
    data2 <- fread(fileb)

    ## Obtain variance components
    if (opt$var == "pheno") {
      ## rename
      if (is.null(opt$pcol)) {
        int_index <- apply(data1, 2, function(x) all(floor(x) == x))
        opt$pcol <- which(!int_index)[1]
      } else if (opt$pcol > ncol(data1)) {
        stop("pcol (", opt$pcol, ") cannot be bigger than file columns (", ncol(data1), ")")
      }
      names(data1)[opt$pcol] <- names(data2)[opt$pcol] <- "phe"

      ## Missing value handling
      data1$phe[data1$phe == opt$miss] <- NA
      data2$phe[data2$phe == opt$miss] <- NA

      ## Calculate phenotype variance
      pvar1 <- var(data1$phe, na.rm = TRUE)
      pvar2 <- var(data2$phe, na.rm = TRUE)

      ## Obtain variance components
      vara1 <- pvar1 * opt$h2
      vara2 <- pvar2 * opt$h2B
      vare1 <- pvar1 * (1 - opt$h2)
      vare2 <- pvar2 * (1 - opt$h2B)
    } else if (opt$var == "predict") {
      if (opt$add_rf >= max(data1$V1)) stop("add_rf cannot be larger than the number of random groups in dmu!")
      if (any(data1$V2 != data2$V2)) stop("The random effects groups are different in the two files!")

      ## Obtain variance components
      vara1 <- data1$V4[opt$add_rf]
      vara2 <- data2$V4[opt$add_rf]
      vare1 <- tail(data1$V4, 1)
      vare2 <- tail(data2$V4, 1)
    }

    ## Calculate covariance
    cova <- opt$rg * sqrt(vara1) * sqrt(vara2)
    cove <- opt$re * sqrt(vare1) * sqrt(vare2)

    ## Ensure the positive definiteness of the covariance matrix
    ## Genetic variance
    vara <- matrix(c(vara1, cova, cova, vara2), 2, 2)
    vara_pd <- nearPD(vara, keepDiag = TRUE) # default
    if (!vara_pd$converged) {
      stop("warning! The positive definiteness of genetic covariance matrix cannot be guaranteed.\n")
    } else if (vara_pd$iterations > 1) {
      cat("add a small value to genetic effect covariance matrix\n")
      cova <- vara_pd$mat[2, 2]
    }
    ## Residual variance
    vare <- Matrix(c(vare1, cova, cova, vare2), 2, 2)
    vare_pd <- nearPD(vare, keepDiag = TRUE) # default
    if (!vare_pd$converged) {
      stop("warning! The positive definiteness of genetic covariance matrix cannot be guaranteed.\n")
    } else if (vare_pd$iterations > 1) {
      cat("add a small value to genetic effect covariance matrix\n")
      cova <- vare_pd$mat[2, 2]
    }

    ## Variance component form
    if ((opt$type == "multi") && (opt$method == "GBLUP")) {
      prior <- data.frame(
        group = c(1, 1, 1, 2, 2, 2),
        rindx = c(1, 2, 2, 1, 2, 2),
        cindx = c(1, 1, 2, 1, 1, 2),
        var = c(vara1, cova, vara2, vare1, cove, vare2)
      )
      ## Constrained residual covariance
      if (!is.null(opt$norec)) {
        prior <- prior[-5, ]
      }
    } else if (opt$method == "mbBayesAB") {
      if (!is.null(opt$rg_local)) {
        ## Check if the file exists
        if (!file.exists(opt$rg_local)) {
          cat(opt$rg_local, "not found!\n")
          quit(status = 1)
        }

        ## Obtain correlation coefficient
        rg_local <- fread(opt$rg_local)
        names(rg_local)[ncol(rg_local)] <- "cor"
        names(rg_local)[3] <- "nsnp"
        nsnp_total <- sum(rg_local$nsnp)

        ## Calculate the variance components of each block
        prior <- data.frame(
          var1 = vara1 * rg_local$nsnp / nsnp_total,
          var2 = vara2 * rg_local$nsnp / nsnp_total
        )
        prior$cov1 <- prior$cov2 <- rg_local$cor * sqrt(prior$var1) * sqrt(prior$var2)
        prior <- prior[, c("var1", "cov1", "cov2", "var2")]

        ## Add residual variance component
        prior[nrow(prior) + 1, ] <- c(vare1, cove, cove, vare2)

        ## Ensure the positive definiteness of the matrix
        for (i in seq_len(nrow(prior))) {
          vari <- matrix(unlist(prior[i, ]), 2, 2)
          vari_pd <- nearPD(vari, keepDiag = TRUE) # default
          if (!vari_pd$converged) {
            cat("warning! The positive definiteness of bin ", i, " matrix cannot be guaranteed.\n")
          } else if (vari_pd$iterations > 1) {
            cat("The covariance matrix of bin ", i, " is set to: \n", vari_pd$mat, "\n")
            prior[i, c("var1", "cov1", "cov2", "var2")] <- as.vector(vari_pd$mat)
          }
        }
      } else {
        prior <- data.frame(matrix(c(
          vara1, cova, cova, vara2,
          vare1, cove, cove, vare2
        ), 2, 4, byrow = TRUE))
      }
    } else if ((opt$type == "single") && (opt$method == "GBLUP")) {
      if (opt$var == "predict") {
        prior <- data1
        prior$V4 <- (prior$V4 + data2$V4) / 2
      } else {
        prior <- data.frame(
          group = c(1, 2),
          rindx = c(1, 1),
          cindx = c(1, 1),
          var = c(mean(vara1, vara2), mean(vare1, vare2))
        )
      }
    } else {
      stop("Unknown type!")
    }

    if (!is.null(opt$out)) {
      ## Output the result
      write.table(prior, out, row.names = FALSE, col.names = FALSE, quote = FALSE)
      # cat("prior file output to:", out, "\n")
    } else {
      prior
    }
  }
}

if (output) cat("variance component file generated.\n")
