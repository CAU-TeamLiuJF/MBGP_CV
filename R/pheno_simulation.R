#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript

########################################################################################################################
## Version: 1.3.0
## Author:    Liweining liwn@jaas.ac.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-08-20
##
## Function：
## phenotype simulation
##
##
## Usage: ./pheno_simulation.R --gt "breedA" ...(Please refer to --help for detailed parameters)
##
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


## Get command-line parameters
spec <- matrix(
  c(
    "gt",       "A", 1, "character", "[Required] qtl plink file prefix of A",
    "h2",       "1", 1, "character", "[Optional] heritability of populations [0.5]",
    "rg",       "r", 1, "character", "[Optional] genetic correlation between breeds [0.2]",
    "win",      "w", 1, "integer",   "[Optional] number of snp in each bins [100]",
    "min",      "M", 1, "integer",   "[Optional] minimum SNPs in the bins where rg exist [nsnp_cor]",
    "bin",      "b", 1, "character", "[Optional] win / chr [win]",
    "binf",     "B", 1, "character", "[Optional] ld block file path [NULL]",
    "binc",     "C", 1, "integer",   "[Optional] column number of nSNP in block file [last col]",
    "miss",     "m", 1, "character", "[Optional] miss value in phenotype [-99]",
    "dist_cor", "c", 1, "character", "[Optional] distribution of genetic correlation of qtl [identical/uniform]",
    "nbin_cor", "q", 1, "integer",   "[Optional] [0.05]",
    "nsnp_cor", "Q", 1, "integer",   "[Optional] [10]",
    "nqtl",     "p", 1, "integer",   "[Optional] Proportion of QTL to total SNPs [0.01]",
    "seed",     "s", 1, "integer",   "[Optional] seed of random effect sampling",
    "mean",     "e", 1, "character", "[Optional] population mean of poopulations [1.0]",
    "out",      "o", 1, "character", "[Optional] output phenotype file name",
    "qtlf",     "t", 1, "character", "[Optional] output QTL effects file name",
    "overlap",  "v", 0, "logical",   "[Optional] allow QTLs and markers to overlap",
    "evenly",   "E", 0, "logical",   "[Optional] Distribute QTLs evenly within the selected bins",
    "fid",      "f", 0, "logical",   "[Optional] whether output family id",
    "help",     "h", 0, "logical",   "This is Help!"
  ),
  byrow = TRUE, ncol = 5
)
opt <- getopt::getopt(spec = spec)

## Check parameters
n_need <- length(c(opt$gt))
if (!is.null(opt$help) || n_need < 1) {
  cat(paste(getopt::getopt(spec = spec, usage = TRUE), "\n"))
  quit()
}

## Custom function
argp_parse <- function(char, num, sep = " ", label = "parameters") {
  argv <- as.numeric(unlist(strsplit(char, sep)))

  if (length(argv) == 1) {
    return(rep(argv, num))
  } else if (length(argv) != num) {
    cat("The number of", label, "should be:", num, "!\n")
    quit()
  } else {
    return(argv)
  }
}

find_local_extreme <- function(x, diff = 0.2, type = "min") { # nolint
  ## Find extreme points
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

    if (type == "min" && left_peak && right_peak) {
      final <- c(final, i)
    } else if (type == "max" && (left_peak || right_peak)) {
      final <- c(final, i)
    }
  }

  ## Exclude endpoints
  final <- final[final != 1 & final != n]

  return(final)
}

## Load required packages (install them in advance if not already installed)
cat("Loading required packages... \n")
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("MASS"))
suppressPackageStartupMessages(library("Matrix"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyr"))

## Default parameters
if (is.null(opt$seed)) opt$seed <- Sys.time()
if (is.null(opt$miss)) opt$miss <- -99
if (is.null(opt$h2)) opt$h2 <- "0.5"
if (is.null(opt$rg)) opt$rg <- "0.2"
if (is.null(opt$nbin_cor)) opt$nbin_cor <- 0.05 ## >1 means correlated regions, <1 means the proportion of correlated regions
if (is.null(opt$nqtl)) opt$nqtl <- 0.005 ## The proportion of QTLs to the total number of SNPs
if (is.null(opt$nsnp_cor)) opt$nsnp_cor <- 10
if (is.null(opt$min)) opt$min <- opt$nsnp_cor
if (is.null(opt$dist_cor)) opt$dist_cor <- "identical"
if (is.null(opt$mean)) opt$mean <- "1.0"
if (is.null(opt$bin)) opt$bin <- "win"
if (is.null(opt$win)) opt$win <- 100
if (!is.null(opt$binf)) opt$bin <- opt$binf

## Random seed
set.seed(as.integer(opt$seed))
cat("random number seed: ", opt$seed, "\n")

## Check if genotype files exist
bfiles <- trimws(opt$gt, which = c("both"), whitespace = "[ \t\r\n]")  ## Trim leading and trailing whitespace
bfiles <- unlist(strsplit(bfiles, " "))
for (prefix in bfiles) {
  for (suffix in c(".map", ".ped")) {
    f <- paste0(prefix, suffix)
    if (!file.exists(f)) {
      cat("file needed not found:", f, "\n")
      quit()
    }
  }
}

## Number of breeds
np <- length(bfiles)

## Number of upper triangular elements in the correlation matrix
nrg <- choose(np, 2)

## Parameter parsing
rgs <- argp_parse(opt$rg, nrg, label = "rg")
h2s <- argp_parse(opt$h2, np, label = "h2")
means <- argp_parse(opt$mean, np, label = "mean")

## Report heritability and other parameters
cat("\ncorrelations of QTL effect: ", rgs, "\n")
cat("heritabilities: ", h2s, "\n")
cat("population means: ", means, "\n")
cat("bins definition: ", opt$bin, "\n")
cat("bins with genetic correlation: ", opt$nbin_cor, "\n")
cat("number of SNPs with genetic correlation in each bins: ", opt$nsnp_cor, "\n")
if (opt$bin == "win") cat("nsnp in bins: ", opt$win, "\n")
cat("distribution of genetic correlation : ", opt$dist_cor, "\n\n")

## Load plink genotype files
cat("Loading plink genotype files... \n")
maps <- list()
peds <- list()
nsnps <- c()
for (i in 1:np) {
  maps[[i]] <- fread(paste0(bfiles[i], ".map"))
  peds[[i]] <- fread(paste0(bfiles[i], ".ped"))
  names(maps[[i]]) <- c("CHR", "SNP", "cM", "POS")
  nsnps <- c(nsnps, nrow(maps[[i]]))
}

## Shared markers
if (length(unique(nsnps)) > 1) {
  cat("The number of markers in the genotype file is inconsistent with that in the ldblock file!\n")
  quit()
}

## Generate SNP regions
map <- maps[[1]]
chrs <- unique(map$CHR)
if (opt$bin == "win") {
  ## Divide regions according to the provided win, where each interval contains win SNPs
  bin <- nsnp_bin <- NULL
  for (i in seq_along(chrs)) {
    nsnp_i <- sum(map$CHR == chrs[i])  ## Number of SNPs on chromosome i
    win_num <- floor(nsnp_i / opt$win) ## Number of regions that can be divided on chromosome i
    if (win_num < 1) win_num <- 1
    nsnp_bini <- rep(opt$win, win_num)
    nsnp_bini[win_num] <- nsnp_bini[win_num] + nsnp_i - sum(nsnp_bini)
    nsnp_bin <- c(nsnp_bin, nsnp_bini)
  }
  bin <- data.frame(bin = seq_len(length(nsnp_bin)), nsnp = nsnp_bin)
} else if (opt$bin == "chr") {
  ## Each chromosome is one interval, and the number of QTLs on each chromosome is equal
  bin <- as.data.frame(table(as.factor(map$CHR)))
  opt$nbin_cor <- nrow(bin)
  names(bin) <- c("bin", "nsnp")
} else if (file.exists(opt$bin)) {
  ## Select QTLs in each interval according to the provided interval division file
  bin <- fread(paste0(opt$bin))
  plink_ld_name <- c("CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2")
  if (all(names(bin) %in% plink_ld_name)) {
    ## Determine SNP regions based on local LD
    cat("calculating mean LD in windows...\n")
    map$index <- FALSE
    for (chr in chrs) {
      snps <- map$SNP[map$CHR == chr]
      ldi <- bin[bin$CHR_A == chr]
      mapi <- map[map$CHR == chr, ]
      for (j in seq_along(snps)) {
        index <- (max(j - opt$win, 1)):(min(j + opt$win, length(snps)))
        mapi$R2[j] <- mean(unlist(subset(ldi, SNP_A %in% snps[index] & SNP_B %in% snps[index], select = "R2")))
      }
      max <- find_local_extreme(mapi$R2, type = "max")
      map$index[map$CHR == chr & map$SNP %in% mapi$SNP[max]] <- TRUE
    }

    index <- unlist(sapply(which(map$index), function(x) (x - opt$win + 1):(x + opt$win + 1)))
    index <- index[index > 0 & index <= nrow(map)]
    map$index <- FALSE
    map$index[index] <- TRUE

    map$index <- c(FALSE, map$index[-length(map)])
    map$bin <- cumsum(!map$index)
  } else {
    if (is.null(opt$binc)) opt$binc <- ncol(bin)
    names(bin)[opt$binc] <- "nsnp"
    bin_nsnp_above <- which(bin$nsnp > opt$min)
  }
} else {
  cat("wrong parameter for bin:", opt$bin, "\n")
  quit()
}

## Indicate the interval each SNP belongs to
if ("bin" %in% names(map)) {
  cat("number of regions: ", max(map$bin), "\n")
} else {
  cat("number of regions: ", nrow(bin), "\n")
  map$bin <- rep(seq_len(nrow(bin)), times = bin$nsnp)
}

### Select QTL ###
## Select one marker locus for each interval
qmap <- map %>%
  subset(!is.na(bin)) %>%
  group_by(bin) %>%
  slice_sample()

## Randomly select nbin_cor regions
bins_sel <- c()
if (opt$nbin_cor > 0) {
  if (opt$nbin_cor <= 1) {
    ## If the provided parameter is a percentage, calculate the number of regions to select based on the percentage
    opt$nbin_cor <- floor(opt$nbin_cor * nrow(bin))
    cat("The number of bins with non-zero correlation is:", opt$nbin_cor, "\n")
  }

  ## regions where the number of SNPs meets the requirement
  if (!exists("bin_nsnp_above")) {
    bin_nsnp_above <- map$bin[duplicated(map$bin)]
  }

  ## regions where QTLs with genetic correlation can be selected
  bin_snps <- unique(map$bin[map$bin %in% bin_nsnp_above])

  ## Check if the number of regions is sufficient to meet the sampling requirements
  if (length(bin_snps) < opt$nbin_cor) {
    cat("The number of candidate regions cannot meet the sampling requirements\n")
    cat(length(bin_snps), "/", opt$nbin_cor, "\n")
    quit()
  }

  ## Select nbin_cor regions where QTLs are associated across breeds
  bins_sel <- sample(qmap$bin[qmap$bin %in% bin_snps], opt$nbin_cor)

  ## To ensure correlation within the interval, select an additional nsnp_cor-1 marker loci
  if (opt$nsnp_cor > 1) {
    if (is.null(opt$evenly)) {
      qmap2 <- subset(map, bin %in% bins_sel & !SNP %in% qmap$SNP) %>%
        subset(!is.na(bin)) %>%
        group_by(bin) %>%
        slice_sample(n = opt$nsnp_cor - 1)
      qmap <- bind_rows(qmap, qmap2)
    } else {
      for (i in seq_along(bins_sel)) {
        bini <- bins_sel[i]
        snp_count <- length(map$bin[map$bin == bini])

        # Calculate evenly distributed regions
        interval <- floor(snp_count / (opt$nsnp_cor))

        # Select QTLs at these regions
        qtl_indices <- seq(from = 1, to = snp_count, by = interval)
        qtl_indices <- qtl_indices[1:(opt$nsnp_cor)]
        if (length(qtl_indices) != opt$nsnp_cor) break

        # Select QTLs from the interval based on the selected indices
        qmap2 <- map %>%
          filter(bin == bini) %>%
          slice(qtl_indices)

        ## Remove the previously selected single SNP
        qmap <- subset(qmap, bin != bini)

        # Add the selected QTLs to the result dataframe
        qmap <- bind_rows(qmap, qmap2)
      }
    }
  }
}

## Number of QTLs
if (opt$nqtl < 1) {
  nqtl <- round(opt$nqtl * nrow(map))
} else {
  nqtl <- opt$nqtl
}

## Ensure the number of QTLs meets the requirement
over <- nrow(qmap) - nqtl
if (over > 0) {
  keep_single <- nqtl - sum(qmap$bin %in% bins_sel)
  keep_qtl <- sample(qmap$SNP[!qmap$bin %in% bins_sel], keep_single)
  qmap <- subset(qmap, SNP %in% keep_qtl | bin %in% bins_sel)
} else if (over < -0.5 * nrow(bin)) {
  over <- abs(over)
  add_each_bin <- ceiling(over / nrow(bin))
  qmap3 <- subset(map, !SNP %in% qmap$SNP) %>%
    subset(!is.na(bin)) %>%
    group_by(bin) %>%
    slice_sample(n = add_each_bin)
  qmap <- bind_rows(qmap, qmap3)
}

## Sort
qmap <- qmap[order(qmap$bin, qmap$POS), ]

## Report the number of QTLs
cat("total number of qtl: ", nrow(qmap), "\n")

## Extract QTL genotype information
qtl <- ebvs <- list()
fids <- c()
for (i in 1:np) {
  qtl_index <- which(maps[[i]]$SNP %in% qmap$SNP)
  col_index <- c(1:6, 2 * qtl_index + 6, 2 * qtl_index - 1 + 6)

  ## QTL genotype matrix
  qtl[[i]] <- subset(peds[[i]], select = col_index)
  ebvs[[i]] <- qtl[[i]][, 2:1]
  names(ebvs[[i]]) <- c("iid", "fid")

  ## breed ID
  fids <- c(fids, unique(qtl[[i]]$V1))
}

## Check the validity of fid
if (length(fids) != np) {
  cat("The 'fid' in the ped file is not a population identifier!\n")
  quit()
}

## Additive effect covariance matrix
cormat <- diag(np)
index <- 1
for (i in 1:(np - 1)) {
  for (j in (i + 1):np) {
    cormat[i, j] <- cormat[j, i] <- rgs[index]
    index <- index + 1
  }
}

## Sample QTL effects
if (opt$nbin_cor > 0) { # nolint
  ## Columns with genetic correlation
  rg_cols <- c()
  for (i in 1:(np - 1)) {
    for (j in (i + 1):np) {
      rg_cols <- c(rg_cols, paste(fids[i], fids[j], sep = "_"))
    }
  }

  ## Marker effects with genetic correlation
  if (opt$dist_cor == "identical") {
    qtl_cor <- which(qmap$bin %in% bins_sel)
    nsample <- length(qtl_cor)
    qmap[qtl_cor, fids] <- mvrnorm(nsample, rep(0, np), cormat)
    qmap[qtl_cor, rg_cols] <- matrix(rep(rgs, times = nsample), nsample, nrg, byrow = TRUE)
  } else {
    ## Genetic correlation in each region
    if (opt$dist_cor == "normal") {
      rg_bins <- rnorm(opt$nbin_cor, mean = 0, sd = 0.5)
    } else if (opt$dist_cor == "uniform") {
      rg_bins <- runif(opt$nbin_cor, -1, 1)
    } else {
      cat("dist_cor can only be identical, normal or uniform.\n")
      quit()
    }

    ## Sampling
    qmap[, fids] <- 0.00
    rgsi <- rep(0, nrg)
    for (i in seq_along(bins_sel)) {
      ## Marker effect covariance matrix
      cormati <- cormat
      ## Assign values
      index <- 1
      for (j in 1:(np - 1)) {
        for (k in (j + 1):np) {
          if (cormat[j, k] == 0) next
          rgij <- rg_bins[i] - mean(rg_bins) + cormat[j, k]

          ## Ensure validity
          if (rgij > 1) rgij <- 1.0
          if (rgij < -1) rgij <- -1.0

          rgsi[index] <- cormati[k, j] <- cormati[j, k] <- rgij
          index <- index + 1
        }
      }

      ## Ensure positive definiteness of the matrix
      cormat_pd <- nearPD(cormati, keepDiag = TRUE) # default
      if (!cormat_pd$converged) {
        stop("warning! The positive definiteness of genetic covariance matrix cannot be guaranteed.\n")
      } else if (cormat_pd$iterations > 1) {
        cat("add a small value to genetic effect covariance matrix\n")
        cormati <- cormat_pd$mat
      }

      ## Sample marker effects
      qtl_cor <- which(qmap$bin == bins_sel[i])
      nsample <- length(qtl_cor)
      qmap[qtl_cor, fids] <- matrix(mvrnorm(nsample, rep(0, np), cormati), nsample, np)

      ## Save correlation coefficients
      qmap[qtl_cor, rg_cols] <- matrix(rep(rgsi, times = nsample), nsample, nrg, byrow = TRUE)
    }
  }

  ## Marker effects without genetic correlation
  qtl_no_cor <- which(qmap[[fids[1]]] == 0.0 | is.na(qmap[[fids[1]]]))
  if (length(qtl_no_cor) > 0) {
    qmap[qtl_no_cor, fids] <- mvrnorm(length(qtl_no_cor), rep(0, np), diag(np))
    qmap[qtl_no_cor, rg_cols] <- 0
  }
}

cat("global genetic correlation is:\n")
cor(subset(qmap, select = fids))

## Reference allele (consistent across populations)
gene <- qtl[[1]][, 7:ncol(qtl[[1]])]
ref_allele <- unlist(subset(gene, c(TRUE, rep(FALSE, nrow(gene) - 1)),
    select = seq(1, ncol(gene), 2)))

## Generate phenotype information
cat("Creating phenotypes files... \n")
pheno <- data.frame()
for (i in seq_len(np)) {
  gene <- qtl[[i]][, 7:ncol(qtl[[i]])]

  ## If genotype is represented by letters, convert to 1, 2
  if (gene[1, 1] %in% LETTERS) {
    cat("Converting genotype data to 012 format...\n")

    ## Convert to 1, 2 based on the reference allele
    gene12 <- as.matrix(gene)
    for (j in 1:(ncol(gene) / 2)) {
      gene12[, 2 * j - 1] <- as.integer(subset(gene12, select = 2 * j - 1) == ref_allele[j]) + 1
      gene12[, 2 * j] <- as.integer(subset(gene12, select = 2 * j) == ref_allele[j]) + 1
    }

    ## Convert characters to numeric values
    gene <- apply(gene12, 2, as.numeric)
  }

  ## Convert gene content to 012 format
  # cat("Converting genotype data to 012 format...\n")
  gene <- sapply(
    seq(1, ncol(gene) - 1, 2),
    function(i) rowSums(gene[, i:(i + 1)]) - 2
  )

  ## Calculate additive effects
  ebvs[[i]]$tbv_raw <- apply(gene, 1, function(x) sum(x * qmap[[fids[i]]]))

  ## Normalize to standard normal distribution
  ebv_mean <- mean(ebvs[[i]]$tbv_raw)
  ebv_sd <- sd(ebvs[[i]]$tbv_raw)
  ebvs[[i]]$tbv_nor <- (ebvs[[i]]$tbv_raw - ebv_mean) / ebv_sd

  ## Residual effects (uncorrelated)
  ve <- ((1 / h2s[i]) - 1) * var(ebvs[[i]]$tbv_nor)
  ebvs[[i]]$envir <- rnorm(nrow(ebvs[[i]]), 0, sqrt(ve))

  ## Additive + Residual + Population mean = Phenotype
  ebvs[[i]]$phe <- ebvs[[i]]$tbv_nor + ebvs[[i]]$envir + means[i]

  ## Combine phenotypes of all breeds
  pheno <- rbind(pheno, ebvs[[i]])
}

## Add breed information and sort
if (is.null(opt$fid)) {
  pheno <- pheno[, c("iid", "breed", "phe", "tbv_raw", "tbv_nor", "envir")]
} else {
  ## Population (breed) effect identifier
  pheno$breed <- as.integer(as.factor(pheno$fid))
  pheno <- pheno[, c("iid", "fid", "breed", "phe", "tbv_raw", "tbv_nor", "envir")]
}

## Output phenotype file
if (is.null(opt$out)) opt$out <- paste(c(fids, "pheno.txt"), collapse = "_")
write.table(pheno, opt$out, row.names = FALSE, quote = FALSE)
cat("phenotypes output to file:", opt$out, "\n")

## True additive effects
if (is.null(opt$qtlf)) opt$qtlf <- paste(c(fids, "qtl.txt"), collapse = "_")
qmap <- qmap[order(qmap$bin), ]
write.table(qmap, opt$qtlf, row.names = FALSE, quote = FALSE, na = "0.0")
cat("true qtl effects output to file:", opt$qtlf, "\n")

## Overlap between markers and QTLs
if (is.null(opt$overlap) && !is.null(opt$bin)) {
  ## Update the number of markers in each region
  qtl_num <- as.data.frame(table(as.factor(qmap$bin)), stringsAsFactors = FALSE)
  qtl_num$Var1 <- as.integer(qtl_num$Var1)
  if (nrow(qtl_num) < nrow(bin)) {
    bins <- seq_len(nrow(bin))
    add_qtl_num <- data.frame(Var1 = bins[!bins %in% qmap$bin], Freq = 0)
    qtl_num <- bind_rows(qtl_num, add_qtl_num)
  }
  qtl_num <- qtl_num[order(qtl_num$Var1), ]
  bin$nsnp <- bin$nsnp - qtl_num$Freq
  newbinf <- paste0(basename(opt$bin), ".noqtl")
  write.table(bin$nsnp, newbinf, row.names = FALSE, quote = FALSE, col.names = FALSE)
  cat("new bins file output to:", newbinf, "\n")
}
