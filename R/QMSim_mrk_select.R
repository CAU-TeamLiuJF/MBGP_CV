#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript

########################################################################################################################
## Version: 1.3.0
## Author:    Liweining liwn@cau.edu.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-08-20
##
## Function：
## Select a certain number of tags from QMSim tags based on the condition settings
##
##
## Usage: ./QMSim_mrk_select.R --freqf "/path/to/freq/file" ...(Please refer to --help for detailed parameters)
##
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


## Load packages required
# cat("Loading required packages... \n")
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("getopt"))

## Command-line parameters
spec <- matrix(
  c("freqf",    "F", 2, "character",  "[Required] gene frequency files from QMSim, eg. './%pop%_freq_mrk_001.txt'",
    "popN",     "A", 2, "character",  "[Required] populations name, e.g. 'breedA breedB' ",
    "mapf",     "M", 2, "character",  "[Required] map file from QMSim",
    "gen",      "G", 2, "character",  "[Optional] generation for obtaining gene frequency",
    "nsel",     "N", 2, "integer",    "[Optional] number of markers needed [50000]",
    "seg",      "S", 2, "character",  "[Optional] one/both [one]",
    "binDiv",   "D", 2, "character",  "[Optional] pos/frq [NULL]",
    "binThr",   "T", 2, "double",     "[Optional] frequency or distance (Mbp) [0.1]",
    "binPdf",   "P", 2, "character",  "[Optional] eql/nor [eql]",
    "maf",      "m", 2, "double",     "[Optional] minor allele frequency [NULL]",
    "frqDif",   "d", 2, "character",  "[Optional] high/low [NULL]",
    "difp",     "p", 2, "character",  "[Optional] [0.14]",
    "dif2pq",   "q", 2, "character",  "[Optional] [1.0]",
    "format",   "f", 2, "character",  "[Optional] id/index [index]",
    "out",      "o", 2, "character",  "[Optional] output file name [mrk.txt]",
    "outmapf",  "a", 2, "character",  "[Optional] [NULL]",
    "inOut",    "n", 2, "character",  "[Optional] include [in] or exclude [ex]",
    "preSel",   "e", 2, "character",  "[Optional] [NULL]",
    "preSelC",  "c", 2, "integer",    "[Optional] [1]",
    "Mapf",     "i", 2, "character",  "[Optional] map file 2 from QMSim",
    "Out",      "I", 2, "character",  "[Optional] output file name 2",
    "fixPer",   "x", 2, "integer",    "[Optional] percentage [NULL]",
    "seed",     "s", 2, "integer",    "[Optional] [NULL]",
    "plotfrq",  "l",  2, "character", "[Optional] plot scatter plot of gene frequency",
    "cid",      "C",  0, "logical",   "[Optional] change ID to chr_pos",
    "help",     "h", 0, "logical",    "This is Help!"),
  byrow = TRUE, ncol = 5)
opt <- getopt(spec = spec)

## Check parameters
if (!is.null(opt$help) || is.null(opt$freqf) || is.null(opt$mapf)) {
  cat(paste(getopt(spec = spec, usage = TRUE), "\n"))
  quit()
}

## Default parameters
if (is.null(opt$nsel))     opt$nsel <- 50000
if (is.null(opt$nbin))     opt$nbin <- 500
if (is.null(opt$binThr))   opt$binThr <- 0.1
if (is.null(opt$binPdf))   opt$binPdf <- "eql"
if (is.null(opt$difp))     opt$difp <- "0.05"
if (is.null(opt$dif2pq))   opt$dif2pq <- "0.1"
if (is.null(opt$format))   opt$format <- "index"
if (is.null(opt$out))      opt$out <- "mrk.txt"
if (is.null(opt$inOut))    opt$inOut <- "ex"

## Random Number Seed
if (is.null(opt$seed)) {
  opt$seed <- Sys.time()
  write.table(opt$seed, paste0(out, ".seed"))
}
set.seed(opt$seed)

## Functions
insufficient_snps <- function(left, need, step) {
  if (left < need) {
    cat("Insufficient number of SNP (", left, ") after '", step, "' processing!\n", sep = "")
    quit(status = 2)
  }
}

## Mark coordinate file
if (!file.exists(opt$mapf)) {
  cat("Map file '", opt$mapf, "' does not exist!\n", sep = "")
  quit(status = 2)
}
map <- fread(opt$mapf)
names(map) <- c("ID", "Chr", "Pos")
# ## Change coordinates from cM to physical location(x 1e6)
# if (max(map$Pos < 1000)) {
#   names(map)[3] <- "cM"
#   map$Pos <- map$cM * 1e6
# }
## Change marker ID
if (!is.null(opt$cid)) {
  map$ID <- paste(map$Chr, map$Pos * 1e6, sep = "_")
}
## NUmber of SNP
snp_total <- nrow(map)
cat("Number of SNPs before filtering:", snp_total, "\n")

## breed name extraction
pop_name <- unlist(strsplit(opt$popN, " "))
np <- length(pop_name)

## Gene frequency related size
opt$difp <- as.numeric(unlist(strsplit(opt$difp, " ")))
if (length(opt$difp) == 1) opt$difp <- rep(opt$difp, choose(np, 2))
opt$dif2pq <- as.numeric(unlist(strsplit(opt$dif2pq, " ")))
if (length(opt$dif2pq) == 1) opt$dif2pq <- rep(opt$dif2pq, choose(np, 2))

## Read gene frequency file
for (i in 1:np) {
  cat("Processing the gene frequency file of the", i, "population...\n")

  ## File name
  freqf <- gsub("%pop%", pop_name[i], opt$freqf)

  ## Check if the file exists
  if (!file.exists(freqf)) {
    cat(paste0("file not found: ", freqf, "\n"))
    quit()
  }

  ## Read the file
  frqi <- fread(freqf, fill = TRUE)
  names(frqi) <- c("ID", "Gen", "Chr", "Allele_Freq", "Freq2")

  ## Change ID
  if (!is.null(opt$cid)) {
    frqi$ID <- paste(frqi$Chr, frqi$Pos * 1e6, sep = "_")
  }

  ## Source of gene frequency generation
  if (is.null(opt$gen)) {
    opt$gen <- max(frqi$Gen)
  }

  ## Extract the gene frequency of a specified generation
  frqi <- subset(frqi, Gen == opt$gen)

  ## Analyzing gene frequency from QMSim file
  ref_frq <- str_split_fixed(frqi$Allele_Freq, ":", 2)
  frqi$ref <- as.numeric(ref_frq[, 1])
  frqi$frq <- as.numeric(ref_frq[, 2])

  datai <- frqi[, c("ID", "Chr", "ref", "frq")]
  names(datai)[3:4] <- c(paste0("ref_", pop_name[i]), paste0("frq_", pop_name[i]))

  ## Merge
  if (i == 1) {
    datas <- datai
    reffirst <- paste0("ref_", pop_name[1])
    frqfirst <- paste0("frq_", pop_name[1])
  } else {
    datas <- left_join(datas, datai, by = c("ID", "Chr"))

    ## Ensure that the reference bases of the two breeds are consistent (with breed A as the reference), 
    ## and note that only biallelic genes are assumed to exist here by default
    refnow <- paste0("ref_", pop_name[i])
    frqnow <- paste0("frq_", pop_name[i])

    not_same <- datas[[refnow]] != datas[[reffirst]]
    if (sum(not_same) > 0) {
      cat("Ensure that reference allele (", sum(not_same), ") are consistent...\n")
      datas[[refnow]][not_same] <-  datas[[reffirst]][not_same]
      datas[[frqnow]][not_same] <-  1 - datas[[frqnow]][not_same]
    }
    datas <- subset(datas, select = !grepl(refnow, colnames(datas)))
  }
}

## Filter tags based on the marker ID provided by the user
if (!is.null(opt$preSel)) {
  if (is.null(opt$preSelC)) opt$preSelC <- 2
  exc <- fread(opt$preSel)
  names(exc)[opt$preSelC] <- "ID"
  include_idx <- datas$ID %in% exc$ID
  if (opt$inOut == "ex") include_idx <- !include_idx

  ## Check if the number of filtered marker is sufficient
  insufficient_snps(sum(include_idx), opt$nsel, "user specified exclusion")

  ## filter
  datas <- datas[include_idx, ]
  cat(sum(!include_idx), "variants removed due to user specified exclusion.\n")
}

## Delete intermediate files
rm(datai, frqi, ref_frq)

## Gene frequency column
frqcols <- paste0("frq_", LETTERS[1:np])

## Match marker coordinates
datas <- left_join(datas, map, by = c("ID", "Chr"))

## Mean gene frequency
cat("Calculating average gene frequency...\n")
datas$frq_mean <- apply(datas[, (ncol(datas) - np + 1):(ncol(datas))], 1, mean)

## Screening based on the separation status of markers in different breeds
if (!is.null(opt$seg)) {
  if (opt$seg == "one") {
    include_idx <- datas$frq_mean > 0 & datas$frq_mean < 1
  } else if (opt$seg == "both") {
    cat("Markers are being screened according to the segregation of alleles...")
    include_idx <- apply(subset(datas, select = frqcols), 1, function(x) all(x > 0 & x < 1))
  } else {
    cat("seg can only be one or both, please provide correct parameters.\n")
    quit()
  }

  ## Check if the number of filtered SNP is sufficient
  insufficient_snps(sum(include_idx), opt$nsel, "segregation status")

  ## filter
  datas <- datas[include_idx, ]
  cat(sum(!include_idx), "variants removed due to allele segregation.\n")
}

## Filter based on MAF
if (!is.null(opt$maf)) {
  include_idx <- apply(subset(datas, select = frqcols), 1, function(x) all(x > opt$maf & x < (1 - opt$maf)))

  ## Check if the number of filtered SNPs is sufficient
  insufficient_snps(sum(include_idx), opt$nsel, "maf")

  ## filter
  datas <- datas[include_idx, ]
  cat(sum(!include_idx), "variants removed due to minor allele threshold(s).\n")
}

## Filter based on gene frequency differences
if (!is.null(opt$frqDif)) {
  index <- 1
  for (i in 1:np) {
    for (j in 1:np) {
      if (i >= j) next
      ## Location
      frqi <- paste0("frq_", pop_name[i])
      frqj <- paste0("frq_", pop_name[j])

      ## Correlation of frequencies before filtering
      corij <- cor(datas[[frqi]], datas[[frqj]])
      pop_name_ij <-  paste0(pop_name[i], " and ", pop_name[j])
      cat("correlation of the gene frequencies between", pop_name_ij, "before filtering:", corij, "\n")

      ## Difference magnitude (absolute value)
      frq_mean <- (datas[[frqi]] + datas[[frqj]]) / 2
      difp <- abs(datas[[frqi]] - datas[[frqj]])
      pqi <- 2 * datas[[frqi]] * (1 - datas[[frqi]])
      pqj <- 2 * datas[[frqj]] * (1 - datas[[frqj]])
      pqij <- 2 * frq_mean * (1 - frq_mean)
      dif2pq <- abs(pqi - pqj) / pqij
      dif2pq[pqij == 0] <- 0

      ## Markers that meet the criteria
      if (opt$frqDif == "high") {
        diffp <- difp >= opt$difp[index]
        diff2pq <- dif2pq >= opt$dif2pq[index]
      } else if (opt$frqDif == "low") {
        diffp <- difp < opt$difp[index]
        diff2pq <- dif2pq < opt$dif2pq[index]
      } else {
        cat("frqDif can only be high or low, please provide correct parameters.\n")
        quit() 
      }

      ## Markers that meet the criteria
      include_idx <- diffp & diff2pq

      ## Report
      cat("number of markers excluded due to gene frequency between ", pop_name_ij, ":", sum(!include_idx), "\n")
      # cat("number of markers excluded due to gene (p) frequency between ", popNij, ":", sum(!diffp), "\n")
      # cat("number of markers excluded due to genotype (2pq) frequency between ", popNij, ":", sum(!diff2pq), "\n")

      if (i == 1 && j == 2) {
        include_idxs <- include_idx
      } else {
        include_idxs <- include_idxs & include_idx
      }

      ## Check if the number of markers after filtering is sufficient
      insufficient_snps(sum(include_idxs), opt$nsel, "gene frequency difference")

      index <- index + 1
    }
  }
  ## Filtering
  datas <- datas[include_idxs, ]
  cat(sum(!include_idxs), "variants removed due to gene frequency difference.\n")
}

cat("correlation of the gene frequencies between breeds after filtering:\n")
cor(subset(datas, select = frqcols))

## Assign candidate markers to different regions
if (!is.null(opt$binDiv)) {
  cat("assigning markers to different bins...\n")
  if (opt$binDiv == "pos") {
    chrs <- unique(datas$Chr)
    star <- 0
    bin <- c()
    for (i in chrs) {
      pos_chri <- subset(datas, Chr == i, select = "Pos")
      nbin <- (max(pos_chri$Pos) - min(pos_chri$Pos)) %/% opt$binThr
      bini <- as.integer(cut(pos_chri$Pos, breaks = nbin)) + star
      bin <- c(bin, bini)
      star <- star + nbin
    }

    ## Total number of regions
    nbin <- star

    ## Assign regions
    if (length(bin) == nrow(datas)) {
      datas$bins <- bin
    } else {
      cat("error in bins determination!\n")
      quit()
    }

    cat("The number of bins:", nbin, "\n")
  } else if (opt$binDiv == "frq") {
    nbin <- (max(datas$frq_mean) - min(datas$frq_mean)) %/% opt$binThr
    datas$bins <- cut(datas$frq_mean, breaks = nbin, labels = 1:nbin)
    cat("The number of bins:", nbin, "\n")
  } else {
    cat("binDiv can only be pos or frq, please provide correct parameters.\n")
    quit()
  }
} else {
  datas$bins <- 1
}

## Number of candidate markers
cat("Number of candidate SNPs:", nrow(datas), "\n")

## Select candidate markers in different regions
if (opt$binPdf == "eql") {
  ## Number of markers needed in each region
  nsnp_bini <- opt$nsel %/% length(unique(datas$bins)) + 1

  ## Check if the number of markers in each region is sufficient for sampling
  bin_count <- table(datas$bins)
  insuffix <- sum(bin_count < nsnp_bini)
  if (insuffix > 0) {
    if (!is.null(opt$fixPer)) {
      if (insuffix <= opt$fixPer / 100 * nbin) {
        cat("The number of markers in some bins is not enough, adjust the bins.\n")

        ## Regions with insufficient markers
        insuff_bins <- as.numeric(names(bin_count)[bin_count < nsnp_bini])

        ## Select markers from adjacent regions
        allow_miss <- opt$fixPer / 100 * nsnp_bini
        for (i in insuff_bins) {
          bini_insuffix <- nsnp_bini - sum(datas$bins == i)
          if (bini_insuffix < allow_miss) {
            index <- c()
            if (i != 1) index <- c(index, which(datas$bins == (i - 1)))
            if (i != nbin) index <- c(index, which(datas$bins == (i + 1)))
            add <- sample(index, bini_insuffix)
            datas$bins[add] <- i
          } else {
            cat("bin ", i, " is missing more than ", opt$fixPer, "% of markers ",
              sum(datas$bins == i), "/", nsnp_bini, "\n",
              sep = ""
            )
            quit()
          }
        }
      } else {
        cat("More than ", opt$fixPer, "% of bins (", insuffix, ") with insufficient labels!\n", sep = "")
        quit()
      }

      ## Recheck after adjustment
      bin_count <- table(datas$bins)
      insuffix <- sum(bin_count < nsnp_bini)
      if (insuffix > 0) {
        cat("The number of markers in", insuffix, "bins are less than the number of samples after adjusting!\n")
        quit()
      }
    } else {
      cat("The number of markers in", insuffix, "bins are less than the number of samples!\n")
      quit()
    }
  }

  ## Sampling within regions
  mrk_sel <- datas %>%
    group_by(bins) %>%
    slice_sample(n = nsnp_bini)
} else if (opt$binPdf == "nor") {
  cat("This parameter is not supported now!\n")
  quit()
} else {
  cat("binPdf can only be eql or nor, please provide correct parameters.\n")
  quit()
}

## Output map file
options(scipen = 10) ## Prevent scientific notation when outputting SNP physical positions
if (!is.null(opt$outmapf)) {
  map_out <- mrk_sel[, c("Chr", "ID", "Pos")]
  map_out$pos <- map_out$Pos * 1e6

  ## Sorting
  map_out <- map_out[order(map_out$Chr, map_out$Pos), ]

  ## Handle markers with the same position
  pos_dup <- duplicated(map_out[, c("Chr", "Pos")])
  if (any(pos_dup)) {
    map_out$Pos[pos_dup] <- map_out$Pos[pos_dup] + 1
  }

  ## Output file
  write.table(map_out, opt$outmapf, row.names = FALSE, col.names = FALSE, quote = FALSE)
}

## Output file format
if (opt$format == "id") {
  out <- mrk_sel$ID
} else if (opt$format == "index") {
  out <- as.numeric(map$ID %in% mrk_sel$ID)
} else {
  cat("format can only be id or index, please provide correct parameters.\n")
  quit()
}

## Gene frequency scatter plot
if (!is.null(opt$plotfrq)) {
  # png(paste0(opt$plotfrq, ".png"), width = 600 * 1.5, height = 600, bg = "white")
  CairoPNG(paste0(opt$plotfrq, ".png"), width = 600 * 1.5, height = 600, bg = "white")
  plot(mrk_sel$frq_A, mrk_sel$frq_B)
  hide_message <- dev.off()
  cat("gene frequency plot output to:", paste0(opt$plotfrq, ".png"), "\n")
}

## Report
cat("The number of markers finally selected:", sum(out), "\n")
cat("correlation of the gene frequencies between breeds:\n")
cor(subset(datas, select = frqcols))

## Output selected markers
write.table(out, opt$out, row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("File output to:", opt$out, "\n")

## Index for another breed
if (!is.null(opt$Mapf)) {
  out2 <- as.numeric(map2$ID %in% mrk_sel$ID)
  write.table(out2, opt$Out, row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat("File 2 output to:", opt$Out, "\n")
}
