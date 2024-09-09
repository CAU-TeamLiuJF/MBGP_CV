#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript


########################################################################################################################
## Version: 1.3.0
## Author:    Liweining liwn@jaas.ac.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-08-20
##
## Function：
## Used to calculate the accuracy of cross validation
##
##
## Usage: ./accuracy_bias_calculation.R --dir_all "/path/to/phenotype/file" ...(Please refer to --help for detailed parameters)
##
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


# Load the required packages
suppressPackageStartupMessages(library("getopt"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dplyr"))

## Command line parameters
spec <- matrix(
  c(
    "dir_val", "R", 1, "character", "[Required] reduce(val) dataset parameter card prefix/valX",
    "phe_all", "P", 1, "character", "[Required] Full dataset phenotype file name",
    "tbvf",    "B", 1, "character", "[Required] file contain true breeding value",
    "ebvf",    "e", 1, "character", "[Required] file contain Estimated breeding value",
    "dir_all", "A", 1, "character", "[Optional] Full dataset parameter card prefix",
    "label",   "l", 1, "character", "[Optional] breed label",
    "famf",    "f", 1, "character", "[Optional] fam file whose ID order is consistent with that in EBV file",
    "ebvidc",  "i", 1, "integer",   "[Optional] id columns in Estimated breeding value file",
    "id_col",  "C", 1, "character", "[Optional] id columns in phenotype file [1]",
    "ebv_id",  "c", 1, "integer",   "[Optional] id columns in true breeding value file [1]",
    "fold",    "F", 1, "integer",   "[Optional] fold/1",
    "val_idf", "W", 1, "character", "[Optional] Validation ids file/valX_id_2col.txt",
    "ebv_col", "t", 1, "integer",   "[Optional] trait rank [1]",
    "add_rnd", "n", 1, "integer",   "[Optional] Random effects group number for additive genetic effects [1]",
    "tbv_col", "T", 1, "integer",   "[Optional] Position of TBV in real variable columns of phe_all",
    "add_sol", "a", 1, "integer",   "[Optional] code for type of additive effect in sol first column",
    "rep",     "r", 1, "integer",   "[Optional] number of traits [4]",
    "digit",   "d", 1, "integer",   "[Optional] digit of output accuracy [4]",
    "fid",     "I", 1, "character", "[Optional] use to index valid in column 1 of vidf  [4]",
    "out",     "O", 1, "character", "[Optional] output file name prefix/NULL\n",
    "rmNeg",   "g", 0, "logical",   "[Optional] rm negative value in accuracy",
    "mean",    "m", 0, "logical",   "[Optional] output mean of cor",
    "help",    "h", 0, "logical",   "This is Help!"
  ),
  byrow = TRUE, ncol = 5
)
opt <- getopt(spec = spec)

## Check parameters
if (!is.null(opt$help) || (is.null(opt$phe_all) && is.null(opt$tbvf)) ||
    (is.null(opt$dir_val) && is.null(opt$ebvf))) {
  cat(paste(getopt(spec = spec, usage = TRUE), "\n"))
  quit()
}

## Default parameters
if (is.null(opt$add_sol)) opt$add_sol <- 4
if (is.null(opt$fold)) opt$fold <- 1
if (is.null(opt$rep)) opt$rep <- 1
if (is.null(opt$digit)) opt$digit <- 5
if (is.null(opt$dir_val)) opt$dir_val <- "./val#val#/rep#rep#/dmu"
if (is.null(opt$val_idf)) opt$val_idf <- "./val#val#/rep#rep#/val.id"
if (is.null(opt$tbvf) && !is.null(opt$tbv_col)) opt$tbvf <- opt$phe_all
if (is.null(opt$id_col)) {
  opt$id_col <- 1
} else {
  opt$id_col <- as.numeric(opt$id_col)
}


##############################################################################
##################  True breeding value/corrected phenotype  #################
##############################################################################
if (!is.null(opt$tbvf)) {
  if (file.exists(opt$tbvf)) {
    if (is.null(opt$tbv_col)) opt$tbv_col <- 2
    y_adj <- fread(opt$tbvf)
    names(y_adj)[opt$id_col] <- "id"
    names(y_adj)[opt$tbv_col] <- "yhat"
  } else {
    cat(opt$tbvf, "not found!\n")
    quit()
  }
} else if (!is.null(opt$dir_all)) {
  ## Obtain corrected phenotypes from estimating EBV and residue from a complete dataset
  phe_all <- fread(opt$phe_all)
  names(phe_all)[opt$id_col] <- "id"
  sol_all <- fread(paste0(opt$dir_all, ".SOL"))
  res_all <- fread(paste0(opt$dir_all, ".RESIDUAL"))

  ## Extract ebv (in within-breed prediction)
  ebv_all <- subset(sol_all, V1 == opt$add_sol & V2 == 1, c(5, 8))
  names(ebv_all) <- c("id", "ebv")

  ## Extract residuals (the corrected phenotype used for multi-trait models is calculated
  ## based on within-breed prediction, so here is the residual file for single-trait models)
  re <- subset(res_all, select = c(1, 4))
  names(re) <- c("rows", "re")
  re$id <- phe_all$id[re$rows]

  ## Calculate the corrected phenotype
  y_adj <- left_join(re, ebv_all, by = "id")
  y_adj$yhat <- y_adj$ebv + y_adj$re

  ## Output corrected phenotype
  if (!is.null(opt$tbvf)) {
    fwrite(y_adj[, c("id", "yhat")], opt$tbvf, col.names = FALSE, sep = "\t")
  }
}


## Accurate result storage
accs <- data.frame(
  rep = rep(1:opt$rep, each = opt$fold),
  fold = rep(1:opt$fold, times = opt$rep),
  cor = NA, bias = NA, rank_cor = NA, num = NA
)

## Calculation of accuracy
for (r in 1:opt$rep) { # nolint r=1;f=1
  for (f in 1:opt$fold) {
    if (!is.null(opt$ebvf)) {
      ## File name
      ebvf <- gsub("#val#", f, opt$ebvf)
      ebvf <- gsub("#rep#", r, ebvf)

      ## Check if the file exists
      if (!file.exists(ebvf)) {
        ## There is no result file for this situation, skip it
        accs <- accs[!(accs$rep == r & accs$fold == f), ]
        cat("warning:", ebvf, "does not exist!\n")
        next
      } else if (file.info(ebvf)$size == 0) {
        accs <- accs[!(accs$rep == r & accs$fold == f), ]
        next
        cat("warning:", ebvf, "has a size of 0!\n")
      }

      ## load estimated breeding values
      ebv_val <- fread(ebvf)
      if (!is.null(opt$famf)) {
        ## Column where ebv is located
        if (is.null(opt$ebv_col)) opt$ebv_col <- 1
        names(ebv_val)[opt$ebv_col] <- "ebv_val"

        ## filename
        famf <- gsub("#val#", f, opt$famf)
        famf <- gsub("#rep#", r, famf)

        ## The estimated breeding value ID in Bayes is consistent with the ID in the fam file
        ebv_id <- fread(famf)
        ebv_id <- ebv_id[, 1:2]
        names(ebv_id) <- c("fid", "id")
        ebv_val <- cbind(ebv_id, ebv_val)
      } else {
        if (is.null(opt$ebv_id)) opt$ebv_id <- 1
        if (is.null(opt$ebv_col)) opt$ebv_col <- 2
        names(ebv_val)[opt$ebv_id] <- "id"
        names(ebv_val)[opt$ebv_col] <- "ebv_val"
      }
    } else {
      ## DMU result file name
      dir_val <- gsub("#val#", f, opt$dir_val)
      dir_val <- gsub("#rep#", r, dir_val)
      solf <- paste0(dir_val, ".SOL")

      ## Check if the DMU result file is normal
      if (!file.exists(solf)) {
        cat(solf, "not found!\n")
        accs <- accs[!(accs$rep == r & accs$fold == f), ]
        next
      } else if (file.info(solf)$size == 0) {
        cat(solf, "zero size!\n")
        accs <- accs[!(accs$rep == r & accs$fold == f), ]
        next
      }

      ## trait codes in SOL
      if (is.null(opt$ebv_col)) opt$ebv_col <- 1

      ## Which group of random effects is the genetic effect in SOL
      if (is.null(opt$add_rnd)) opt$add_rnd <- 1

      ## Extract the EBV of the validation
      sol_val <- fread(solf)
      ebv_val <- subset(sol_val, V1 == opt$add_sol & V2 == opt$ebv_col & V4 == opt$add_rnd, c(5, 8))
      names(ebv_val) <- c("id", "ebv_val")
    }

    ## Match reference population ebv
    y_adj_ebv <- inner_join(y_adj, ebv_val, by = "id")

    ## reference individual's id
    val_idf <- gsub("#val#", f, opt$val_idf)
    val_idf <- gsub("#rep#", r, val_idf)
    if (file.exists(val_idf)) {
      val_id <- fread(val_idf)
    } else {
      cat(val_idf, " not found!\n")
      next
    }

    ## rename
    if (ncol(val_id) > 1) {
      names(val_id)[1:2] <- c("fid", "id")
      if (!is.null(opt$fid)) {
        val_id <- subset(val_id, fid == opt$fid)
      }
    } else {
      names(val_id) <- "id"
    }

    ## Filter out individuals in the validation group
    ids_index <- y_adj_ebv$id %in% val_id$id
    if (sum(ids_index) < 1) {
      cat("warn: the number of validation individuals in val", f, "rep", r, "is 0, please check\n")
      next
    }
    y_adj_ebv_val <- subset(y_adj_ebv, ids_index)

    ## check the standard deviation
    if (sd(y_adj_ebv_val$ebv_val) == 0) {
      cat("warn: the standard deviation of EBVs in val", f, "rep", r, "is 0, please check\n")
      next
    }

    ## Accuracy (Pearson correlation coefficient)
    accur_rf <- cor(y_adj_ebv_val$yhat, y_adj_ebv_val$ebv_val)

    ## Unbiasedness (regression coefficient)
    bias_rf <- lm(yhat ~ ebv_val, data = y_adj_ebv_val)$coefficients[2]

    ## Rank correlation (Spearman correlation coefficient)
    rank_cor <- cor(y_adj_ebv_val$yhat, y_adj_ebv_val$ebv_val, method = "spearman")

    ## Save results
    accs$cor[accs$rep == r & accs$fold == f] <- round(rank_cor, opt$digit)
    accs$bias[accs$rep == r & accs$fold == f] <- round(bias_rf, opt$digit)
    accs$rank_cor[accs$rep == r & accs$fold == f] <- round(accur_rf, opt$digit)
    accs$num[accs$rep == r & accs$fold == f] <- sum(ids_index)
  }
}

## Drop the NA line
accs <- subset(accs, !is.na(accs$cor))

if (any(!is.na(accs$cor))) {
  ## Check if there are negative value in the results
  if (any(accs$cor < 0)) {
    if (is.null(opt$rmNeg)) {
      cat("Negative values founded in result. Provide the --rmNeg parameter can delete these values\n")
    } else {
      cat("delete", sum(accs$cor < 0), "values in results\n")
      accs <- subset(accs, cor >= 0)
    }
  }

  ## mean
  mean_cor <- paste(round(mean(accs$cor, na.rm = TRUE), opt$digit),
    round(sd(accs$cor, na.rm = TRUE), opt$digit),
    sep = "±"
  )
  mean_bias <- paste(round(mean(accs$bias, na.rm = TRUE), opt$digit),
    round(sd(accs$bias, na.rm = TRUE), opt$digit),
    sep = "±"
  )
  mean_rank_cor <- paste(round(mean(accs$rank_cor, na.rm = TRUE), opt$digit),
    round(sd(accs$rank_cor, na.rm = TRUE), opt$digit),
    sep = "±"
  )

  mean_cor <- data.frame(
    rep = "mean",
    fold = NA,
    cor = mean_cor,
    bias = mean_bias,
    rank_cor = mean_rank_cor,
    num = mean(accs$num, na.rm = TRUE)
  )
  accs$cor <- as.character(accs$cor)
  accs$bias <- as.character(accs$bias)
  accs$rep <- as.character(accs$rep)
  accs$rank_cor <- as.character(accs$rank_cor)

  accs <- bind_rows(accs, mean_cor)

  ## Print the results on the screen
  print(accs)

  ## Output result file
  if (!is.null(opt$out)) fwrite(accs, opt$out, sep = "\t")
} else {
  cat("No outcome document in all cases!\n")
}
