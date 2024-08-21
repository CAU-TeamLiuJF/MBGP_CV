#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript

########################################################################################################################
## Version: 1.3.0
## Author:    Liweining liwn@cau.edu.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-08-20
##
## Function：
## Check if all the required R packages have been installed
##
##
## Usage: ./package_required.R
##
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


## List of R packages to be checked
packages <- c("data.table", "LAVA", "getopt", "MASS", "RcppEigen", "Matrix", "ggplot2", "reshape2",
  "pedigree", "RcppArmadillo", "dplyr", "coda", "lattice", "stringr", "Cairo", "Rcpp", "reshape")

## Check and install missing packages
missing_packages <- packages[!packages %in% installed.packages()]
if (length(missing_packages) > 0) {
  message("The following R packages are not installed: ")
  cat(missing_packages, sep = "\n")
  message("Please run the following command in R language to install packages: ")
  command <- paste0("install.package(\"", paste(missing_packages, collapse = "\", \""), "\")")
  cat(command, "\n")
}
