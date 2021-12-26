#' Synthetic dataset containing single nucleotide polymorphisms (SNP)
#'
#' The raw.data is the simulated dataset which consists of 2,000 independent SNPs
#' and 103 individuals belonging to one of 2 populations (50 individuals
#' each) and 3 outlying individuals. The matrix \code{raw.data} contains the number 0, 1,
#' and 2 representing SNP in additive coding.
#'
#' @name raw.data
#' @docType data
#' @format A matrix with 2,000 columns and 103 rows
#' @seealso \code{\link{label}} and \code{\link{PC}}
#' @usage data(ipcaps_example)
#' @keywords raw.data
#' @md
#' @references
#' Balding, D.J., and Nichols, R.A. (1995). A method for quantifying
#' differentiation between populations at multi-allelic loci and its
#' implications for investigating identity and paternity. Genetica 96, 3-12.
"raw.data"





#' Synthetic dataset containing population labels for the dataset \code{raw.data}
#'
#' A dataset contains a character vector of 103 elements containing labels or
#' populations of 103 individuals which they belong. Two populations and
#' outliers were labeled as "pop1", "pop2", and "outlier3".
#'
#' @name label
#' @docType data
#' @format A vector with 103 elements.
#' @seealso \code{\link{raw.data}} and \code{\link{PC}}
#' @usage data(ipcaps_example)
#' @keywords label
"label"





#' Synthetic dataset containing the top 10 principal components (PC) from the
#' dataset \code{raw.data}
#'
#' A dataset contains a numeric matrix of 103 rows and 10 columns of top 10
#' PCs calculated from the dataset \code{raw.data}. The PCs were calculated
#' using linear principal component analysis (PCA), see more datails at
#' \code{KRIS::cal.pc.linear}
#'
#' @name PC
#' @docType data
#' @format A matrix with 10 columns and 103 rows
#' @seealso \code{\link{raw.data}} and \code{\link{label}}
#' @usage data(ipcaps_example)
#' @keywords PC
"PC"

NULL
