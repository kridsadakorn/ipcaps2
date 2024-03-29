#' IPCAPS2 : Iterative Pruning to CApture Population Structure
#'
#' An unsupervised clustering algorithm based on iterative pruning is for capturing
#' population structure. This version supports ordinal data which can be applied
#' directly to SNP data to identify fine-level population structure and it is
#' built on the iterative pruning Principal Component Analysis (ipPCA) algorithm
#' (Intarapanich et al., 2009; Limpiti et al., 2011). The IPCAPS2 involves an
#' iterative process using multiple splits based on multivariate Gaussian
#' mixture modeling of principal components and Clustering EM estimation as explained in
#' Lebret et al. (2015). In each iteration, rough clusters and outliers are also
#' identified using the function rubikclust() from the R package \pkg{KRIS}.
#'
#' The R package \pkg{IPCAPS2} requires the package \pkg{KRIS}.
#'
#' Here is the list of functions in the R package \pkg{IPCAPS2}:
#' \itemize{
#' \item \code{\link{export.groups}}
#' \item \code{\link{get.node.info}}
#' \item \code{\link{ipcaps2}}
#' \item \code{\link{save.eigenplots.html}}
#' \item \code{\link{save.html}}
#' \item \code{\link{save.plots.cluster.html}}
#' \item \code{\link{save.plots.label.html}}
#' \item \code{\link{save.plots}}
#' \item \code{\link{top.discriminator}}
#' }
#'
#' Moreover, here is the list of example datasets in the R package \pkg{IPCAPS2}:
#' \itemize{
#' \item \code{\link{raw.data}}
#' \item \code{\link{label}}
#' \item \code{\link{PC}}
#' }
#' @keywords internal
#' @references
#'
#' Intarapanich, A., Shaw, P.J., Assawamakin, A., Wangkumhang, P., Ngamphiw, C.,
#' Chaichoompu, K., Piriyapongsa, J., and Tongsima, S. (2009). Iterative pruning
#' PCA improves resolution of highly structured populations. BMC Bioinformatics
#' 10, 382.
#'
#' Lebret, R., Iovleff, S., Langrognet, F., Biernacki, C., Celeux, G., and
#' Govaert, G. (2015). Rmixmod: TheRPackage of the Model-Based Unsupervised,
#' Supervised, and Semi-Supervised ClassificationMixmodLibrary. J. Stat. Softw.
#' 67.
#'
#' Limpiti, T., Intarapanich, A., Assawamakin, A., Shaw, P.J., Wangkumhang, P.,
#' Piriyapongsa, J., Ngamphiw, C., and Tongsima, S. (2011). Study of large and
#' highly stratified population datasets by combining iterative pruning
#' principal component analysis and structure. BMC Bioinformatics 12, 255.

"_PACKAGE"
#> [1] "_PACKAGE"

NULL

## quiets concerns of R CMD check
if (getRversion() >= "2.15.1")
    utils::globalVariables(c("thread_idx"))
