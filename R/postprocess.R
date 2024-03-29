

#' (Internal function) Perform the post-processing step of IPCAPS2
#'
#' @param result.dir A result directory as the \code{$output} object returned
#' from the \code{\link{ipcaps2}} function.
#' @param reanalysis (Unimplemented) To specify whether it is re-analysis or not. If TRUE, it is
#' re-analysis, otherwise it is not. Default = FALSE.
#'
#' @return A data frame of clustering result containing 4 columns;
#' \code{group}, \code{node}, \code{label}, \code{row.number}, as described
#' below for more details:
#' \itemize{
#' \item \code{group} represents group membership of IPCAPS2 result.
#' \item \code{node} represents node numbers of IPCAPS2 result.
#' \item \code{label} represents labels of rows in orginal input data.
#' \item \code{row.number} represents row numbers of orginal input data.
#' }
#'
#' @include export.groups.R
#' @include save.plots.R
#' @include save.html.R
#'
postprocess <- function(result.dir, reanalysis = FALSE)
{
  no.plot <- NULL
  silence.mode <- NULL

  file.name <- file.path(result.dir, "RData", "leafnode.RData")
  load(file = file.name)
  file.name <- file.path(result.dir, "RData", "condition.RData")
  load(file = file.name)

  if (length(leaf.node) == 0)
  {
    leaf.node <- c(1)
    save(leaf.node, file = file.name, compress = "bzip2")
  }

  # Generate HTML output file
  cluster.tab <-
    export.groups(result.dir, silence.mode = silence.mode)
  save.html(result.dir)
  if (no.plot == FALSE)
  {
    save.plots(result.dir)
  }
  if (!silence.mode)
    cat("The result files were saved at: ", result.dir, "\n")
  return(cluster.tab)

}
