#' Plot 3 views from the result of IPCAPS
#'
#' @description Plot the first three components of PCA using in 3 views.
#'
#' @param ipcaps_output the returned object from the \code{\link{ipcaps2}} function.
#'
#' @return \code{NULL}
#'
#' @details This function generates a plot in the standard R graphic. Therefore,
#' it can be saved to a file using png(), pdf(),  jpeg(), and tiff().
#'
#' @export
#'
#' @import KRIS
#'
#' @examples
#'
#' # Importantly, bed file, bim file, and fam file are required
#' # Use the example files embedded in the package
#'
#' BED.file <- system.file("extdata","ipcaps_example.bed",package="IPCAPS2")
#' LABEL.file <- system.file("extdata","ipcaps_example_individuals.txt.gz",package="IPCAPS2")
#'
#' my.cluster <- ipcaps2(bed=BED.file,label.file=LABEL.file,lab.col=2,out=tempdir(),silence=TRUE,no.rep=1)
#'
#' # Show a plot
#' plot3views.ipcaps(my.cluster)

plot3views.ipcaps <- function(ipcaps_output)
{
    # check if the parameter is a correct format
    if ( typeof(ipcaps_output) != "list"){
        cat(paste0("Error: the parameter 'ipcaps_output' is not correct.\n",
            "Use the returned value from the function 'ipcaps2' as an input to this function\n"))
        return(NULL)
    }

    if ( length(ipcaps_output) != 2){
        cat(paste0("Error: the parameter 'ipcaps_output' is not correct.\n",
                   "Use the returned value from the function 'ipcaps2' as an input to this function\n"))
        return(NULL)
    }

    load(file.path(ipcaps_output$output.dir,'RData','node1.RData'))
    KRIS::plot3views(X = PCs,
                     labels = ipcaps_output$cluster[order(ipcaps_output$cluster$row.number),]$group)
}
