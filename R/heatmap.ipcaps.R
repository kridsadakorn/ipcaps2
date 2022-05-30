#' Plot heatmap from the result of IPCAPS
#'
#' @description Create a heatmap plot base on the top differentiators among groups
#' and the input data.
#'
#' @param ipcaps_output the returned object from the \code{\link{ipcaps2}} function.
#' @param groups the number of groups or clusters to be used in a heatmap.
#' Default = NULL.
#' @param min.member.in.group a minimum number of member in the groups that will
#' be used for top-marker selection and will be used for heatmap. Default = 10.
#' @param top.fst a number to indicate the top markers according to Fst.
#' The value can be set from 1 until the number of input markers. Defalt = 100.
#' @param fst.threshold a threshold to cut off for Fst. First, fst.threshold is considered,
#' then top.fst is considered. When the number of filtered markers is less than top.fst,
#' top.fst is ignored. The value can be set between 0 and 1. Default = 0.
#' @param is.cluster.snp the logic to indicate whether the rows of heatmap will be clustered or not.
#' If is.cluster.snp is TRUE, the rows will be clustered to K groups and K represents
#' the number of IPCAPS results by ignoring outliers (set by min.member.in.group).
#' Default = TRUE.
#'
#' @return \code{pheatmap}
#'
#' @details This function generates a graphical object which is compatible with
#' the function ggsave from the R package ggplot2. Therefore you can save a plot to
#' file via ggplot2::ggsave() from the return object. In order to see more
#' parameters for ploting, also see pheatmap::pheatmap().
#'
#' @import pheatmap
#' @import RColorBrewer
#'
#' @export
#'
#' @examples
#'
#' # Load the required packages
#' library(ggplot2)
#' library(gridExtra)
#'
#' # Importantly, bed file, bim file, and fam file are required
#' # Use the example files embedded in the package
#'
#' BED.file <- system.file("extdata","ipcaps_example.bed",package="IPCAPS2")
#' LABEL.file <- system.file("extdata","ipcaps_example_individuals.txt.gz",package="IPCAPS2")
#'
#' my.cluster <- ipcaps2(bed=BED.file,
#'                       label.file=LABEL.file,
#'                       lab.col=2,
#'                       out=tempdir(),
#'                       silence=TRUE,
#'                       no.rep=1)
#'
#' # Show a plot
#' my.plot <- heatmap.ipcaps(my.cluster)
#'
#' # To display gaps after each cluster, we can check the distribution of clustering result
#' table(my.cluster$cluster$group)
#' #1  2  3  4  5
#' #1  1  1 50 50
#' # You can calculate what to put the gabs by creating a vector for the location,
#' # by adding up the numbers according to the order of IPCAPS2_cluster in the heatmap
#' gaps_position = c(50, 100, 101, 102)
#' my.plot <- heatmap.ipcaps(my.cluster, gaps_col = gaps_position)
#'
#' # Save the plot to file.
#' ggplot2::ggsave("heatmap.pdf", my.plot)
#'
#' # To create multiple heatmap in one plot, for example
#' my.plot1 <- heatmap.ipcaps(my.cluster, main = "The first heatmap")
#' my.plot2 <- heatmap.ipcaps(my.cluster, main = "The second heatmap")
#' combined.plot <- gridExtra::grid.arrange(my.plot1$gtable,
#'                                          my.plot2$gtable,
#'                                          nrow=2)
#'
#' # To save the combined plots to file
#' ggplot2::ggsave("combined_heatmap.pdf", combined.plot)

heatmap.ipcaps <- function(ipcaps_output, groups = NULL,
                           min.member.in.group = 10,
                           top.fst = 1000,
                           fst.threshold = 0.01,
                           is.cluster.snp = TRUE,
                           ...)
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

    if ( (fst.threshold < 0) || (fst.threshold > 1) ){
        cat(paste0("Error: the parameter 'fst.threshold' is not correct.\n",
                   "The value needs to be between 0 and 1\n"))
        return(NULL)
    }

    # PCs is stored in node1.RData
    PCs = NULL
    load(file.path(ipcaps_output$output.dir,'RData','rawdata.RData'))

    if ( (top.fst < 1) || (top.fst > ncol(raw.data)) ){
        cat(paste0("Error: the parameter 'top.fst' is not correct.\n",
                   "The value needs to be between 1 and the number of input markers\n"))
        return(NULL)
    }



    # check for frequencies for clusters
    cluster_freq = as.data.frame(table(ipcaps_output$cluster$group))

    # filter out group with low members
    cluster_freq = cluster_freq[which(cluster_freq$Freq > min.member.in.group),]

    # select top SNPs based on Fst
    iter_group = cluster_freq$Var1
    if (length(iter_group) <2){
        print(paste0("Cannot create a heatmap because there is only one target cluster\ns"))
        return(NULL)
    }
    snp_to_plot = NULL
    for (i in 1:(length(iter_group)-1)){
        for (j in (i+1):length(iter_group)){
            # Here is the output from ipcaps_output$cluster
            #     group node    label row.number
            #1       1    3 outlier3        101
            #2       2    4 outlier3        102

            idx.p1 = ipcaps_output$cluster[which(ipcaps_output$cluster$group == iter_group[i]), c('row.number')]
            idx.p2 = ipcaps_output$cluster[which(ipcaps_output$cluster$group == iter_group[j]), c('row.number')]

            fst = fst.each.snp.hudson(raw.data,
                                      idx.p1 = idx.p1,
                                      idx.p2 = idx.p2)

            fst = as.data.frame(fst)
            snp_df = cbind(snp.info, fst)
            snp_df = snp_df[order(snp_df$fst, decreasing = TRUE),]

            # filter Fst
            snp_df = snp_df[which(snp_df$fst >= fst.threshold),]

            # in case that top.fst less than the number of filtered SNPs.
            if (nrow(snp_df) < top.fst){
                top.fst = nrow(snp_df)
            }
            snp_df = snp_df[1:top.fst,]
            snp_to_plot = rbind(snp_to_plot, snp_df)
        }
    }
    # in case that the top SNPs are duplicated, so check only the unique SNPs
    snp_to_plot_unique = unique(snp_to_plot$ID)

    # get SNP matrix to plot
    colnames(raw.data) = snp.info[,2]
    rownames(raw.data) = ind.info[,2]
    to_plot_df = t(raw.data[,which(snp.info$ID %in% snp_to_plot_unique)])

    annotation_col = ipcaps_output$cluster
    rownames(annotation_col) = annotation_col$row.number
    annotation_col = annotation_col[,c('group','label')]
    colnames(annotation_col) = c('IPCAPS2_cluster','label')
    annotation_col$IPCAPS2_cluster = as.character(annotation_col$IPCAPS2_cluster)

    # scale the values
    to_plot_df = t(scale(t(to_plot_df), center = TRUE, scale = TRUE))

    phm = NULL
    if (is.cluster.snp == TRUE){
        phm = pheatmap(to_plot_df,
                       cluster_rows = TRUE,
                       cluster_cols = FALSE,
                       annotation_col = annotation_col,
                       show_rownames = TRUE,
                       show_colnames = TRUE,
                       legend = TRUE,
                       kmeans_k = length(iter_group),
                       color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                       ...)
    }else{
        phm = pheatmap(to_plot_df,
                       cluster_rows = TRUE,
                       cluster_cols = FALSE,
                       annotation_col = annotation_col,
                       show_rownames = TRUE,
                       show_colnames = TRUE,
                       legend = TRUE,
                       color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                       ...)
    }

    return(phm)
}

