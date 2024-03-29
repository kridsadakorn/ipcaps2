% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.node.info.R
\name{get.node.info}
\alias{get.node.info}
\title{Get the information for specified node}
\usage{
get.node.info(cluster.obj, node)
}
\arguments{
\item{cluster.obj}{A list returned from the \code{\link{ipcaps2}} function.}

\item{node}{An integer representing a node number to enquire information as
shown in the HTML output files.}
}
\value{
The return value is \code{NULL} if node's information does not exist or a list containing
\code{PCs}, \code{eigen.fit}, \code{index}, and \code{label} as explained
below:
\itemize{
\item \code{$PCs} is a matrix of pricipal components of this node.
\item \code{$eigen.fit} is a number represeting the EigenFit value of this
node.
\item \code{$index} is a vector of row number (individuals) of \code{raw.data} (input data).
\item \code{$label} is the vector of labels of all individuals that belongs
to this node.
}
}
\description{
Obtain the information for specified
node from the output list of \code{\link{ipcaps2}}.
}
\examples{

# Importantly, bed file, bim file, and fam file are required
# Use the example files embedded in the package

BED.file <- system.file("extdata","ipcaps_example.bed",package="IPCAPS2")
LABEL.file <- system.file("extdata","ipcaps_example_individuals.txt.gz",package="IPCAPS2")

my.cluster <- ipcaps2(bed=BED.file,
                      label.file=LABEL.file,
                      lab.col=2,
                      out=tempdir(),
                      silence=TRUE,
                      no.rep=1)

#Here, to obtain the information of specified node, for example, node 3
node.info <- get.node.info(my.cluster,3)
ls(node.info)
}
