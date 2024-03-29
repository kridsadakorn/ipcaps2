% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.mode.R
\name{clustering.mode}
\alias{clustering.mode}
\title{(Internal function) Select a clustering method to be used for the IPCAPS2
process.}
\usage{
clustering.mode(
  node,
  work.dir,
  method,
  silence.mode = FALSE,
  seed = NULL,
  no.rep = 5,
  cutoff.confident = 0.5
)
}
\arguments{
\item{node}{An integer representing the current node number which is being
analyzed.}

\item{work.dir}{A working directory.}

\item{method}{A clustering method selected from the \code{\link{ipcaps2}}
function. See \code{\link{ipcaps2}} for available methods.}

\item{silence.mode}{To enable or disable silence mode. If silence mode is
enabled, the fuction is processed without printing any message on the screen,
and it is slightly faster. Default = TRUE.}

\item{seed}{To specify a seed number for random generator. Default = NA,
which means that a seed number is automatically chose.}

\item{no.rep}{To specify a  number of time to replicate the internal clustering. Default = 5,}

\item{cutoff.confident}{To specify a cutoff for confident values. The confident values
are calculated based the clustering results of all replicates. The confident values
represent the best of average values (agreed clusters) among the member of clusters.
Default = 0.5.}
}
\value{
A vector of cluster assignment, for which cluster each individual
belongs.
}
\description{
(Internal function) Select a clustering method to be used for the IPCAPS2
process.
}
\seealso{
\code{\link{ipcaps2}}
}
