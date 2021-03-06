% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output.template.R
\docType{data}
\name{output.template}
\alias{output.template}
\title{(Internal object) The HTML output template for IPCAPS2}
\format{
A list with with 5 objects
}
\usage{
output.template
}
\description{
\code{output.template} contains \code{$lno_title}, \code{$lno_data},
\code{$lno_leafnode}, \code{$lno_body}, and \code{$template} as explained
below:
\describe{
\item{lno_title}{A index number of \code{$template} that a title of HTML
file is replaced.}
\item{lno_data}{A index number of \code{$template} that a data section in
the HTML file is replaced.}
\item{lno_leafnode}{A index number of \code{$template} that a leaf-node
section in the HTML file is replaced.}
\item{lno_body}{A index number of \code{$template} that a body section in
the HTML file is replaced.}
\item{template}{A vector of characters for HTML file.}
}
}
\keyword{output.template}
