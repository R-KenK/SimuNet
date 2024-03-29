% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/edgeProb.R
\name{reconstruct_edgeProb}
\alias{reconstruct_edgeProb}
\title{Convenience function to reconstruct an \code{edgeProb} object from what has been stored in \code{scanList} object's attributes
Avoid \code{edgeProb} stored in \code{scanList} to store duplicates of \code{Adj}, \code{samp.effort}, \code{mode}, and \code{Adj.subfun}}
\usage{
reconstruct_edgeProb(scan.list)
}
\arguments{
\item{scan.list}{a \code{scanList} object (see \code{\link[=simunet]{simunet()}})}
}
\value{
an \code{edgeProb} object, i.e. a list containing:
\itemize{
\item \code{P}: the edge presence probability matrix
\item \code{Adj}: the inputted \code{Adj}
\item \code{samp.effort}: the inputted \code{samp.effort}
\item \code{mode}: the inputted \code{mode}
\item \code{Adj.subfun}: the inputted \code{Adj.subfun}
}
}
\description{
Convenience function to reconstruct an \code{edgeProb} object from what has been stored in \code{scanList} object's attributes
Avoid \code{edgeProb} stored in \code{scanList} to store duplicates of \code{Adj}, \code{samp.effort}, \code{mode}, and \code{Adj.subfun}
}
\keyword{internal}
