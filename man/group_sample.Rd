% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expDesign_sampling.R
\name{group_sample}
\alias{group_sample}
\title{Performs a group-scan sampling over a \code{scanList} object
Internal.}
\usage{
group_sample(
  scan.list,
  sampling = c("constant", "matrix", "random", "function"),
  all.sampled = TRUE
)
}
\arguments{
\item{scan.list}{a \code{scanList} object}

\item{sampling}{for \verb{method = "group}, users should input either:
\itemize{
\item a numeric scalar (\code{"constant"}): the constant probability of observing an edge for all edges
\item a numeric matrix (\code{"matrix"}): the probabilities of observing an edge for each edges
\item a character scalar: for common sampling regimes:
\itemize{
\item \code{"random"}: random edge observation probabilities
}
\item a user-defined function (\code{"function"}): a function of the adjacency matrix \code{Adj} (can be
named anything) that returns a matrix of the probabilities of observing an edge for each edges
}}

\item{all.sampled}{logical scalar, should all nodes be sampled at least once? (TO CHECK: does it
work with group-scan sampling?)}
}
\value{
an empirical \code{scanList} object in which, compared to the \code{theoretical.scanList} (added
to \code{attrs}), unobserved edges are replaced by \code{NA}s (regardless of them being 0 or 1).

Returned \code{scanList} has new attributes added to attrs:
\itemize{
\item \code{obs.P}: matrix of probabilities of observing an edge (whether it is 0 or 1)
\item \code{theoretical.scanList}: the original theoretical \code{scanList} from which some edges have not been
observed
}
}
\description{
Performs a group-scan sampling over a \code{scanList} object
Internal.
}
\seealso{
\code{\link[=design_sampling]{design_sampling()}}, \code{\link[=determine_obsProb]{determine_obsProb()}}.
}
\keyword{internal}
