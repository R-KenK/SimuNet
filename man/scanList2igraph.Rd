% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scanList_matList_conversion.R
\name{scanList2igraph}
\alias{scanList2igraph}
\title{Convert scanList 3D array into \code{igraph} object (with \code{attrs})}
\usage{
scanList2igraph(scan.list)
}
\arguments{
\item{scan.list}{a 3 dimensional array representing adjacency matrices (first 2 dimensions)
throughout the different scans (3rd dimension)}
}
\value{
an \code{igraphSN} object:
\itemize{
\item the \code{igraph} network object obtained from the weighted adjacency matrix corresponding to the
inputted \code{scanList}
\item the \code{attrs} attributes list carried over from the inputted \code{scanList}
}
}
\description{
To back-transform an \code{igraphSN} object into a \code{scanList} object, users can use
\code{igraph2scanList()}, wrapper to \code{\link[=igraph2array]{igraph2array()}} that
handles attributes list \code{attrs}
}
\examples{
set.seed(42)
n <- 5L
samp.effort <- 100L

# Adjacency matrix import
## random directed adjacency matrix
Adj <- sample(1:samp.effort,n * n) |>
  matrix(nrow = 5,dimnames = list(letters[1:n],letters[1:n]))
Adj[lower.tri(Adj,diag = TRUE)] <- 0L
Adj

sL <- simunet(Adj,samp.effort,"upper",10)
G <- sL |> scanList2igraph()
G
G |> igraph2scanList()
G |> igraph2scanList(format = "scanList")
identical(G |> igraph2scanList(format = "scanList"),sL)
}
