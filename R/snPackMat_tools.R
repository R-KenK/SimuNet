# SimuNet packed matrix tools ---------------------------------------------

#' snPackMat (SimuNet Packed Matrix) object constructor
#'
#' @param M a regular matrix to pack
#' @param Adj.subfun Subsetting function of the adjacency matrix. Driven by
#'   igraph `"mode"` argument.
#' @param mode Character scalar, specifies what type of igraph network `mode`
#'   should be used to convert the supplied matrix. Ignored if `sampling.param` is
#'   provided. Possible values are:
#' \itemize{
#'   \item{`"directed"` (the default): for non-symmetrical adjacency matrix where
#'   `Adj[i,j]` doesn't have the same meaning as `Adj[j,i]`}
#'   \item{`"undirected"`: same as `"max"`}
#'   \item{`"upper"`: undirected matrix focusing only on the upper triangle of
#'   `Adj` (relying on `upper.tri`). Either `"upper"` or `"lower"` could be
#'   favor if only one of `Adj[i,j]` and `Adj[j,i]` should be randomized}
#'   \item{`"lower"`: undirected matrix focusing only on the lower triangle of
#'   `Adj` (relying on `lower.tri`)}
#'   \item{`"max"`: from a `"directed"` randomization process (both `Adj[i,j]`
#'   and `Adj[j,i]` will be drawn at each scan), `max(Adj[i,j],Adj[j,i])` will
#'   be kept for both}
#'   \item{`"min"`: from a `"directed"` randomization process (both `Adj[i,j]`
#'   and `Adj[j,i]` will be drawn at each scan), `min(Adj[i,j],Adj[j,i])` will
#'   be kept for both}
#'   \item{`"plus."`:  from a `"directed"` randomization process (both
#'   `Adj[i,j]` and `Adj[j,i]` will be drawn at each scan), `Adj[i,j] +
#'   Adj[j,i]` will be kept for both}
#'   \item{`"vector"`: experimental. To consider adjacency matrices as flat
#'   vectors to be randomized. Relevance unexplored yet.}
#'   \item{See details \link[igraph]{graph_from_adjacency_matrix}}
#' }
#'
#' @return a `snPackMat` object with:
#' \itemize{
#'   \item{`M` a vector containing the data from M}
#'   \item{`mode` character string indicating chosen igraph's `mode`}
#'   \item{`n` integer, number of nodes in the graph}
#'   \item{`node_names` character vector (or `NULL`) containing the names of the nodes.}
#' }
#' @noRd
generate_snPackMat <- function(M,Adj.subfun,mode) {
  M <- list(
    M = M[Adj.subfun(M)],
    mode = mode,
    n = nrow(M),
    node_names = rownames(M)
  )
  class(M) <- "snPackMat"
  M
}

#' Wrapper for the snPackMat constructor
#' Choose the right mode
#'
#' @param M a regular matrix to pack
#' @param mode Character scalar, specifies what type of igraph network `mode`
#'   should be used to convert the supplied matrix. Ignored if `sampling.param` is
#'   provided. Possible values are:
#' \itemize{
#'   \item{`"directed"` (the default): for non-symmetrical adjacency matrix where
#'   `Adj[i,j]` doesn't have the same meaning as `Adj[j,i]`}
#'   \item{`"undirected"`: same as `"max"`}
#'   \item{`"upper"`: undirected matrix focusing only on the upper triangle of
#'   `Adj` (relying on `upper.tri`). Either `"upper"` or `"lower"` could be
#'   favor if only one of `Adj[i,j]` and `Adj[j,i]` should be randomized}
#'   \item{`"lower"`: undirected matrix focusing only on the lower triangle of
#'   `Adj` (relying on `lower.tri`)}
#'   \item{`"max"`: from a `"directed"` randomization process (both `Adj[i,j]`
#'   and `Adj[j,i]` will be drawn at each scan), `max(Adj[i,j],Adj[j,i])` will
#'   be kept for both}
#'   \item{`"min"`: from a `"directed"` randomization process (both `Adj[i,j]`
#'   and `Adj[j,i]` will be drawn at each scan), `min(Adj[i,j],Adj[j,i])` will
#'   be kept for both}
#'   \item{`"plus."`:  from a `"directed"` randomization process (both
#'   `Adj[i,j]` and `Adj[j,i]` will be drawn at each scan), `Adj[i,j] +
#'   Adj[j,i]` will be kept for both}
#'   \item{`"vector"`: experimental. To consider adjacency matrices as flat
#'   vectors to be randomized. Relevance unexplored yet.}
#'   \item{See details \link[igraph]{graph_from_adjacency_matrix}}
#' }
#'
#' @return a `snPackMat` object with:
#' \itemize{
#'   \item{`M` a vector containing the data from M}
#'   \item{`mode` character string indicating chosen igraph's `mode`}
#'   \item{`n` integer, number of nodes in the graph}
#'   \item{`node_names` character vector (or `NULL`) containing the names of the nodes.}
#' }
#' @noRd
pack_snPackMat <- function(M,mode) {
  switch(mode,
         "upper" = generate_snPackMat(M,upper.tri,mode),
         "lower" = generate_snPackMat(M,lower.tri,mode),
         "undirected" = ,
         "max" = ,
         "min" = ,
         "plus" = ,
         "directed" = ,
         generate_snPackMat(M,non.diagonal,mode)
  )
}

#' Reconstruct a regular matrix from data stored in snPackMat object
#'
#' @param M an snPackMat object
#'
#' @return a regular matrix
#' @noRd
unpack_snPackMat <- function(M) {
  if (!is.snPackMat(M)) {return(M)}
  switch(M$mode,
         "upper" = {
           unpacked <- matrix(
             0L,nrow = M$n,ncol = M$n,
             dimnames = list(M$node_names,M$node_names)
           )
           unpacked[upper.tri(unpacked)] <- M$M
           unpacked
         },
         "lower" = {
           unpacked <- matrix(
             0L,nrow = M$n,ncol = M$n,
             dimnames = list(M$node_names,M$node_names)
           )
           unpacked[lower.tri(unpacked)] <- M$M
           unpacked
         },
         "undirected" = ,
         "max" = ,
         "min" = ,
         "plus" = ,
         "directed" = ,
         {
           unpacked <- matrix(
             0L,nrow = M$n,ncol = M$n,
             dimnames = list(M$node_names,M$node_names)
           )
           unpacked[non.diagonal(unpacked)] <- M$M
           unpacked
         }
  )
}

#' Print method for `snPackMat` objects
#' @export
#' @noRd
print.snPackMat <- function(x,...) {
  use_printSpMatrix(unpack_snPackMat(x),...)
}

#' Operators method for `snPackMat` objects
#' @export
#' @noRd
Ops.snPackMat <- function(e1,e2) {
  e1.snPM <- is.snPackMat(e1)
  e2.snPM <- is.snPackMat(e2)
  is.boolean <- switch(.Generic,`<` = , `>` = , `==` = , `!=` = ,
                       `<=` = , `>=` = TRUE, FALSE)
  if (e1.snPM & e2.snPM) {
    ans <- e1
    e1 <- e1$M
    e2 <- e2$M
  } else if (e1.snPM) {
    ans <- e1
    e1 <- e1$M
  } else {
    ans <- e2
    e2 <- e2$M
  }

  ans$M <- NextMethod()

  if (is.boolean) {
    matrix(as.logical(unpack_snPackMat(ans)),
           nrow = ans$n,ncol = ans$n,
           dimnames = list(ans$node_names,ans$node_names)
    )
  } else {
    ans
  }
}

#' Math functions method for `snPackMat` objects
#' @export
#' @noRd
Math.snPackMat <- function(x,...) {
  ans <- x
  x <- x$M # functional, but does not retain snPackMat class (e1 + e2 becomes a "matrix" "array')
  ans$M <- NextMethod()
  ans
}

#' Summary functions method for `snPackMat` objects
#' @export
#' @noRd
Summary.snPackMat <- function(...,na.rm = FALSE) {
  args <- list(...)
  args <- lapply(args, function(x) {
    x <- x$M
  })
  do.call(.Generic, c(args, na.rm = na.rm))
}

#' Subsetting functions method for `snPackMat` objects
#' @export
#' @noRd
`[.snPackMat` <- function(x,i,j = NULL){
  if (!is.null(j)) {
    unpack_snPackMat(x)[i,j]
  } else if (length(dim(i)) == 2L) {
    unpack_snPackMat(x)[i]
  } else {
    x$M[i]
  }
}

#' Subsetting functions method for `snPackMat` objects
#' @export
#' @noRd
`[<-.snPackMat` <- function(x,i,j = NULL,value){
  if (!is.null(j)) {
    mode <- x$mode
    x <- unpack_snPackMat(x)
    x[i,j] <- value
    pack_snPackMat(x,mode)
  } else {
    x$M[i] <- value
    x
  }
}



#' Test if object if a `snPackMat` object
#'
#' @param x an object to test.
#'
#' @return logical, `TRUE` if the inputted object is a `snPackMat` object.
#'
#' @noRd
is.snPackMat <- function(x) {
  inherits(x, "snPackMat")
}

#'‘Not Available’ / Missing Values of a `snPackMat` object
#' @export
#' @noRd
# is.na.snPackMat <- function(x) {
#   is.na(x$M)
# }
