# SimuNet packed matrix tools ---------------------------------------------


#' TO WRITE
#'
#' @param M TO WRITE
#' @param mode TO WRITE
#'
#' @return TO WRITE
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

#' TO WRITE
#'
#' @param M TO WRITE
#' @param mode TO WRITE
#'
#' @return TO WRITE
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
         M
  )
}

#' TO WRITE
#'
#' @param M TO WRITE
#'
#' @return TO WRITE
#' @noRd
unpack_snPackMat <- function(M) {
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
         M
  )
}

#' Print method for `snPackMat` objects
#' @export
#' @noRd
print.snPackMat <- function(x,...) {
  print(unpack_snPackMat(x),...)
}

#' Print method for `snPackMat` objects
#' @export
#' @noRd
Ops.snPackMat <- function(e1,e2) {
  M1 <- unpack_snPackMat(e1)
  M2 <- unpack_snPackMat(e2)
  Ops.(M1,M2)
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

