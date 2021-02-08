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
         M
  )
}

#' Print method for `snPackMat` objects
#' @export
#' @noRd
print.snPackMat <- function(x,...) {
  print(unpack_snPackMat(x),...)
}

#' Operators method for `snPackMat` objects
#' @export
#' @noRd
Ops.snPackMat <- function(e1,e2) {
  e1 <- unpack_snPackMat(e1) # functional, but does not retain snPackMat class (e1 + e2 becomes a "matrix" "array')
  e2 <- unpack_snPackMat(e2)
  NextMethod()
}

#' Math functions method for `snPackMat` objects
#' @export
#' @noRd
Math.snPackMat <- function(x,...) {
  x <- unpack_snPackMat(x) # functional, but does not retain snPackMat class (e1 + e2 becomes a "matrix" "array')
  NextMethod()
}

#' Summary functions method for `snPackMat` objects
#' @export
#' @noRd
Summary.snPackMat <- function(...,na.rm = FALSE) {
  args <- list(...)
  args <- lapply(args, function(x) {
    x <- unpack_snPackMat(x) # functional, but does not retain snPackMat class (e1 + e2 becomes a "matrix" "array')
  })
  do.call(.Generic, c(args, na.rm = na.rm))
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

