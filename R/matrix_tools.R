# Matrix subset tools -----------------------------------------------------

#' Non Diagonal Part of a Matrix
#'
#' Similarly to `upper.tri` and `lower.tri`, returns a matrix of `logical`s to identify the diagonal of a square matrix
#'
#' @param M square matrix.
#' @param output choose to either return a matrix of logical values (`TRUE`s, `FALSE`s on the diagonal, i.e. comparable behavior than `upper.tri` for instance) or a vector of the subsetted values of the matrix.
#'
#' @return either return a matrix of logical values (`TRUE`s, `FALSE`s on the diagonal) or a vector of the subsetted values of the matrix.
#' @export
#'
#' @examples
#' M<- matrix(sample(1:10,16,replace = TRUE),4,4)
#' non.diagonal(M)
non.diagonal <- function(M,output=c("matrix.logical","vector.values")) {
  output <- match.arg(output)

  if (dim(M)[1] == dim(M)[2]) {
    logicals <- upper.tri(M,diag = FALSE) |
      lower.tri(M,diag = FALSE)
  } else stop("Matrix provided is not a square matrix.")

  switch(output,
         "matrix.logical" = logicals,
         "vector.values" = M[logicals])
}

#' Diagonal Part of a Matrix
#'
#' Similarly to `upper.tri` and `lower.tri`, returns a matrix of `logical`s to identify the diagonal of a square matrix
#'
#' @param M matrix or other R object with `length(dim(x)) == 2`. For back compatibility reasons, when the above is not fulfilled, `as.matrix(x)` is called first
#' @param output choose to either return a matrix of logical values (`TRUE`s on the diagonal, i.e. comparable behavior than `diag` for instance) or a vector of the subsetted values of the matrix
#'
#' @return square logical matrix with diagonal of `TRUE`s
#' @export
#'
#' @examples
#' M<- matrix(sample(1:10,16,replace = TRUE),4,4)
#' diagonal(M)
diagonal <- function(M,output=c("matrix.logical","vector.values")) {
  output <- match.arg(output)
  if(dim(M)[1]==dim(M)[2]) logicals<- upper.tri(M,diag = TRUE)&!upper.tri(M,diag = FALSE) else stop("Matrix provided is not a square matrix.")
  switch(output,
         "matrix.logical" = logicals,
         "vector.values" = M[logicals])
}

#' Reassign values in the non diagonal part of a matrix
#' similar use as its `diag(x)` and `diag(x) <- value` equivalent
#'
#' @param x matrix or other R object with `length(dim(x)) == 2`
#' @param value scalar or vector of values to replace the non diagonal part of `x` with. Filled in by increasing row then by increasing column,
#'
#' @return a matrix with similar dimensions as `x` where the non-diagonal parts have been replaced by the value in `value`
#' @export
#'
#' @examples
#' set.seed(42)
#' M <- matrix(sample(1:10,16,replace = TRUE),4,4)
#' M
#' non.diagonal(M) <- round(runif(12,0,1))
#' M
#' non.diagonal(M) <- 42
#' M
#' non.diagonal(M) <- 1:12
#' M
`non.diagonal<-` <- function(x,value) {
  dx <- dim(x)
  if (length(dx) != 2L)
    stop("only matrix non-diagonals can be replaced")
  len.i <- prod(dx) - min(dx)
  len.v <- length(value)
  if (len.v != 1L && len.v != len.i)
    stop("replacement non-diagonal has wrong length")
  x[non.diagonal(x)] <- value
  x
}

#' Identify coordinates of non null non diagonal elements of a matrix
#'
#' @param M a square matrix
#'
#' @return an array of row and col indices of the positive and non diagonal elements
#' @noRd
non.zero.non.diag<- function(M) {which(M>0&!diagonal(M),arr.ind = TRUE,useNames = TRUE)}
