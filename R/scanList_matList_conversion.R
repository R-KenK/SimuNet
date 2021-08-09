# scanList convenience and compatibility functions ----

#' Convert scanList 3D array into list of 2D matrices
#'
#' Convenience function to allow the use of functions designed to affect list of matrices instead of
#' 3D array (specifically `scanList` objects) where the 3rd dimension is the scan index. Keeps track
#' of the attributes list `attrs`
#'
#' To back-transform a list of matrices into a `scanList` object, users can use
#' `matList2array()`, wrapper to [`simplify2array()`][simplify2array()] that
#' handles attributes list `attrs`
#'
#' @param scan.list a 3 dimensional array representing adjacency matrices (first 2 dimensions)
#'   throughout the different scans (3rd dimension)
#'
#' @return a list of matrices where each element is the 2D matrix for each scan index (3rd dimension of `array.3D`)
#' @export
#'
#' @examples
#' set.seed(42)
#' n <- 5L
#' samp.effort <- 100L
#'
#' # Adjacency matrix import
#' ## random directed adjacency matrix
#' Adj <- sample(1:samp.effort,n * n) |>
#'   matrix(nrow = 5,dimnames = list(letters[1:n],letters[1:n]))
#' Adj[lower.tri(Adj,diag = TRUE)] <- 0L
#' Adj
#' mL <- simunet(Adj,samp.effort,"upper",10) |> scanList2matList()
#' mL
#' mL |> matList2scanList()
scanList2matList <- function(scan.list) {
  mat.list <- array2matList(scan.list)
  mat.list <- copy_attrs_to(scan.list,mat.list)
  # replaces the scanList class by matList, but keeps the rest
  class(mat.list)[length(class(mat.list))] <- "matList"
  mat.list
}

#' Convert `matList` into `scanList` objects
#'
#' Convenience function to back-transform a `matList` object (list of matrices) into a `scanList`
#' object
#'
#' @param mat.list a 3 dimensional array representing adjacency matrices (first 2 dimensions)
#'   throughout the different scans (3rd dimension)
#'
#' @return scanList object. See [`simunet()`][simunet()]
#' @export
#'
#' @examples
#' set.seed(42)
#' n <- 5L
#' samp.effort <- 100L
#'
#' # Adjacency matrix import
#' ## random directed adjacency matrix
#' Adj <- sample(1:samp.effort,n * n) |>
#'   matrix(nrow = 5,dimnames = list(letters[1:n],letters[1:n]))
#' Adj[lower.tri(Adj,diag = TRUE)] <- 0L
#' Adj
#' mL <- simunet(Adj,samp.effort,"upper",10) |> scanList2matList()
#' mL |> matList2scanList()
matList2scanList <- function(mat.list) {
  scan.list <- matList2array(mat.list)
  scan.list <- copy_attrs_to(mat.list,scan.list)
  # replaces the scanList class by matList, but keeps the rest
  class(scan.list)[length(class(scan.list))] <- "scanList"
  scan.list
}

#' Convert scanList 3D array into list of 2D matrices
#'
#' Convenience function to allow the use of functions designed to affect list of matrices instead of
#' 3D array (specifically `scanList` objects) where the 3rd dimension is the scan index. Keeps track
#' of the attributes list `attrs`
#'
#' To back-transform a list of matrices into a `scanList` object, users can use
#' `matList2array()`, wrapper to [`simplify2array()`][simplify2array()] that
#' handles attributes list `attrs`
#'
#' @param array.3D a 3 dimensional array representing adjacency matrices (first 2 dimensions)
#'   throughout the different scans (3rd dimension)
#'
#' @return a list of matrices where each element is the 2D matrix for each scan index (3rd dimension of `array.3D`)
#' @export
#'
#' @keywords internal
array2matList <- function(array.3D) {
  lapply(1:dim(array.3D)[3],\(s) array.3D[,,s])
}

#' Convert list of 2D matrices into scanList 3D array
#'
#' Convenience function to back-transform a list of matrices into a 3D array
#'
#' @param mat.list a 3 dimensional array representing adjacency matrices (first 2 dimensions)
#'   throughout the different scans (3rd dimension)
#'
#' @return a 3 dimensional array representing adjacency matrices (first 2 dimensions) throughout the
#'   different scans (3rd dimension)
#' @export
#'
#' @keywords internal
matList2array <- function(mat.list) {
  simplify2array(mat.list)
}

#' Print method for `matList` objects
#' @export
#' @noRd
print.matList <- function(x,...) {
  print.default(without_attrs(x))
  cat("\n\nHidden attributes:",names(get_attrs(x)))
}

