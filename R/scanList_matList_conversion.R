# scanList convenience and compatibility functions ----

## To and from matList ----

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
#'
#' sL <- simunet(Adj,samp.effort,"upper",10)
#' mL <- sL |> scanList2matList()
#' mL
#' mL |> matList2scanList()
#' identical(mL |> matList2scanList(),sL)
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
#'
#' sL <- simunet(Adj,samp.effort,"upper",10)
#' mL <- sL |> scanList2matList()
#' mL
#' mL |> matList2scanList()
#' identical(mL |> matList2scanList(),sL)
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
#' @return a list of matrices where each element is the 2D matrix for each scan index (3rd dimension
#'   of `array.3D`)
#' @export
#'
#' @keywords internal
array2matList <- function(array.3D) {
  lapply(1:dim(array.3D)[3],function(s) array.3D[,,s])
}

#' Convert list of 2D matrices into scanList 3D array
#'
#' Convenience function to back-transform a list of matrices into a 3D array
#'
#' @param mat.list a list of matrices where each element is the 2D matrix for each scan index (3rd
#'   dimension of `array.3D`)
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
  format_attributes(x,...)
  invisible(x)
}

## To and from igraph ----

#' Convert scanList 3D array into `igraph` object (with `attrs`)
#'
#'
#' To back-transform an `igraphSN` object into a `scanList` object, users can use
#' `igraph2scanList()`, wrapper to [`igraph2array()`][igraph2array()] that
#' handles attributes list `attrs`
#'
#' @param scan.list a 3 dimensional array representing adjacency matrices (first 2 dimensions)
#'   throughout the different scans (3rd dimension)
#'
#' @return an `igraphSN` object:
#'
#' * the `igraph` network object obtained from the weighted adjacency matrix corresponding to the
#' inputted `scanList`
#' * the `attrs` attributes list carried over from the inputted `scanList`
#'
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
#'
#' sL <- simunet(Adj,samp.effort,"upper",10)
#' G <- sL |> scanList2igraph()
#' G
#' G |> igraph2scanList()
#' G |> igraph2scanList(format = "scanList")
#' identical(G |> igraph2scanList(format = "scanList"),sL)
scanList2igraph <- function(scan.list) {
  if (!inherits(scan.list,"weightedAdj")) {
    scan.list <- sum_scans(scan.list)
  }
  G <- convert_to_igraph(scan.list)
  G <- copy_attrs_to(scan.list,G,copy.class = FALSE)
  attrs(G,"original.class") <- class(scan.list)
  class(G) <- c("igraphSN",class(G))
  G
}

#' Convert scanList 3D array into `igraph` object
#'
#' @param scan.list a 3 dimensional array representing adjacency matrices (first 2 dimensions)
#'   throughout the different scans (3rd dimension)
#'
#' @return an `igraph` object
#' @export
#'
#' @importFrom igraph graph.adjacency
#'
#' @keywords internal
convert_to_igraph <- function(scan.list) {
  mode <- scan.list$mode
  igraph::graph.adjacency(scan.list,mode = mode,weighted = TRUE)
}

#' Convert `igraphSN` into `scanList` objects
#'
#' Convenience function to back-transform a `igraphSN` object (`igraph` object with `attrs`) into a
#' `scanList` object
#'
#' @param G an `igraphSN` object (see [`scanList2igraph()`][scanList2igraph()])
#' @param format character, option to specify if a `weightedAdj` or a 3D binary `scanList` should be
#'   returned
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
#'
#' sL <- simunet(Adj,samp.effort,"upper",10)
#' G <- sL |> scanList2igraph()
#' G
#' G |> igraph2scanList()
#' G |> igraph2scanList(format = "scanList")
#' identical(G |> igraph2scanList(format = "scanList"),sL)
igraph2scanList <- function(G,format = c("weightedAdj","scanList")) {
  format <- match.arg(format)

  scan.list <- igraph2array(G)
  scan.list <- copy_attrs_to(G,scan.list,copy.class = FALSE)
  class(scan.list) <- attrs(G,"original.class")
  attrs(scan.list,"original.class") <- NULL

  switch(format,
         "weightedAdj" = {
           if (!inherits(scan.list,"weightedAdj"))
             sum_scans(scan.list)
           else
             scan.list
         },
         "scanList" = {
           if (inherits(scan.list,"weightedAdj"))
             weightedAdj2scanList(scan.list)
           else
             scan.list
         }
  )
}

#' Convert list of 2D matrices into scanList 3D array
#'
#' Convenience function to back-transform a list of matrices into a 3D array
#'
#' @param G a list of matrices where each element is the 2D matrix for each scan index (3rd
#'   dimension of `array.3D`)
#'
#' @return a 3 dimensional array representing adjacency matrices (first 2 dimensions) throughout the
#'   different scans (3rd dimension)
#' @export
#'
#' @importFrom igraph get.adjacency
#'
#' @keywords internal
igraph2array <- function(G) {
  type <- switch(attrs(G,"mode"),"upper" = "upper","lower" = "lower","both")
  igraph::get.adjacency(G,type = type,attr = "weight",sparse = FALSE)
}

#' Print method for `igraphSN` objects
#' @export
#' @noRd
print.igraphSN <- function(x,...) {
  NextMethod()
  format_attributes(x,...)
  invisible(x)
}


## From weightedAdj to binary scanList ----

#' Convert `weightedAdj` into 3D binary array `scanList` objects
#' Back-transform `weightedAdj`into the `scanList` it was summed from (see [`sum_scans()`][sum_scans()])
#'
#' @param sum a `weightedAdj` object
#'
#' @return `scanList` object. See [`simunet()`][simunet()]
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
#'
#' sL <- simunet(Adj,samp.effort,"upper",10)
#' wAdj <- sL |> sum_scans()
#' wAdj
#' wAdj |> weightedAdj2scanList()
#' identical(wAdj |> weightedAdj2scanList(),sL)
weightedAdj2scanList <- function(sum) {
  scan.list <- sum$summed.scanList
  scan.list <- copy_attrs_to(sum,scan.list,copy.class = FALSE)
  class(scan.list) <- class(sum)[-1]
  attrs(scan.list,"summed.scanList") <- attrs(scan.list,"sampled") <- NULL
  scan.list
}
