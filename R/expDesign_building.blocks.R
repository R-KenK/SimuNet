# expDesign "building-block" functions --------------------------------------------------------
. <- NULL

# scanList manipulations ----

## sum ----

#' Sum list of binary scans into a weighted adjacency matrix
#' Written as a S3 method to be applied to `scanList` or `sLlist` (list of `scanList`) objects.
#'
#' @param scan.list a `scanList` or `sLlist` object. See objects returned by
#'   [`simunet()`][simunet()]
#' @param ... additional arguments to be passed.
#'
#' `sum_scans()` can use argument:
#' * `which`: character scalar, either:
#'   * `"auto"`: determine automatically which `scanList` to sum
#'   * `"theoretical"`: sum the `theoretical.scanList` (useful in the case of empirical
#'   `scanList`s)
#'   * `"raw"`: sum the `raw.scanList` (useful to see the impact of chosen `mode`)
#'
#' @return a `weightedAdj` object, or list of such, consisting mainly on a weighted adjacency matrix
#'   where each edge weight is equal to the sum of all binary edges. Inherits from the previous
#'   `scanList` class (theoretical or empirical, inheriting from `scanList`), and keeps track of the
#'   `scan.list`'s list of attributes `attrs`.
#'
#'   Also adds these attributes to `attrs`:
#'   * `summed.scanList`: the original `scanList` (3D array) that has been summed
#'   * `sampled`: an integer matrix representing how many time each edge has been sampled (i.e. was
#'   not `NA`). Determined via [`count_nonNA()`][count_nonNA()]
#'
#' @export
#'
#' @seealso [simunet()], [design_exp()], [perform_exp()].
#'
#' @examples
#' set.seed(42)
#' n <- 5L
#' samp.effort <- 241L
#'
#' # Adjacency matrix import
#' ## random directed adjacency matrix
#' Adj <- matrix(sample(1:samp.effort,n * n),nrow = 5,dimnames = list(letters[1:n],letters[1:n]))
#' diag(Adj) <- 0L
#' Adj
#'
#' # social network simulations
#' ## theoretical scans
#' sL <- simunet(Adj = Adj,samp.effort = samp.effort,mode = "directed",n.scans = 120L)
#' sL
#' sL |> sum_scans()
#'
#' ## group-scan sampling
#' sL |> perform_exp(design_sampling("group",.6)) |> sum_scans()
#'
#' ## comparing group-scan and focal sampling
#' sL |> perform_exp(design_sampling("group",.6),
#'                   design_sampling("focal","even")
#'                   ) |> sum_scans()
sum_scans <- function(scan.list,...) {
  UseMethod("sum_scans")
}


#' sum_scans method for `scanList` objects
#' @export
#' @noRd
sum_scans.scanList <- function(scan.list,which = c("auto","theoretical","raw"),...) {
  which <- match.arg(which)
  sf <- attrs(scan.list,"Adj.subfun")
  sL.ori <- scan.list
  switch(which,
         "auto" = {},
         "theoretical" = scan.list <- attrs(scan.list,"theoretical.scanList"),
         "raw" = scan.list <- attrs(scan.list,"raw.scanList.type")
  )
  summed <- rowSums(scan.list,na.rm = TRUE,dims = 2L)
  summed <- copy_attrs_to(sL.ori,summed)
  attrs(summed,"summed.scanList") <- without_attrs(sL.ori)
  attrs(summed,"sampled") <- count_nonNA(scan.list)
  class(summed) <- c("weightedAdj",class(scan.list))
  summed
}

#' sum_scans method for `sLlist` objects
#' @export
#' @noRd
sum_scans.sLlist <- function(scan.list,which = c("auto","theoretical","raw"),...) {
  sLlapply(scan.list,sum_scans,which = which)
}

#' sum_scans method for `empirical` (`scanList`) objects
#' @export
#' @noRd
sum_scans.empirical <- function(scan.list,which = c("auto","theoretical","raw"),...) {
  NextMethod()
}

## scale ----

#' Scale list of binary scans into a weighted adjacency matrix
#' Scaling here is dividing by the sum of 1s by the number of time an edge has been observed
#' (whether it was 0 or 1).
#'
#' Written as a S3 method to be applied to `scanList` or `sLlist` (list of `scanList`) objects.
#'
#' @param scan.list a `scanList` or `sLlist` object. See objects returned by
#'   [`simunet()`][simunet()]
#' @param ... additional arguments to be passed. At the moment `scale_scans()` does not use
#'
#' At the moment `scale_scans()` does not use additional argument, arguments passed will be ignored.
#'
#' @return a `scaled` object, or list of such, consisting mainly on a weighted adjacency matrix
#'   where each edge weight is equal to the sum of all binary edges divided by the number of times
#'   they have been sampled (determined via [`count_nonNA()`][count_nonNA()]). Inherits from
#'   `weightedAdj` and the previous `scanList` class (theoretical or empirical, inheriting from
#'   `scanList`), and keeps track of the `scan.list`'s list of attributes `attrs`
#'
#' @export
#'
#' @seealso [simunet()], [design_exp()], [perform_exp()], [sum_scans()].
#'
#' @examples
#' set.seed(42)
#' n <- 5L
#' samp.effort <- 241L
#'
#' # Adjacency matrix import
#' ## random directed adjacency matrix
#' Adj <- sample(1:samp.effort,n * n) |>
#'   matrix(nrow = 5,dimnames = list(letters[1:n],letters[1:n]))
#' diag(Adj) <- 0L
#' Adj
#'
#' # social network simulations
#' ## theoretical scans
#' sL <- simunet(Adj = Adj,samp.effort = samp.effort,mode = "directed",n.scans = 120L)
#' sL
#' # scale_scans() can scale `weightedAdj` objects...
#' sL |> sum_scans() |> scale_scans()
#'
#'
#' # ... or `scanList` object directly
#' ## group-scan sampling
#' sL |> perform_exp(design_sampling("group",.6)) |> scale_scans()
#'
#' ## comparing group-scan and focal sampling
#' sL |> perform_exp(design_sampling("group",.6),
#'                   design_sampling("focal","even")
#'                   ) |> scale_scans()
scale_scans <- function(scan.list,...) {
  UseMethod("scale_scans")
}

#' scale_scans method for `weightedAdj` objects
#' @export
#' @noRd
scale_scans.weightedAdj <- function(scan.list,...) {
  sf <- attrs(scan.list,"Adj.subfun")
  sampled <- attrs(scan.list,"sampled")
  scaled <- scan.list
  scaled[sf(scaled)] <- scan.list[sf(scan.list)] / sampled[sf(sampled)]
  ifelse(!is.infinite(scaled),sampled,NA)
  scaled <- copy_attrs_to(scan.list,scaled)
  class(scaled) <- c("scaled",class(scan.list))
  scaled
}

#' scale_scans method for `scanList` objects
#' @export
#' @noRd
scale_scans.scanList <- function(scan.list,...) {
  scale_scans(sum_scans(scan.list))
}

#' scale_scans method for `empirical` (`scanList`) objects
#' @export
#' @noRd
scale_scans.empirical <- function(scan.list,...) {
  NextMethod()
}

#' scale_scans method for `sLlist` objects
#' @export
#' @noRd
scale_scans.sLlist <- function(scan.list,...) {
  sLlapply(scan.list,scale_scans)
}

## count NAs ----
#' Count number of times edges have not been observed
#' For empirical `scanList`, an edge that is not observed during a scan is `NA`. This function
#' counts how many time this was the case in the inputted `scan.list` for each edge.
#'
#' Written as a S3 method to be applied to `scanList` or `sLlist` (list of `scanList`) objects.
#'
#' @param scan.list a `scanList` or `sLlist` object, where an unobserved edge (whether it is 0 or 1)
#'   is `NA`. See objects returned by [`simunet()`][simunet()]
#' @param ... additional arguments to be passed. At the moment `scale_scans()` does not use
#'
#' At the moment `count_NA()` does not use additional argument, arguments passed will be ignored.
#'
#' @return an integer matrix, or list of such, representing how many time each edge has been
#'   unobserved (i.e. was `NA`). Inherits from `weightedAdj` and the previous `scanList` class (theoretical
#'   or empirical, inheriting from `scanList`), and keeps track of the `scan.list`'s list of
#'   attributes `attrs`.
#' @export
#'
#' @seealso [simunet()], [design_exp()], [perform_exp()], [count_nonNA()].
#'
#' @examples
#' set.seed(42)
#' n <- 5L
#' samp.effort <- 241L
#'
#' # Adjacency matrix import
#' ## random directed adjacency matrix
#' Adj <- sample(1:samp.effort,n * n) |>
#'   matrix(nrow = 5,dimnames = list(letters[1:n],letters[1:n]))
#' diag(Adj) <- 0L
#' Adj
#'
#' # social network simulations
#' ## theoretical scans
#' sL <- simunet(Adj = Adj,samp.effort = samp.effort,mode = "directed",n.scans = 120L)
#' sL
#'
#' ## group-scan sampling
#' sL |> perform_exp(design_sampling("group",.6)) |> count_NA()
count_NA <- function(scan.list,...) {
  UseMethod("count_NA")
}

#' count_NA method for `scanList` objects
#' @export
#' @noRd
count_NA.scanList <- function(scan.list,...) {
  sf <- attrs(scan.list,"Adj.subfun")
  scan.sampled <-  copy_attrs_to(from = scan.list,ifelse(is.na(scan.list),1L,0L))
  scan.sampled <- rowSums(scan.sampled,na.rm = TRUE,dims = 2L)
  scan.sampled[!sf(scan.sampled)] <- 0L
  scan.sampled <- copy_attrs_to(scan.list,scan.sampled)
  class(scan.sampled) <- c("weightedAdj",class(scan.list))
  scan.sampled
}

#' count_NA method for `empirical` (`scanList`) objects
#' @export
#' @noRd
count_NA.empirical <- function(scan.list,...) {
  NextMethod()
}

#' count_NA method for `sLlist` objects
#' @export
#' @noRd
count_NA.sLlist <- function(scan.list,...) {
  sLlapply(scan.list,count_NA)
}

## count non NAs ----

#' Count number of times edges have not been observed
#' For empirical `scanList`, an edge that is observed during a scan is *not* `NA` (it will be either
#' 0 or 1). This function counts how many time this was the case in the inputted `scan.list` for
#' each edge.
#'
#' Written as a S3 method to be applied to `scanList` or `sLlist` (list of `scanList`) objects.
#'
#' @param scan.list a `scanList` or `sLlist` object, where an unobserved edge (whether it is 0 or 1)
#'   is `NA`. See objects returned by [`simunet()`][simunet()]
#' @param ... additional arguments to be passed. At the moment `scale_scans()` does not use
#'
#' At the moment `count_nonNA()` does not use additional argument, arguments passed will be ignored.
#'
#' @return an integer matrix, or list of such, representing how many time each edge has been sampled
#'   (i.e. was *not* `NA`). Inherits from `weightedAdj` and the previous `scanList` class (theoretical or
#'   empirical, inheriting from `scanList`), and keeps track of the `scan.list`'s list of attributes
#'   `attrs`.
#' @export
#'
#' @seealso [simunet()], [design_exp()], [perform_exp()]], [count_NA()].
#'
#' @examples
#' set.seed(42)
#' n <- 5L
#' samp.effort <- 241L
#'
#' # Adjacency matrix import
#' ## random directed adjacency matrix
#' Adj <- sample(1:samp.effort,n * n) |>
#'   matrix(nrow = 5,dimnames = list(letters[1:n],letters[1:n]))
#' diag(Adj) <- 0L
#' Adj
#'
#' # social network simulations
#' ## theoretical scans
#' sL <- simunet(Adj = Adj,samp.effort = samp.effort,mode = "directed",n.scans = 120L)
#' sL
#'
#' ## group-scan sampling
#' sL |> perform_exp(design_sampling("group",.6)) |> count_nonNA()
count_nonNA <- function(scan.list,...) {
  UseMethod("count_nonNA")
}

#' count_nonNA method for `scanList` objects
#' @export
#' @noRd
count_nonNA.scanList <- function(scan.list,...) {
  sf <- attrs(scan.list,"Adj.subfun")
  scan.sampled <- copy_attrs_to(from = scan.list,to = ifelse(is.na(scan.list),0L,1L))
  scan.sampled <- rowSums(scan.sampled,na.rm = TRUE,dims = 2L)
  scan.sampled[!sf(scan.sampled)] <- 0L
  scan.sampled <- copy_attrs_to(scan.list,scan.sampled)
  class(scan.sampled) <- c("weightedAdj",class(scan.list))
  scan.sampled
}

#' count_nonNA method for `empirical` (`scanList`) objects
#' @export
#' @noRd
count_nonNA.empirical <- function(scan.list,...) {
  NextMethod()
}

#' count_nonNA method for `sLlist` objects
#' @export
#' @noRd
count_nonNA.sLlist <- function(scan.list,...) {
  sLlapply(scan.list,count_nonNA)
}

## add scans ----

#' Perform additional scans and add them to the `scanList`
#' New scans rely on edge probability matrix previously drawn from beta distributions. This can be
#' used to compare sampling regime when it is expected that one will sample edges less in a
#' predictable fashion (e.g. group-scan vs focal sampling)
#'
#' Written as a S3 method to be applied to `scanList` or `sLlist` (list of `scanList`) objects.
#'
#' @param scan.list a `scanList` or `sLlist` object. See objects returned by
#'   [`simunet()`][simunet()]
#' @param ... additional arguments to be passed. At the moment `scale_scans()` does not use
#'
#' `add_scans()` can use argument:
#' * `new.scans`: integer scalar, the number of additional new scans to be performed and added to
#' the `scan.list`
#'
#' @return a `scanList`, or list of such, to which additional scans have been added to the 3D array. See also [`simunet()`][simunet()].
#'
#'   Also modifies these attributes to `attrs`:
#'   * `n.scans`: integer indicating how many scans have been performed in total. `n.scans` itself
#'   has the attribute `scans.performed` that keeps track of what "batches" of scans have been
#'   performed, by growing a vector of "scan batch" sizes
#'
#' @export
#'
#' @seealso [simunet()], [design_exp()], [perform_exp()].
#'
#' @examples
#' set.seed(42)
#' n <- 5L
#' samp.effort <- 241L
#'
#' # Adjacency matrix import
#' ## random directed adjacency matrix
#' Adj <- sample(1:samp.effort,n * n) |>
#'   matrix(nrow = 5,dimnames = list(letters[1:n],letters[1:n]))
#' diag(Adj) <- 0L
#' Adj
#'
#' # social network simulations
#' ## theoretical scans
#' sL <- simunet(Adj = Adj,samp.effort = samp.effort,mode = "directed",n.scans = 120L)
#' sL
#'
#' ## focal-scan sampling
#' sL |> perform_exp(design_sampling("focal","even")) |> add_scans(50)
add_scans <- function(scan.list,...) {
  UseMethod("add_scans")
}

#' add_scans method for `scanList` objects
#' @export
#' @noRd
add_scans.scanList <- function(scan.list,new.scans,...) {
  edge.Prob <- reconstruct_edgeProb(scan.list)
  n.scans <- attrs(scan.list,"n.scans")

  new.scan.list <- generate_scanList(edge.Prob = edge.Prob,n.scans = new.scans)

  new.scan.list <- rbind(scan.list,new.scan.list)
  new.scan.list <- copy_attrs_to(scan.list,new.scan.list)

  total.scans <- attrs(new.scan.list,"n.scans") + new.scans
  attr(total.scans,"scans.performed") <-
    c(new.scans,
      if (is.null(attr(n.scans,"scans.performed"))) n.scans
      else attr(n.scans,"scans.performed")
    )
  attrs(new.scan.list,"n.scans") <- total.scans
  new.scan.list
}

#' add_scans method for `empirical` (`scanList`) objects
#' @export
#' @noRd
add_scans.empirical <- function(scan.list,new.scans,...) {
  NextMethod()
}

#' add_scans method for `sLlist` objects
#' @export
#' @noRd
add_scans.sLlist <- function(scan.list,new.scans,...) {
  sLlapply(scan.list,add_scans,new.scans = new.scans)
}

## remove peripheral individual ----

#' Remove from all scans the (overall) most peripheral individual
#' Individual centrality based on eigen vectors. Mostly given as an example of experimental
#' manipulations that could be performed on `scanList` as `expDesign`, even as user-defined
#' functions
#'
#' Written as a S3 method to be applied to `scanList` or `sLlist` (list of `scanList`) objects.
#'
#' @param scan.list a `scanList` or `sLlist` object. See objects returned by
#'   [`simunet()`][simunet()]
#' @param ... additional arguments to be passed. At the moment `scale_scans()` does not use
#'
#' At the moment `remove_mostPeripheral()` does not use additional argument, arguments passed will be ignored.
#'
#' @return a `scanList`, or list of such, in which the most peripheral node has been removed. See
#'   also [`simunet()`][simunet()].
#'
#' @seealso [simunet()], [design_exp()], [perform_exp()].
#'
#' @export
#'
#' @examples
#' set.seed(42)
#' n <- 5L
#' samp.effort <- 241L
#'
#' # Adjacency matrix import
#' ## random directed adjacency matrix
#' Adj <- sample(1:samp.effort,n * n) |>
#'   matrix(nrow = 5,dimnames = list(letters[1:n],letters[1:n]))
#' Adj[lower.tri(Adj,diag = FALSE)] <- 0L
#' Adj
#'
#' # social network simulations
#' ## theoretical scans
#' sL <- simunet(Adj = Adj,samp.effort = samp.effort,mode = "upper",n.scans = 120L)
#' sL
#'
#' ## focal-scan sampling
#' sL |> perform_exp(design_sampling("focal","even")) |> remove_mostPeripheral()
remove_mostPeripheral <- function(scan.list,...) {
  UseMethod("remove_mostPeripheral")
}

#' remove_mostPeripheral method for `scanList` objects
#' @export
#' @noRd
remove_mostPeripheral.scanList <- function(scan.list,...) {
  mode <- attrs(scan.list,"mode")
  directed <- switch(mode,"directed" = TRUE,FALSE)
  G <- igraph::graph.adjacency(sum_scans(scan.list),weighted = TRUE)
  m <- which.min(igraph::eigen_centrality(G,directed = directed)$vector)
  copy_attrs_to(from = scan.list,to = scan.list[-c(m),-c(m),])
}

#' remove_mostPeripheral method for `sLlist` objects
#' @export
#' @noRd
remove_mostPeripheral.sLlist <- function(scan.list,...) {
  sLlapply(remove_mostPeripheral,scan.list)
}

## remove central individual ----

#' Remove from all scans the (overall) most central individual
#' Individual centrality based on eigen vectors. Mostly given as an example of experimental
#' manipulations that could be performed on `scanList` as `expDesign`, even as user-defined
#' functions
#'
#' Written as a S3 method to be applied to `scanList` or `sLlist` (list of `scanList`) objects.
#'
#' @param scan.list a `scanList` or `sLlist` object. See objects returned by
#'   [`simunet()`][simunet()]
#' @param ... additional arguments to be passed. At the moment `scale_scans()` does not use
#'
#' At the moment `remove_mostCentral()` does not use additional argument, arguments passed will be
#' ignored.
#'
#' @return a `scanList`, or list of such, in which the most central node has been removed. See
#'   also [`simunet()`][simunet()].
#'
#' @seealso [simunet()], [design_exp()], [perform_exp()].
#'
#' @export
#'
#' @examples
#' set.seed(42)
#' n <- 5L
#' samp.effort <- 241L
#'
#' # Adjacency matrix import
#' ## random directed adjacency matrix
#' Adj <- sample(1:samp.effort,n * n) |>
#'   matrix(nrow = 5,dimnames = list(letters[1:n],letters[1:n]))
#' Adj[lower.tri(Adj,diag = FALSE)] <- 0L
#' Adj
#'
#' # social network simulations
#' ## theoretical scans
#' sL <- simunet(Adj = Adj,samp.effort = samp.effort,mode = "upper",n.scans = 120L)
#' sL
#'
#' ## focal-scan sampling
#' sL |> perform_exp(design_sampling("focal","even")) |> remove_mostCentral()
remove_mostCentral <- function(scan.list,...) {
  UseMethod("remove_mostCentral")
}

#' remove_mostCentral method for `scanList` objects
#' @export
#' @noRd
remove_mostCentral.scanList <- function(scan.list,...) {
  mode <- attrs(scan.list,"mode")
  directed <- switch(mode,"directed" = TRUE,FALSE)
  G <- igraph::graph.adjacency(sum_scans(scan.list),weighted = TRUE)
  m <- which.max(igraph::eigen_centrality(G,directed = directed)$vector)
  copy_attrs_to(from = scan.list,to = scan.list[-c(m),-c(m),])
}

#' remove_mostCentral method for `sLlist` objects
#' @export
#' @noRd
remove_mostCentral.sLlist <- function(scan.list,...) {
  sLlapply(remove_mostCentral,scan.list)
}

