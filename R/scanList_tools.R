#' Generator for `scanList` objects
#' Internal use. The user should rather rely on `simunet` as a wrapper for the different steps
#' needed to perform the simulations from inputted data.
#'
#' @param edge.Prob an `edgeProb` object, i.e. a list containing:
#' * `P`: the edge presence probability matrix
#' * `Adj`: the inputted `Adj`
#' * `samp.effort`: the inputted `samp.effort`
#' * `mode`: the inputted `mode`
#' * `Adj.subfun`: the inputted `Adj.subfun`
#' @param samp.effort integer scalar, the sampling effort, or number of scans, that led to obtaining
#'   of `Adj`
#'
#' @return a `theoretical` inheriting from `scanList` object, primarily a 3 dimensional array
#'   representing the (binary) adjacency matrices (coded within the first two dimensions of the
#'   3D-array) obtained at each simulated scan (coded as the 3rd dimension of the 3D-array), and a
#'   list of attributes, `attrs`.
#'
#'   The list of attributes `attrs` contains:
#'  * `scanList.type`: character scalar, `"theoretical"` at first and `"empirical"` after a non-`NULL`
#'  experimental manipulation has been applied to the `scanList` (via [`perform_exp()`][perform_exp()]
#'  and a `expDesign` object)
#'  * `raw.scanList`: the 3D binary array, directed, before potential symmetrization attempt by
#'  applying the igraph's mode via [`apply_mode()`][apply_mode()]
#'  * `Adj`: integer matrix, `Adj` contained in `edge.Prob`
#'  * `samp.effort`: integer, `samp.effort` contained in `edge.Prob`
#'  * `n.scans`: inputted `n.scans`
#'  * `mode`: character scalar, `mode` contained in `edge.Prob`
#'  * `Adj.subfun`: function, `Adj.subfun` contained in `edge.Prob`
#'  * `edge.Prob`: numeric matrix,`edge.Prob$P` (only the probability matrix) data contained in
#'  `edge.Prob`
#'
#' @seealso [simunet()], [generate_edgeProb()], [draw_edgeProb()], [generate_empiscanList()].
#'
#' @export
#'
#' @keywords internal
generate_scanList <- function(edge.Prob,n.scans){
  raw.scanList <- draw_raw_scanList(edge.Prob = edge.Prob,n.scans = n.scans)
  scanList <- apply_mode(raw.scanList = raw.scanList,mode = edge.Prob$mode)
  attr(scanList,"attrs") <-
    list(
      scanList.type = "theoretical",
      raw.scanList = raw.scanList,
      Adj = edge.Prob$Adj,
      samp.effort = edge.Prob$samp.effort,
      n.scans = n.scans,
      mode = edge.Prob$mode,
      Adj.subfun = edge.Prob$Adj.subfun,
      edge.Prob = edge.Prob$P # in here it is not a edgeProb object anymore, to avoid storing redundant variables,
    )
  class(scanList) <- c("theoretical","scanList")
  scanList
}

#' Generator for *empirical* `scanList` objects
#' Internal use. The user should rather rely on [`simunet()`][simunet()] and/or
#' [`perform_exp()`][perform_exp()] as a wrapper for the different steps needed to perform
#' simulations and experimental manipulations.
#'
#' @param scan.list a `scanList` object (see [`simunet()`][simunet()])
#' @param exp.design an `expDesign` object. See objects returned by [`design_exp()`][design_exp()]
#'
#' @return an `empirical` inheriting from `scanList` object, primarily a 3 dimensional array
#'   representing the (binary) adjacency matrices (coded within the first two dimensions of the
#'   3D-array) obtained at each simulated scan (coded as the 3rd dimension of the 3D-array), and a
#'   list of attributes, `attrs`.
#'
#'   The list of attributes `attrs` contains:
#'   * all the previous attributes contained in `scan.list`'s `attrs` attributes list, as well as:
#'     * `scanList.type`: character scalar, changed from `"theoretical"` to `"empirical"`
#'     * `theoretical.scanList`: the 3D array _before_ the experimental manipulations contained in
#'     the inputted `exp.design` have been applied
#'
#' @export
#'
#' @seealso [simunet()], [generate_edgeProb()], [draw_edgeProb()], [generate_scanList()].
#'
#' @keywords internal
generate_empiscanList <- function(scan.list,exp.design) {
  empiscanList <- exp.design$FUN.seq(scan.list)
  attrs(empiscanList,"scanList.type") <- "empirical"
  attrs(empiscanList,"theoretical.scanList") <- without_attrs(scan.list)
  class(empiscanList)<- c("empirical","scanList")
  empiscanList
}

#'  Draw edge presence according to the edge presence probability matrix
#'
#' @param edge.Prob an `edgeProb` object, i.e. a list containing:
#' * `P`: the edge presence probability matrix
#' * `Adj`: the inputted `Adj`
#' * `samp.effort`: the inputted `samp.effort`
#' * `mode`: the inputted `mode`
#' * `Adj.subfun`: the inputted `Adj.subfun`
#' @param n.scans integer scalar, number of scans to generate in the simulation
#'
#' @return a 3 dimensional array
#'   representing the (binary) adjacency matrices (coded within the first two dimensions of the
#'   3D-array) obtained at each simulated scan (coded as the 3rd dimension of the 3D-array)
#'
#' @seealso [simunet()], [generate_edgeProb()], [draw_edgeProb()], [generate_scanList()].
#'
#' @keywords internal
draw_raw_scanList <- function(edge.Prob,n.scans) {
  sL <- vapply(
    1:n.scans,
    \(s) {
      stats::rbinom(edge.Prob$P,1L,edge.Prob$P)
    },edge.Prob$Adj
  )
  class(sL) <- "scanList"
  sL
}

# scanList tools ------------------------------------------------------------------------------

#' `scanList`'s `attrs` attributes related convenience functions: retrieve all attributes
#'
#' @param scan.list a `scanList` object (see [`simunet()`][simunet()])
#'
#' @return list, the list of attributes stored in `scan.list`'s `attrs` attribute
#' @export
#'
#' @keywords internal
get_attrs <- function(scan.list) {
  attr(scan.list,"attrs")
}

#' `scanList`'s `attrs` attributes related convenience functions: output the 3D array only
#'
#' @param scan.list a `scanList` object (see [`simunet()`][simunet()])
#'
#' @return the 3D array without its `attrs` argument
#' @noRd
without_attrs <- function(scan.list) {
  attr(scan.list,"attrs") <- NULL
  scan.list
}

#' `scanList`'s `attrs` attributes related convenience functions: retrieve or modify attributes
#' `attrs()` and `attrs()<-` can be used to retrieve the named attributes contained in the
#' attributes list `attrs` of a `scanList` object
#'
#' @param scan.list a `scanList` object (see [`simunet()`][simunet()])
#' @param a character (scalar or vector), the name(s) of the attribute(s) to retrieve, modify or add
#'
#' @return  the attribute(s) requested, or the `scan.list` which `attrs` attribute has been modified
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
#' sL <- simunet(Adj = Adj,samp.effort = samp.effort,mode = "upper",n.scans = 120L)
#'
#' # retrieve all attributes in `attrs`
#' sL |> attrs()
#'
#' # retrieve a specific attribute from `attrs`
#' sL |> attrs("edge.Prob")
#'
#' # modify a specific attribute from `attrs` (internal use)
#' attrs(sL,"scanList.type") <- "empirical"
#' attrs(sL,"scanList.type")
attrs <- function(scan.list,a = NULL) {
  if (is.null(a)) return(get_attrs(scan.list))
  get_attrs(scan.list)[[a]]
}

#' `scanList`'s `attrs` attributes related convenience functions: retrieve or modify attributes
#' `attrs()` and `attrs()<-` can be used to retrieve the named attributes contained in the
#' attributes list `attrs` of a `scanList` object
#'
#' @param x a `scanList` object (see [`simunet()`][simunet()])
#' @param which character (scalar or vector), the name(s) of the attribute(s) to retrieve, modify or
#'   add
#' @param value object to replace the requested attribute with
#'
#' @return the `scan.list` which `attrs` attribute has been modified
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
#' sL <- simunet(Adj = Adj,samp.effort = samp.effort,mode = "upper",n.scans = 120L)
#'
#' # retrieve all attributes in `attrs`
#' sL |> attrs()
#'
#' # retrieve a specific attribute from `attrs`
#' sL |> attrs("edge.Prob")
#'
#' # modify a specific attribute from `attrs` (internal use)
#' attrs(sL,"scanList.type") <- "empirical"
#' attrs(sL,"scanList.type")
`attrs<-` <- function(x,which,value) {
  new <- get_attrs(x)
  new[[which]] <- value
  attr(x,"attrs") <- new
  x
}

#' `scanList`'s `attrs` attributes related convenience functions: copy attrs from one `scanList` to
#' another
#'
#' @param from a `scanList` object which `attrs` attribute to copy (see [`simunet()`][simunet()])
#' @param to a `scanList` object to which `attrs` attribute should be pasted (see
#'   [`simunet()`][simunet()])
#'
#' @return a `scanList` object, the 3D array containted in `to` with `from`'s `attrs`
#' @export
#'
#' @keywords internal
copy_attrs_to <- function(from,to) {
  if (!inherits(to,"scanList")) {class(to) <- class(from)}
  attr(to,"attrs") <- attrs(from)
  to
}

#' Shortcut to a `lapply` equivalent to apply a function to each 2D matrix contained in a `scanList`
#' Written analogously to [vapply()]. Values returned by `.f` should be a similarly dimensionned
#' matrix as the first one contained in the 3D array
#'
#' @param sL a `scanList` object (see [`simunet()`][simunet()])
#' @param .f a function,to apply a function to each 2D matrix contained in `sL`
#' @param ... extra argument to be passed, notably named arguments used by `.f` (see [lapply()])
#' @param USE.NAMES logical; if `TRUE` and if `X` is character, use `X` as names for the result
#'   unless it had names already (see [vapply()])
#'
#' @return a 3D array onto which the function has been applied to each scan
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
#' sL <- simunet(Adj = Adj,samp.effort = samp.effort,mode = "upper",n.scans = 120L)
#' sL |> sLapply(\(scan) {scan[1,2] <- NA;scan})
sLapply <- function(sL,.f,...,USE.NAMES = TRUE) {
  vapply(
    X = 1:(dim(sL)[3]),
    FUN = function(x) .f(sL[,,x]),
    FUN.VALUE = sL[,,1],
    ... = ...,
    USE.NAMES = USE.NAMES
  )
}

#' Shortcut to a `lapply` equivalent to apply a function to a list of `scanList`: a `sLlist` object
#' Written analogously to [lapply()]
#'
#' @param sLlist a `sLlist` object, a list of `scanList` objects (see
#'   [`perform_exp()`][perform_exp()])
#' @param FUN function, to be applied to each `scanList` objects in `sLlist`
#' @param ... extra argument to be passed, notably named arguments used by `FUN` (see [lapply()])
#'
#' @return a `sLlist` object, a list of `scanList` objects on which the function `FUN` has been
#'   applied (see [`perform_exp()`][perform_exp()])
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
#' # Designing the experiments:
#' ## setting a constant probability of not observing edges
#' group.scan <- design_exp(customize_sampling(method = "group",sampling = 0.8))
#'
#' ## setting an even focal sampling
#' focal.scan <- design_exp(customize_sampling(method = "focal",sampling = "even"))
#'
#' sL <- simunet(Adj = Adj,samp.effort = samp.effort,mode = "upper",n.scans = 120L)
#'
#' sL |> perform_exp(group.scan,focal.scan) |> sLlapply(attrs,a = "edge.Prob")
sLlapply <- function(sLlist,FUN,...) {
  sLlist <- lapply(sLlist,FUN,...)
  class(sLlist) <- "sLlist"
  sLlist
}

#' Print method for `scanList` objects
#' @export
#' @noRd
print.scanList <- function(x,...) {
  print.default(without_attrs(x))
  cat("\n\nHidden attributes:",names(get_attrs(x)))
}

#' transpose method for `scanList` objects
#' @export
#' @noRd
t.scanList <- function(x) {
  aperm(x,c(2,1,3))
}

rbind_2scanList <- function(sL1,sL2) {
  if (!identical(dim(sL1)[1:2],dim(sL2)[1:2]))
    stop("Incompatible dimensions (not the same number of nodes?")
  if (!identical(dimnames(sL1)[1:2],dimnames(sL2)[1:2]))
    warning("scanLists have different node names.")

  if (is.null(dimnames(sL1)[[3]]) | is.null(dimnames(sL2)[[3]]))
    dn <- c(dimnames(sL1)[1:2],list(NULL))
  else
    dn <- c(dimnames(sL1)[1:2],list(c(dimnames(sL1)[[3]],dimnames(sL2)[[3]])))

  array(c(sL1,sL2),
        dimnames = dn,
        dim = c(dim(sL1)[1:2],attrs(sL1,"n.scans") + attrs(sL2,"n.scans"))
  )
}

#' rbind method for `scanList` objects
#' @export
#' @noRd
rbind.scanList <- function(...,deparse.level = 1) {
  Reduce(rbind_2scanList,list(...))
}