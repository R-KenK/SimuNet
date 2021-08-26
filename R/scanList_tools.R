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
  attrs(empiscanList,"theoretical.scanList") <- scan.list
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
  class(sL) <- c("raw","scanList")
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
#' @param copy.class logical, should the class be copied too if `to` is not a `scanList`?
#'
#' @return a `scanList` object, the 3D array containted in `to` with `from`'s `attrs`
#' @export
#'
#' @keywords internal
copy_attrs_to <- function(from,to,copy.class = TRUE) {
  if (!inherits(to,"scanList") & copy.class) {class(to) <- class(from)}
  attr(to,"attrs") <- attrs(from)
  to
}

#' Shortcut to a `lapply` equivalent to apply a function to each 2D matrix contained in a `scanList`
#' Written analogously to [vapply()]. Values returned by `.f` should be a similarly dimensionned
#' matrix as the first one contained in the 3D array
#'
#' @param sL a `scanList` object (see [`simunet()`][simunet()])
#' @param FUN a function, to apply to each 2D matrix contained in `sL`
#' @param ... extra argument to be passed, notably named arguments used by `.f` (see [lapply()])
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
sLapply <- function(sL,FUN,...) {
  sL.ori <- sL
  sL <-
    lapply(
      X = 1:(dim(sL)[3]),
      FUN = function(x) FUN(sL[,,x],...)
    ) |> matList2array()
  sL <- copy_attrs_to(sL.ori,sL)
  sL
}

#' Shortcut to a `lapply` equivalent to apply a function to each 2D matrix contained in a `scanList`
#' Written analogously to [vapply()]. Values returned by `.f` should be a similarly dimensionned
#' matrix as the first one contained in the 3D array
#'
#' @param sL a `scanList` object (see [`simunet()`][simunet()])
#' @param .f a function, to apply to each 2D matrix contained in `sL`
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
#' sL |> sLvapply(\(scan) {scan[1,2] <- NA;scan})
sLvapply <- function(sL,.f,...,USE.NAMES = TRUE) {
  sL.ori <- sL
  sL <-
    vapply(
    X = 1:(dim(sL)[3]),
    FUN = function(x) .f(sL[,,x]),
    FUN.VALUE = sL[,,1],
    ... = ...,
    USE.NAMES = USE.NAMES
  )
  sL <- copy_attrs_to(sL.ori,sL)
  sL
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
#' group.scan <- design_sampling(method = "group",sampling = 0.8)
#'
#' ## setting an even focal sampling
#' focal.scan <- design_sampling(method = "focal",sampling = "even")
#'
#' sL <- simunet(Adj = Adj,samp.effort = samp.effort,mode = "upper",n.scans = 120L)
#'
#' sL |> perform_exp(group.scan,focal.scan) |> sLlapply(attrs,a = "edge.Prob")
sLlapply <- function(sLlist,FUN,...) {
  sLlist <- lapply(sLlist,FUN,...)
  class(sLlist) <- "sLlist"
  sLlist
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

#' Subsetting (dollard-sign) method for `scanList` objects
#' @export
#' @noRd
`$.scanList` <- function(x,name) {
  attrs(x,name)
}

#' Subsetting (single brackets) method for `scanList` objects
#' @export
#' @noRd
`[.scanList` <- function(x,i,j,...) {
  x.ori <- x
  x <- without_attrs(x)
  x.sub <- NextMethod()
  x.sub <- copy_attrs_to(x.ori,x.sub)
  x.sub
}

# printing related functions ----

## printing methods ----

#' Print method for `scanList` objects
#' @export
#' @noRd
print.scanList <- function(x,...) {
  print_sLarray(x)
  format_attributes(x,...)
  invisible(x)
}

#' Print method for `weightedAdj` objects
#' @export
#' @noRd
print.weightedAdj <- function(x,...) {
  to.print <- without_attrs(x)
  class(to.print) <- NULL
  print_clean_scan(to.print,"Weighted adjacency matrix",...)
  format_attributes(x,...)
  invisible(x)
}

#' Print method for `edgeProb` objects
#' @export
#' @noRd
print.edgeProb <- function(x,...) {
  print(x$P)
  invisible(x)
}

#' Print method for `edgeProbMat` objects
#' @export
#' @noRd
print.edgeProbMat <- function(x,digits = 3,...) {
  to.print <- x |> round(digits = digits)
  class(to.print) <- NULL
  print_clean_scan(to.print,"Edge presence probability matrix",...)
  format_attributes(x,...)
  invisible(x)
}

#' Print method for `scaled` objects
#' @export
#' @noRd
print.scaled <- function(x,digits = 3,...) {
  to.print <- without_attrs(x) |> round(digits = digits)
  class(to.print) <- NULL
  print_clean_scan(to.print,"Scaled weighted adjacency matrix",...)
  format_attributes(x,...)
  invisible(x)
}

## printing tools ----

#' Cleaner 3D array print
#'
#' @param sL a `scanList` object
#' @param ... additional arguments to be passed to `Matrix::printSpMatrix()`
#'
#' @return `sL` invisibly, but print a cleaner 3D array via `Matrix::printSpMatrix()`
#' @noRd
print_sLarray <- function(sL,...) {
  if (is.na(dim(sL)[3])) {
    print_clean_scan(sL,"scan:",...)
  } else {
    scan.ind <- choose_scan_to_print(sL)
    truncated <- attr(scan.ind,"truncated")
    # prints all but the last
    lapply(scan.ind,\(s) print_clean_scan(sL[,,s],s,...))
    if (truncated) cat("\n... (",dim(sL)[3] - 3," more scans)\n")
    print_clean_scan(sL[,,dim(sL)[3]],dim(sL)[3],...)
  }
  invisible(sL)
}

#' Choose what scan index to display, truncate if too many
#'
#' @param sL a `scanList` object
#'
#' @return integer vector, indices of scan to print
#' @noRd
choose_scan_to_print <- function(sL) {
  truncated <- dim(sL)[3] > 3
  scan.ind <-
    if (truncated) c(1,2) else 1:(dim(sL)[3] - 1)
  attr(scan.ind,"truncated") <- truncated
  attr(scan.ind,"last.scan") <- dim(sL)[3]
  scan.ind
}

#' Cleaner adjacency matrix print
#'
#' @param mat numeric matrix, a scan
#' @param s integer, scan index
#' @param ... additional arguments to be passed to
#'   [`Matrix::printSpMatrix()`][Matrix::printSpMatrix()]
#' @param mode character, igraph's mode
#' @param col.names logical, see [`Matrix::printSpMatrix()`][Matrix::printSpMatrix()]
#' @param note.dropping.colnames logical, see [`Matrix::printSpMatrix()`][Matrix::printSpMatrix()]
#'
#' @return `mat` invisibly, but print a cleaner scan via
#'   [`Matrix::printSpMatrix()`][Matrix::printSpMatrix()]
#'
#' @importFrom Matrix printSpMatrix
#' @importFrom methods as
#'
#' @noRd
print_clean_scan <- function(mat,s,
                             col.names = FALSE,
                             note.dropping.colnames = FALSE,
                             ...) {
  if (is.numeric(s))
    cat("\nscan: ",s,sep = "")
  else
    cat("\n",s,sep = "")
  m <- mat
  class(m) <- NULL
  methods::as(m,"dgCMatrix") |>
    Matrix::printSpMatrix(col.names = col.names,note.dropping.colnames = note.dropping.colnames,
                          ...)
  invisible(mat)
}

#' Display and format attribute names in attrs if they exist
#'
#' @param x a `scanList` or `weightedAdj` object
#' @param ... ignored
#'
#' @noRd
format_attributes <- function(x,...) {
  if (!is.null(get_attrs(x))) {
    attrs.names <-
      get_attrs(x) |>
      names() |>
      split_returnCarriage_attributes()
    cat("\n\nHidden attributes:\n",attrs.names,"\n",sep = "")
  }
  if (inherits(x,"edgeProbMat")) {
    bet <- attr(x,"Beta priors")
    cat("\n","alpha.prior =",bet[1],"-","beta.prior =",bet[2])
  }
  invisible(x)
}

#' Split character vectors in chunks with at max n.attrs elements
#' separate them with " - ", add a return carriage and makes it into a printable character scalar
#'
#' @param attrs.names character vector, names of a `scanList`'s `attrs`
#'
#' @noRd
split_returnCarriage_attributes <- function(attrs.names,n.attrs = 6) {
  split(attrs.names, ceiling(seq_along(attrs.names) / n.attrs)) |>
    lapply(paste,collapse = " - ") %>%
    {do.call(paste,list(.,collapse = "\n"))}
}
