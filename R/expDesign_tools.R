# expDesign object functions ------------------------------------------------------------------

#' Design experiment to perform on theoretical `scanList`
#'
#' Generate `expDesign` objects that mainly consist on a sequence of functions to apply to
#' `scanList` objects. `expDesign` object generation relies on `purrr`'s
#' [`compose()`][purrr::compose()] function.
#'
#' It is best practice to use a combination of:
#'
#' * provided "building blocks": functions included in the `SimuNet` package, such as
#' [`customize_sampling()`][customize_sampling()], [`add_scans()`][add_scans()],
#' [`remove_mostPeripheral()`][remove_mostPeripheral()], or [`sum_scans()`][sum_scans()]
#' * user-defined functions: designed to take a `scanList` object as a first argument, which is is
#' in essence a 3 dimensional array with adjacency matrices on the first 2 dimensions, and the
#' successive scan numbers as the 3rd dimension. User-defined function should returned a modified
#' `scanList` object, therefore allowing for function chaining in a fashion similar to
#' [`tidy`](https://www.tidyverse.org/) functions (see also the notion of grammar of data
#' manipulation with the [`dplyr`](https://dplyr.tidyverse.org/) package).
#'
#' `expDesign` objects are used to store sequences of manipulations to apply to a theoretical
#' `scanList`, in order to obtain a (simulated) empirical `scanList`. Such manipulations can
#' represent empirical sampling (e.g. group-scan sampling, focal sampling, biased sampling),
#' observation error, node removal, etc.
#'
#' `expDesign` objects can be applied to `scanList` objects via the function
#' [`perform_exp()`][perform_exp()], or directly as [`simunet()`][simunet()] `exp.design` argument
#' for `simunet()` to automatically generate a theoretical `scanList` and apply the inputted
#' `expDesign` to it.
#'
#' Providing more than one `expDesign` object to either [`perform_exp()`][perform_exp()] or
#' [`simunet()`][simunet()] functions (i.e. effectively passing them as their `...` argument) will
#' make them output a *list* of `scanList` objects - i.e. a `sLlist` object - which are convenient
#' ways to apply different experimental manipulation sequences to a given theoretical `scanList` and
#' allow for comparison (e.g. of the impact of sampling regime).
#'
#' @param ... functions of `scanList` objects. Can be a combination of the provided "building
#'   blocks" and user-defined function (cf. the example section).
#' @param .dir character scalar, either `"forward"` (default) or `"backward"` (see `purrr`'s
#'   [`compose()`][purrr::compose()]) to indicate in which order inputted functions should be
#'   applied (first function first, or _vice versa_).
#'
#' @return an `expDesign` objects, a list consisting in the following elements:
#' * `FUN.seq`: function sequence created by `purrr`'s [`compose()`][purrr::compose()] function
#' * WIP: more to come, notably to include more descriptive function names.
#'
#' @export
#'
#' @seealso [perform_exp()], [simunet()], [customize_sampling()].
#'
#' @importFrom purrr compose
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
#' ## setting a biased focal sampling favoring central individual (node strength)
#' focal.scan <- design_exp(
#'   customize_sampling(
#'     method = "focal",
#'     sampling = function(Adj) Adj |>
#'       igraph::graph.adjacency("upper",weighted = TRUE) |>
#'       igraph::strength()
#'   )
#' )
#'
#' ## Adding more scans, removing the most peripheral individual, before performing an even focal
#' ## sampling
#' focal.periRemoved <- design_exp(
#'   function(Adj) add_scans(Adj,42),     # users can use anonymous function to specify arguments
#'   remove_mostPeripheral,               # ... or pass functions as arguments directly
#'   customize_sampling(method = "focal",sampling = "even")    # customize_sampling: special case
#'                                                             # that returns sampling functions
#' )
#'
#' # Apply the experimental design
#' ## on previously obtained theoretical scans
#' sL <- simunet(Adj = Adj,samp.effort = samp.effort,mode = "upper",n.scans = 120L)
#'
#' perform_exp(sL,group.scan)
#' perform_exp(sL,focal.periRemoved) |> sum_scans()
#'
#' ## on newly simulated scanList within the simunet() wrapper
#' simunet(Adj = Adj,samp.effort = samp.effort,mode = "upper",n.scans = 120L,
#'         focal.scan
#' )
#'
#' ## performing a list of experiments
#' perform_exp(sL,group.scan,focal.scan)
#' simunet(Adj = Adj,samp.effort = samp.effort,mode = "upper",n.scans = 120L,
#'         focal.scan,focal.periRemoved
#' )
design_exp <- function(...,.dir = c("forward", "backward")) {
  .dir <- match.arg(.dir)

  dots.call <- as.list(substitute(...()))
  dots.list <- list(...)

  dots.fun <- replace_expDesign_by_funs(...)

  named.list <- do.call(namedList,dots.fun)
  # return(named.list)
  FUN.seq <- do.call(\(...) purrr::compose(...,.dir = .dir),dots.fun)
  generate_expDesign(FUN.seq = FUN.seq,
                     fun.input = named.list,
                     input = dots.call
  )
}

#' Replace ellipsis mixing `expDesign` and functions into flat list of functions
#' Internal
#'
#' @param ... ellipsis argument passed to `design_exp()`
#'
#' @return a flat-list of functions
#'
#' @noRd
replace_expDesign_by_funs <- function(...) {
  dots.call <- as.list(substitute(...()))
  dots.content <- list(...)
  is.expDesign <- dots.content |>
    sapply(inherits,what = "expDesign")
  if (any(is.expDesign)) {
    dots.call[is.expDesign] <-
      lapply(dots.content[is.expDesign],"[[","input")
  }
  unlist(dots.call)
}

#' Perform an experimental design on theoretical `scanList`
#'
#' @param scan.list a `scanList` object. See objects returned by [`simunet()`][simunet()]
#' @param exp.design an `expDesign` object. See objects returned by [`design_exp()`][design_exp()].
#'   If `NULL`, the inputted scan.list is returned as is.
#' @param ... additional `expDesign` object.
#'
#'   If not `NULL`, the different expDesign will be applied to `scan.list` in "parallel": the
#'   returned value will be a list of empirical `scanList`, i.e. a `sLlist` object
#'
#' @return an empirical `scanList` object representing the simulated theoretical scan on which the
#'   experimental manipulations have been applied. Such objects contain:
#'
#'   * the 3 dimensional array representing adjacency matrices (first 2 dimensions) throughout the
#'   different scans (3rd dimension)
#'   * the same `attrs` attribute than the inputted `scan.list` (a list of attributes), in which
#'   `scanList.type = "empirical"`
#'   * a class `empirical`, which inherits from `scanList`
#'   * and other attributes might have been added to `attrs` or modifications depending on
#'   `exp.design` (e.g. [`sum_scans()`][sum_scans()] returns an object of class sum, that inherits
#'   from classes `empirical` or `theoretical`, and further from `scanList`)
#'
#'   If more than one `expDesign` has been inputted via `...`, returns a list of empirical
#'   `scanList`, i.e. a `sLlist` object
#' @export
#'
#' @seealso [design_exp()], [simunet()], [customize_sampling()].
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
#' ## setting a biased focal sampling favoring central individual (node strength)
#' focal.scan <- design_exp(
#'   customize_sampling(
#'     method = "focal",
#'     sampling = function(Adj) Adj |>
#'       igraph::graph.adjacency("upper",weighted = TRUE) |>
#'       igraph::strength()
#'   )
#' )
#'
#' ## Adding more scans, removing the most peripheral individual, before performing an even focal
#' ## sampling
#' focal.periRemoved <- design_exp(
#'   function(Adj) add_scans(Adj,42),     # users can use anonymous function to specify arguments
#'   remove_mostPeripheral,               # ... or pass functions as arguments directly
#'   customize_sampling(method = "focal",sampling = "even")    # customize_sampling: special case
#'                                                             # that returns sampling functions
#' )
#'
#' # Apply the experimental design
#' ## on previously obtained theoretical scans
#' sL <- simunet(Adj = Adj,samp.effort = samp.effort,mode = "upper",n.scans = 120L)
#'
#' perform_exp(sL,group.scan)
#' perform_exp(sL,focal.periRemoved) |> sum_scans()
#'
#' ## performing a list of experiments
#' perform_exp(sL,group.scan,focal.scan)
perform_exp <- function(scan.list,exp.design = NULL,...){
  if (!inherits(scan.list,"scanList")) {
    stop("scan.list inputted is not a scanList object.")
  }

  if (is.null(exp.design)) {
    return(scan.list)
  } else if (!inherits(exp.design,"expDesign")) {
    stop("exp.design inputted is not a expDesign object.")
  }

  if (missing(...)) {
    generate_empiscanList(scan.list,exp.design)
  } else {
    expD.list <- list(exp.design,...)
    sL.list <- lapply(expD.list,\(expD) generate_empiscanList(scan.list = scan.list,exp.design = expD))
    class(sL.list) <- "sLlist"
    sL.list
  }
}

#' Generate `expDesign` objects
#'
#' @param FUN.seq function sequence created by `purrr`'s [`compose()`][purrr::compose()] function
#' @param fun.input functions of `scanList` objects. Can be a combination of the provided "building
#'   blocks" and user-defined function (cf. the example section).
#'
#' @return an `expDesign` objects, a list consisting in the following elements:
#' * `FUN.seq`: function sequence created by `purrr`'s [`compose()`][purrr::compose()] function
#' * WIP: more to come, notably to include more descriptive function names.
#'
#' @keywords internal
generate_expDesign <- function(FUN.seq,fun.input,input) {
  expD <-
    list(
      FUN.seq = FUN.seq,
      # FUN.names = if (!is.null(names(fun.input))) names(fun.input) else attr(fun.input,"FUN.names"),
      FUN.names = names(fun.input),
      rest = fun.input,
      input = input
    )
  class(expD) <- "expDesign"
  expD
}

#' Helper to manage function names in list of functions
#' modified from <https://stackoverflow.com/a/16951524>
#'
#' @param ... functions of `scanList` objects. Can be a combination of the provided "building
#'   blocks" and user-defined function (cf. the example section).
#'
#' @return a named list of functions
#' @export
#'
#' @importFrom stats setNames
#'
#' @keywords internal
namedList <- function(...) {
  L <- list(...)
  snm <- sapply(substitute(list(...)),deparse)[-1]
  if (is.null(nm <- names(L))) nm <- snm
  if (any(nonames <- nm=="")) nm[nonames] <- snm[nonames]
  nm <- rename_sampling(nm,...)
  nm <- rename_customFUN(nm)
  stats::setNames(L,nm)
}

#' Rename sampling functions' names
#' @noRd
rename_sampling <- function(nm,...) {
  L <- list(...)
  is.sampling <- sapply(L,\(fn) !is.null(attr(fn,"FUN.names")))
  nm[is.sampling] <- sapply(L[is.sampling],\(fn) attr(fn,"FUN.names"))
  nm
}

#' Rename custom functions' names
#' @noRd
rename_customFUN <- function(nm) {
  custom.f <- substr(nm,1,9) %>% {. %in% c("function(","c(\"functi")}
  if (any(custom.f)) {
    unique.nm <- data.frame(
      fun = as.character(unique(nm[custom.f])),
      nm = paste0("custom.FUN_",1:length(unique(nm[custom.f])))
    )
    nm[custom.f] <- unique.nm$nm[
      match(nm[custom.f],unique.nm$fun)
    ]
  }
  nm
}

#' Print method for `expDesign` objects
#' @export
#' @noRd
print.expDesign <- function(x,...) {
  format_FUN.names(x$FUN.names)
  cat("\nSee `$FUN.seq` for the functions source code.\n")
}

#' Format function sequence names
#' @noRd
format_FUN.names <- function(FUN.names) {
  paste0(FUN.names," ->") |>
    paste(collapse = "\n ") %>%
    c("Theoretical scanList ->\n ------------\n ",
      .,
      "\n ------------\n empirical scanList\n") |>
    cat(sep = "")
}
