# empirical sampling methods ----

## Wrapper to create sampling functions ----

#'  Customize a sampling regime to perform on a theoretical `scanList`
#'
#'  This function returns a function tailored to be used as part of an `expDesign` object (see
#'  [`design_exp()`][design_exp()]). It is written as a convenient wrapper for commonly used
#'  sampling methods: group-scan sampling and focal-scan sampling.
#'
#'  the function accepts as `sampling` parameter:
#'  * character scalar: common options like random edge observation probability or even focal
#'  sampling
#'  * for `method = "group"`: numeric scalar (constant) or matrix representing edge observation
#'  probabilities
#'  * user-defined functions of the adjacency matrix `Adj` that returns either an edge observation
#'  probability matrix, or a vector of the probabilities of selecting a given node at each focal
#'  scan. If the user-defined function returns invalid probabilities e.g.:
#'    * a value > 1 for `method = "group"`: the function tries to rescale values via `scales`
#'    package's [`rescale_max()`][scales::rescale_max()] function
#'    * some probabilities of being a focal = 0 for `method = "focal"`: the function adds the
#'    non-null minimum probability to all probabilities (values > 1 should be handled correctly as
#'    the `prob` argument of the [`sample()`][sample()] function)
#'
#'  The empirical sampling works by replacing unobserved edges by `NA`s in the 3D array, either:
#'  * because a given edge hasn't been observed during the group-scan sampling
#'  * or because the masked edge was not involving the focal node during the scan
#'
#'  Convenience "building blocks" functions - respectively [`count_NA()`][count_NA()] and
#'  [`count_nonNA()`][count_nonNA()] - can be used to count masked and sampled edges throughout the
#'  whole simulation.
#'
#'  New attributes are added to `attrs`:
#'  * in the case of `method = "group"`:
#'    * `obs.P`: matrix of probabilities of observing an edge (whether it is 0 or 1)
#'  * in the case of `method = "focal"`:
#'    * `focalList`: named integer vector representing the node's index (row/column) to be sampled
#'    for each scan. Names are obtain from the adjacency matrix `Adj`, the vector's length is equal
#'    to `n.scans`
#'
#' @param method character scalar, either `"group"` or `"focal"`
#' @param sampling depending on chosen `method`, users should input either:
#'   * a numeric scalar (`"constant"`): the constant probability of observing an edge for all edges
#'   * a numeric matrix (`"matrix"`): the probabilities of observing an edge for each edges
#'   * a character scalar: for common sampling regimes:
#'     * `"even"`: in the case of `method = "focal"`: select focals as evenly as possible, and the
#'     extra scans uniformly
#'     * `"random"`: random edge observation probabilities or uniform probability of choosing a
#'     focal at each scan
#'   * a user-defined function (`"function"`): a function of the adjacency matrix `Adj` (can be
#'   named anything) that:
#'     * in the case of `method = "group"`: returns a matrix of the probabilities of observing an
#'     edge for each edges
#'     * in the case of `method = "focal"`: returns a vector of the probabilities of choosing a
#'     focal node at each scan
#'   * WIP: more option to be added, like with the possibility to pass a `focalList` object directly
#' @param all.sampled logical scalar, should all nodes be sampled at least once? (TO CHECK: does it
#'   work with group-scan sampling?)
#'
#' @return a function of a theoretical `scan.list` simulating the empirical sampling of the
#'   network's edges at each scan. To be used as part of an `expDesign` object (see
#'   [`design_exp()`][design_exp()])
#'
#' @export
#'
#' @seealso [simunet()], [design_exp()], [perform_exp()], [group_sample()], [determine_obsProb()],
#'   [focal_sample()], [draw_focalList()], [mask_non.focals()], [count_NA()], [count_nonNA()].
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
#' # Designing sampling regimes:
#' ## setting a constant probability of not observing edges
#' group.constant <- design_sampling(method = "group",sampling = 0.8)
#'
#' ## setting a random probability of not observing edges
#' group.random <- design_sampling(method = "group",sampling = "random")
#'
#' ## setting probability of not observing edges via user-defined functions
#' g.fun1 <- function(Adj) Adj     # observation proportional to the network's weights,
#'                                 # will be rescaled as probabilities internally
#' group.fun1 <- design_sampling(method = "group",sampling = g.fun1)
#'
#' ### user-defined functions can also be passed as anonymous functions
#' group.fun2 <- design_sampling(method = "group",sampling = function(Adj) Adj^2)
#'
#' ## evenly select focals
#' focal.even <- design_sampling(method = "focal",sampling = "even")
#'
#' ## randomly select focals
#' focal.random <- design_sampling(method = "focal",sampling = "random")
#'
#' ## setting probability of selecting focals via user-defined functions
#' f.fun1 <- function(Adj) 1:nrow(Adj)       # linear increase of probability of being focal,
#'                                           # akin to a linear trait
#' focal.fun1 <- design_sampling(method = "focal",sampling = f.fun1)
#'
#' ### user-defined functions can also be passed as anonymous functions
#' focal.fun2 <- design_sampling(method = "focal",sampling = function(Adj) Adj |>
#'                                    igraph::graph.adjacency(mode = "upper",weighted = TRUE) |>
#'                                    igraph::eigen_centrality() |> {\(x) x$vector}()
#' )                            # probabilities proportional to nodes' eigen-vector centralities
#'
#' # Design and run experiment based on these sampling regime
#' ## sampling regimes can be included in an `expDesign` object and passed to `simunet()`...
#' g.const.exp <- design_exp(group.constant)
#' simunet(Adj = Adj,samp.effort = samp.effort,mode = "upper",n.scans = 120L,g.const.exp)
#'
#' ## ... or passed to `perform_exp()`...
#' g.rand.periRemoved <- design_exp(remove_mostPeripheral,group.random)
#'
#' sL <- simunet(Adj = Adj,samp.effort = samp.effort,mode = "upper",n.scans = 120L)
#' sL |> perform_exp(g.rand.periRemoved)
#'
#' ## ... or used "in situ" in either `simunet()` or `perform_exp()`,
#' ## but need to be passed to `design_exp()`
#' ## (TO DO: recognize sampling regime and manage this automatically)
#' simunet(Adj = Adj,samp.effort = samp.effort,mode = "upper",n.scans = 120L,design_exp(group.fun2))
#' sL |> perform_exp(focal.even)
#' sL |> perform_exp(design_sampling("focal","random"))
design_sampling <- function(method = c("group","focal"),
                               sampling =  c("constant","matrix","even","random","function"),
                               all.sampled = TRUE) {

  method <- match.arg(method)
  samp.type <- determine_sampling_type(method = method,sampling = sampling)

  FUN <- switch(method,
         "group" = \(scan.list) group_sample(scan.list = scan.list,
                                             sampling = sampling,
                                             all.sampled = all.sampled
         ),
         "focal" = \(scan.list) focal_sample(scan.list = scan.list,
                                             sampling = sampling,
                                             all.sampled = all.sampled
         )
  )
  attr(FUN,"FUN.names") <- paste0(method,"-scan sampling: ",samp.type)
  FUN.seq <- purrr::compose(FUN)
  generate_expDesign(FUN.seq = FUN.seq,
                     fun.input = namedList(FUN),
                     input = FUN
  )
}

#' Checks if the method and sampling parameter combination is adequate
#'
#' @param method inputted `method` parameter, see [`design_sampling()`][design_sampling()]
#' @param sampling inputted `sampling` parameter, see [`design_sampling()`][design_sampling()]
#'
#' @return character, the type of sampling. Returns an error if incompatible method/sampling
#'   combination
#' @export
#'
#' @keywords internal
determine_sampling_type <-
  function(method,
           sampling =  c("constant","matrix","even","random","function")) {
    switch(
      method,
      "group" =
        switch(
          class(sampling)[1],
          "character" = {
            sampling <- match.arg(sampling)
            switch(
              sampling,
              "even" = stop('Invalid `sampling` option for method = "group" (cf. `?design_sampling`)'),
              "random" =,
              "constant" =, # this is the default choice when no sampling argument is inputted
              "matrix" =,
              "function" = sampling
            )
          },
          "function" = "function","numeric" = "constant","matrix" = "matrix"
        ),
      "focal" =
        switch(
          class(sampling)[1],
          "character" = {
            sampling <- match.arg(sampling)
            switch(
              sampling,
              "constant" =, # this is the default choice when no sampling argument is inputted
              "matrix" = stop('Invalid `sampling` option for method = "focal" (cf. `?design_sampling`)'),
              "random" =,
              "even" =,
              "function" = sampling
            )
          },
          "matrix" =,
          "numeric" = stop('Invalid `sampling` option for method = "focal" (cf. `?design_sampling`)'),
          "function" = "function"
        )
    )
  }


## Group-scan sampling ----

#' Performs a group-scan sampling over a `scanList` object
#' Internal.
#'
#' @param scan.list a `scanList` object
#' @param sampling for `method = "group`, users should input either:
#'   * a numeric scalar (`"constant"`): the constant probability of observing an edge for all edges
#'   * a numeric matrix (`"matrix"`): the probabilities of observing an edge for each edges
#'   * a character scalar: for common sampling regimes:
#'     * `"random"`: random edge observation probabilities
#'   * a user-defined function (`"function"`): a function of the adjacency matrix `Adj` (can be
#'   named anything) that returns a matrix of the probabilities of observing an edge for each edges
#' @param all.sampled logical scalar, should all nodes be sampled at least once? (TO CHECK: does it
#'   work with group-scan sampling?)
#'
#' @return  an empirical `scanList` object in which, compared to the `theoretical.scanList` (added
#'   to `attrs`), unobserved edges are replaced by `NA`s (regardless of them being 0 or 1).
#'
#' Returned `scanList` has new attributes added to attrs:
#' * `obs.P`: matrix of probabilities of observing an edge (whether it is 0 or 1)
#' * `theoretical.scanList`: the original theoretical `scanList` from which some edges have not been
#' observed
#'
#' @export
#'
#' @seealso [design_sampling()], [determine_obsProb()].
#'
#' @importFrom stats rbinom
#'
#' @keywords internal
group_sample <- function(scan.list,sampling = c("constant","matrix","random","function"),
                         all.sampled = TRUE){

  Adj <- attrs(scan.list,"Adj")
  sf <- attrs(scan.list,"Adj.subfun")
  mode <- attrs(scan.list,"mode")

  obs.P <- determine_obsProb(scan.list = scan.list,sampling = sampling,all.sampled = all.sampled)

  groupSampled <-
    sLvapply(scan.list,
            \(s) {
              s[sf(s)] <-
                s[sf(s)] |>
                {\(x) ifelse(stats::rbinom(x,1L,obs.P[sf(obs.P)]) == 1,x,NA)}()
              s |>
                resolve_NA(mode = mode)
            }
    )

  groupSampled <- copy_attrs_to(scan.list,groupSampled)
  attrs(groupSampled,"obs.P") <- obs.P
  groupSampled
}

#' Determine the matrix of probabilities of observing the edges
#' Internal.
#'
#' @param scan.list a `scanList` object
#' @param sampling for `method = "group`, users should input either:
#'   * a numeric scalar (`"constant"`): the constant probability of observing an edge for all edges
#'   * a numeric matrix (`"matrix"`): the probabilities of observing an edge for each edges
#'   * a character scalar: for common sampling regimes:
#'     * `"random"`: random edge observation probabilities
#'   * a user-defined function (`"function"`): a function of the adjacency matrix `Adj` (can be
#'   named anything) that returns a matrix of the probabilities of observing an edge for each edges
#' @param all.sampled logical scalar, should all nodes be sampled at least once? (TO CHECK: does it
#'   work with group-scan sampling?)
#'
#' @return an `obsProb` object, being:
#' * a matrix of probabilities of observing an edge (whether it is 0 or 1)
#' * with attribute `"sampling"`, one of `"constant"`,`"matrix"`,`"random"`,`"function"`,
#' accordingly
#' @export
#'
#' @importFrom scales rescale_max
#' @keywords internal
determine_obsProb <- function(scan.list,sampling = c("constant","matrix","random","function"),all.sampled) {
  Adj <- attrs(scan.list,"Adj")
  sf <- attrs(scan.list,"Adj.subfun")

  if (is.function(sampling)) {
    sampling_fun <- sampling
    sampling <- "function"

    obs.P <- sampling_fun(Adj)
    attr(obs.P,"sampling_fun") <- sampling_fun
  } else {
    switch(typeof(sampling),
           "integer" =,
           "double" = {
             n <- nrow(Adj)
             obs.P <- matrix(sampling,nrow = n,ncol = n,dimnames = dimnames(Adj))
             if (inherits(sampling,"array")) sampling <- "matrix"
             else sampling <- "constant"
           },
           "character" = {
             sampling <- match.arg(sampling)
             switch(sampling,
                    "constant" = stop("Please input a `sampling` variable as a numeric constant/scalar directly (cf. `?group_sample`"),
                    "matrix" = stop("Please input a `sampling` variable as a numeric matrix directly (cf. `?group_sample`"),
                    "function" = stop("Please input a `sampling` variable as a function object directly (cf. `?group_sample`"),
                    "random" = {
                      obs.P <- Adj
                      obs.P[sf(obs.P)] <- runif(length(obs.P[sf(obs.P)]))
                    },
                    # more group-scan sampling "types" can be added here
                    stop("Inputted `sampling` variable not recognized.")
             )
           }
    )
  }
  if (any(obs.P[sf(obs.P)] > 1)) obs.P[sf(obs.P)] <- scales::rescale_max(obs.P[sf(obs.P)])
  obs.P[!sf(obs.P)] <- 0
  attr(obs.P,"sampling") <- sampling
  obs.P
}

## Focal-scan sampling ----

#' Performs a focal-scan sampling over a `scanList` object
#' Internal.
#'
#' @param scan.list a `scanList` object
#' @param sampling for `method = "focal`, users should input either:
#'   * a character scalar: for common sampling regimes:
#'     * `"even"`: select focals as evenly as possible, and the extra scans uniformly
#'     * `"random"`: uniform probability of choosing a focal at each scan
#'   * a user-defined function (`"function"`): a function of the adjacency matrix `Adj` (can be
#'   named anything) that returns a vector of the probabilities of choosing a focal node at each
#'   scan
#'   * WIP: more option to be added, like with the possibility to pass a `focalList` object directly
#' @param all.sampled logical scalar, should all nodes be sampled at least once? (TO CHECK: does it
#'   work with group-scan sampling?)
#'
#' @return  an empirical `scanList` object in which, compared to the `theoretical.scanList` (added
#'   to `attrs`), edges not involving the scan's focal are replaced by `NA`s (regardless of them
#'   being 0 or 1).
#'
#' Returned `scanList` has new attributes added to attrs:
#' * `focalList`: named integer vector representing the node's index (row/column) to be sampled for
#' each scan. Names are obtain from the adjacency matrix `Adj`, the vector's length is equal
#' * `theoretical.scanList`: the original theoretical `scanList` from which some edges have not been
#' observed
#'
#' @export
#' @seealso [design_sampling()], [draw_focalList()].
#'
#' @keywords internal
focal_sample <- function(scan.list,sampling = c("even","random","function"),all.sampled = TRUE){
  focalList <- draw_focalList(scan.list = scan.list,
                              sampling = sampling,
                              all.sampled = all.sampled
  )
  focalSampled <- mask_non.focals(scan.list = scan.list,focalList = focalList)
  focalSampled <- copy_attrs_to(scan.list,focalSampled)
  attrs(focalSampled,"focalList") <- focalList
  focalSampled
}

#' Draw the list of focals to sample at each scan
#'
#' @param scan.list a `scanList` object
#' @param sampling for `method = "focal`, users should input either:
#'   * a character scalar: for common sampling regimes:
#'     * `"even"`: select focals as evenly as possible, and the extra scans uniformly
#'     * `"random"`: uniform probability of choosing a focal at each scan
#'   * a user-defined function (`"function"`): a function of the adjacency matrix `Adj` (can be
#'   named anything) that returns a vector of the probabilities of choosing a focal node at each
#'   scan
#'   * WIP: more option to be added, like with the possibility to pass a `focalList` object directly
#' @param all.sampled logical scalar, should all nodes be sampled at least once? (TO CHECK: does it
#'   work with group-scan sampling?)
#'
#' @return named integer vector representing the node's index (row/column) to be sampled for each
#'   scan. Names are obtain from the adjacency matrix `Adj`, the vector's length is equal to
#'   `n.scans`
#' @export
#'
#' @keywords internal
draw_focalList<- function(scan.list,sampling = c("even","random","function"),all.sampled = TRUE) {
  if (is.function(sampling)) {
    sampling_fun <- sampling
    sampling <- "function"
  } else {
    sampling_fun <- NULL
    sampling <- match.arg(sampling)
  }

  # retrieve variables from scan.list
  Adj <- attrs(scan.list,"Adj")
  n <- Adj |> nrow()
  node_names <- Adj |> rownames()
  n.scans <- attrs(scan.list,"n.scans")

  # shape future focal.list, filling it with NAs
  focal.list <- rep(NA,n.scans);

  switch(sampling,
         # manage the case of an even focal.list
         "even" = {
           focal.list <-
             c(rep(1:n,n.scans %/% n),sample(1:n,n.scans %% n)) |>
             {\(f) f[sample(seq_along(f))]}()
         },
         # if not even, select at least each each node once, and adjust the rest of the sampling effort needed
         "random" = {
           if (all.sampled){
             if (n > n.scans) stop("n.scans is too small to sample all nodes.")
             focal.list[sample(1:n.scans,n)] <- 1:n; n.scans <- n.scans - n
           }

           focal.list[is.na(focal.list)] <- quick_sample(1:n,n.scans)
           focal.list
         },
         "function" = {
           if (is.null(sampling_fun)) {stop("Please input a `sampling` variable as a function object directly (cf. `?focal_sample`")}
           if (all.sampled) {
             if (n > n.scans) stop("n.scans is too small to sample all nodes.")
             focal.list[sample(1:n.scans,n)] <- 1:n
             n.scans <- n.scans - n
           }
           # applies the user-defined function, adjust the minimum probability to be non zero
           P <- sampling_fun(Adj)
           if (any(P == 0)) P <- P + min(P[P > 0])
           # replace remaining NAs for each scan with a node given their probability distribution at each scan
           focal.list[is.na(focal.list)] <- sample(1:n,n.scans,replace = TRUE,prob = P)
           focal.list
         },
         stop("`sampling` class not recognized.")
  )
  names(focal.list) <- node_names[focal.list] # keeps nodes names `NULL` if `Adj` didn't have nodes names
  attr(focal.list,"sampling") <- sampling
  focal.list
}

#' Mask edges not-involving the scan's focal node (applies the focal list)
#' Internal.
#'
#' @param scan.list a `scanList` object
#' @param focalList named integer vector, returned by draw_focalList representing the node's index
#'   (row/column) to be sampled for each scan. Names are obtain from the adjacency matrix `Adj`, the
#'   vector's length is equal to `n.scans`
#'
#' @return a 3D array where edges not-involving the scan's focal node are replaced by`NA`s
#' @export
#'
#' @keywords internal
mask_non.focals <- function(scan.list,focalList) {
  sf <- attrs(scan.list,"Adj.subfun")

  vapply(1:dim(scan.list)[3],
         \(s) {
           scan <- scan.list[,,s]
           foc <- focalList[s]

           scan[-foc,-foc] <- NA
           scan[!sf(scan)] <- 0L
           scan
         },scan.list[,,1]
  )
}
