#' `edgeProb` object generator
#'
#' @param Adj integer matrix, the adjacency matrix (see [`simunet()`][simunet()])
#' @param samp.effort integer scalar, the sampling effort (see [`simunet()`][simunet()])
#' @param mode character scalar, the network igraph's mode (see [`simunet()`][simunet()])
#' @param Adj.subfun function, the matrix subsetting function relevant for the adjacency matrix
#'   `mode` (see [`simunet()`][simunet()])
#'
#' @return an `edgeProb` object, i.e. a list containing:
#' * `P`: the edge presence probability matrix
#' * `Adj`: the inputted `Adj`
#' * `samp.effort`: the inputted `samp.effort`
#' * `mode`: the inputted `mode`
#' * `Adj.subfun`: the inputted `Adj.subfun`
#'
#' @export
#'
#' @keywords internal
generate_edgeProb <- function(Adj,samp.effort,mode,
                                  Adj.subfun = NULL){
  Adj[] <- as.integer(Adj)
  if (is.null(Adj.subfun)) {
    Adj.subfun <- determine_Adj.subfun(mode = mode)
  }

  edge.Prob <- list(
    P = draw_edgeProb(Adj = Adj,samp.effort = samp.effort,Adj.subfun = Adj.subfun),
    Adj = Adj,
    samp.effort = samp.effort,
    mode = mode,
    Adj.subfun = Adj.subfun
  )
  class(edge.Prob)<- "edgeProb"
  edge.Prob
}

#' Determine if an `edgeProb` object should be generated from inputs
#'
#' @param Adj integer matrix, the adjacency matrix (see [`simunet()`][simunet()])
#' @param mode character scalar, the network igraph's mode (see [`simunet()`][simunet()])
#' @param samp.effort integer scalar, the sampling effort (see [`simunet()`][simunet()])
#' @param edge.Prob optional, an `edgeProb` object (see [`generate_edgeProb()`][generate_edgeProb()])
#'
#' @return an `edgeProb` object, i.e. a list containing:
#' * `P`: the edge presence probability matrix
#' * `Adj`: the inputted `Adj`
#' * `samp.effort`: the inputted `samp.effort`
#' * `mode`: the inputted `mode`
#' * `Adj.subfun`: the inputted `Adj.subfun`
#' @noRd
determine_edgeProb <- function(Adj = NULL,mode = NULL,samp.effort = NULL,edge.Prob = NULL) {
  if (!is.null(edge.Prob)) {
    if (!is.null(Adj) & !is.null(samp.effort)) {
      stop("Both Adj/samp.effort and edge.Prob have been provided. Please provide either.")
    }
    edge.Prob
  } else {
    if (is.null(Adj) | is.null(samp.effort)) {
      stop("Adj or samp.effort missing. Either provide both, or provide an edge.Prob.")
    }
    generate_edgeProb(Adj = Adj,mode = mode,samp.effort = samp.effort)
  }
}


#' Draw edge presence probability matrix from posterior Beta distribution
#' (Internal use)
#'
#' Edge presence probabilities are drawn from a posterior Beta distribution
#' \eqn{Beta(\alpha,\beta)}, in which parameters \eqn{\alpha} and \eqn{\beta} correspond to total
#' (pseudo-)counts of the times when an edge was 1 and 0, respectively. By default, an uninformative
#' prior (Jeffrey's prior, i.e. \eqn{\alpha_{prior} = \beta_{prior} = 1 / 2}{\alpha.prior =
#' \beta.prior = 0.5}) is used, and is added to the observed edge weights in `Adj`, so that:
#' \deqn{\alpha = Adj + \alpha_{prior}, \beta = N_{samp.effort} - Adj + \beta_{prior}}{\alpha = Adj
#' + \alpha.prior, \beta = samp.effort - Adj + \beta.prior} where \eqn{N_{samp.effort} -
#' Adj}{samp.effort - Adj} is a positive or null integer.
#'
#' For Bayesian inference and conjugation, a prior beta distribution
#' \eqn{Beta(\alpha_{prior},\beta_{prior})}{Beta(\alpha.prior,\beta.prior)} is used, with default
#' \deqn{\alpha_{prior} = \beta_{prior} = 1 / 2}{\alpha.prior = \beta.prior = 0.5} which corresponds
#' to Jeffrey's prior. Alternative parametrization can rely on \deqn{\alpha_{prior} = \beta_{prior}
#' = 1}{\alpha.prior = \beta.prior = 1} for a prior beta distribution equivalent to a uniform
#' distribution over \eqn{[0,1]}.
#'
#' In [`simunet()`][simunet()], a new `edgeProb` is drawn before simulating binary scans. Two
#' `scanList`s outputted from an identical `Adj` matrix and `samp.effort` will result in different
#' `edgeProb`s, but these `edgeProb`s are drawn from the same posterior Beta distribution
#' \eqn{Beta(\alpha,\beta)}.
#'
#' The `edgeProb` used is stored in the `scanList` outputted by `simunet()` as an attribute in the
#' attributes list `attrs`, and can be retrieved via `attrs(scan.list,"edgeProb")`).
#'
#' This procedure is equivalent to drawing a `scanList` from a Beta-Binomial distribution
#' \eqn{BetaBinom(n_{scans},\alpha,\beta)}{BetaBinom(n.scans,\alpha,\beta)}, but allows
#' "decomposing" the simulated edge weights into a list of binary scans instead of outputting only a
#' new weighted adjacency matrix.
#'
#' @param Adj integer matrix, the adjacency matrix (see [`simunet()`][simunet()])
#' @param Adj.subfun function, the matrix subsetting function relevant for the adjacency matrix
#'   `mode` (see [`simunet()`][simunet()])
#' @param alpha.prior positive numeric scalar, the parameter alpha (added to shape1 in `rbeta()`)
#'   used in the prior beta distribution. See [`rbeta()`][stats::rbeta()]
#' @param beta.prior  positive numeric scalar, the parameter beta (added to shape2 in `rbeta()`)
#'   used in the prior beta distribution. See [`rbeta()`][stats::rbeta()]
#' @param samp.effort integer scalar, the sampling effort (see [`simunet()`][simunet()])
#'
#' @return numeric matrix of edge presence probabilities
#'
#' @importFrom stats rbeta
#'
#' @seealso [simunet()].
#'
#' @export
#'
#' @examples
#' # Internally used in `generate_edgeProb()`, itself used in `simunet()`
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
#' sL <- simunet(Adj = Adj,samp.effort = samp.effort,mode = "directed",n.scans = 20L)
#' attrs(sL,"edge.Prob")

draw_edgeProb <- function(Adj,samp.effort,
                          Adj.subfun = NULL,
                          alpha.prior = 0.5,
                          beta.prior = 0.5){
  if (samp.effort < max(Adj)) {stop("samp.effort provided incompatible with the maximum value found in the provided adjacency matrix.")}
  P <- Adj
  P[Adj.subfun(P)] <-
    Adj[Adj.subfun(Adj)] |>
    {\(x) stats::rbeta(length(x),shape1 = x + alpha.prior,shape2 = samp.effort - x + beta.prior)}()
  attr(P,"Beta priors") <- c(alpha.prior,beta.prior)
  P
}

#' Convenience function to reconstruct an `edgeProb` object from what has been stored in `scanList` object's attributes
#' Avoid `edgeProb` stored in `scanList` to store duplicates of `Adj`, `samp.effort`, `mode`, and `Adj.subfun`
#'
#' @param scan.list a `scanList` object (see [`simunet()`][simunet()])
#'
#' @return an `edgeProb` object, i.e. a list containing:
#' * `P`: the edge presence probability matrix
#' * `Adj`: the inputted `Adj`
#' * `samp.effort`: the inputted `samp.effort`
#' * `mode`: the inputted `mode`
#' * `Adj.subfun`: the inputted `Adj.subfun`
#'
#' @export
#'
#' @keywords internal
reconstruct_edgeProb <- function(scan.list) {
  edge.Prob <- list(
    P = attrs(scan.list,"edge.Prob"),
    Adj = attrs(scan.list,"Adj"),
    samp.effort = attrs(scan.list,"samp.effort"),
    mode = attrs(scan.list,"mode"),
    Adj.subfun = attrs(scan.list,"Adj.subfun")
  )
  class(edge.Prob)<- "edgeProb"
  edge.Prob
}
