#' TO WRITE
#'
#' @param Adj TO WRITE
#' @param samp.effort TO WRITE
#' @param mode TO WRITE
#' @param Adj.subfun TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @examples
#' # TO WRITE
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

#' TO WRITE
#'
#' @param Adj TO WRITE
#' @param mode TO WRITE
#' @param samp.effort TO WRITE
#' @param edge.Prob TO WRITE
#'
#' @return TO WRITE
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


#' Binarize from adjacency matrix
#' Internal use. Provide binary probability for each weight, taking into account the sampling effort.
#'
#' @param Adj square integers matrix of occurrences of dyads. WIP: implement method for association matrices...
#' @param Adj.subfun  TO WRITE
#' @param alpha.prior  TO WRITE
#' @param beta.prior  TO WRITE
#' @param samp.effort integer, sampling effort. Note that 1/samp.effort should be relatively small, increasingly small with increasing precision.
#'
#' @details
#'
#' At the moment the probabilities are calculated to be quasi-proportional to the input weights,
#' with a factor closer to 1/samp.effort as 1/samp.effort become small relatively to 1 and to the weights themselves.
#'
#' The formula used in the function is to constrained output probabilities between:
#' * min_resol = 1/samp.effort instead of 0 (with min_resol supposedly small),
#' * and 1-min_resol instead of 1  (with min_resol supposedly small compared to 1)
#' * cf. Keuk, unpublished yet for further details
#'
#' @return matrix of probability of presence for each dyad
#'
#' @importFrom stats rbeta
#'
#' @noRd
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
