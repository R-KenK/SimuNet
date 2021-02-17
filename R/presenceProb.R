# presence.prob generation and related functions -------------------------------

#' Generator for `presenceProb` objects
#'
#' @param Adj square integers matrix of occurrences of dyads.
#' @param total_scan integer, sampling effort. Note that 1/total_scan should be relatively small, increasingly small with increasing precision. Optional if using presence.prob.
#' @param mode Character scalar, specifies what type of igraph network `mode` should be used to convert the supplied matrix. Possible values are:
#' \itemize{
#'   \item{`"directed"` (the default): for non-symetrical adjacency matrix where `Adj[i,j]` doesn't have the same meaning as `Adj[j,i]`}
#'   \item{`"undirected"`: same as `"max"`}
#'   \item{`"upper"`: undirected matrix focusing only on the upper triangle of `Adj` (relying on `upper.tri`). Either `"upper"` or `"lower"` could be favor if only one of `Adj[i,j]` and `Adj[j,i]` should be randomized}
#'   \item{`"lower"`: undirected matrix focusing only on the lower triangle of `Adj` (relying on `lower.tri`)}
#'   \item{`"max"`: from a `"directed"` randomization process (both `Adj[i,j]` and `Adj[j,i]` will be drawn at each scan), `max(Adj[i,j],Adj[j,i])` will be kept for both}
#'   \item{`"min"`: from a `"directed"` randomization process (both `Adj[i,j]` and `Adj[j,i]` will be drawn at each scan), `min(Adj[i,j],Adj[j,i])` will be kept for both}
#'   \item{`"plus."`:  from a `"directed"` randomization process (both `Adj[i,j]` and `Adj[j,i]` will be drawn at each scan), `Adj[i,j] + Adj[j,i]` will be kept for both}
#'   \item{`"vector"`: experimental. To consider adjacency matrices as flat vectors to be randomized. Relevance unexplored yet.}
#'   \item{See details \link[igraph]{graph_from_adjacency_matrix}}
#' }
#' @param Adj.subfun optional, used internally in scan-related functions. Subsetting function of the adjacency matrix. Driven by igraph `"mode"` argument.
#'
#' @return a `presenceProb` object (S3 class) containing:
#' \itemize{
#'   \item{`P`: a [0,1] numeric matrix of probability of presence of a tie between each dyad}
#'   \item{`Adj`: inputted `Adj`}
#'   \item{`total_scan`: inputted `total_scan`}
#'   \item{`mode`: inputted `mode`}
#'   \item{`Adj.subfun`: inputted `Adj.subfun`}
#' }
#'
#' @export
#'
#' @examples
#' set.seed(42)
#'
#' n<- 5;nodes<- as.character(1:n);total_scan<- 42;
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- sample(0:total_scan,n*(n-1),replace = TRUE)
#' Adj
#'
#' generate_presenceProb(Adj,total_scan,mode = "directed")
generate_presenceProb <- function(Adj,total_scan,mode,
                                  Adj.subfun = NULL){
  if (is.null(Adj.subfun)) {
    Adj.subfun <- determine_Adj.subfun(mode = mode)
  }

  presence.prob <- list(
    P = binary.prob(Adj = Adj,total_scan = total_scan,Adj.subfun = Adj.subfun),
    Adj = Adj,
    total_scan = total_scan,
    mode = mode,
    Adj.subfun = Adj.subfun
  )

  class(presence.prob)<- "presenceProb"
  presence.prob
}

#' Print method for `presenceProb` objects
#' @importFrom Matrix Matrix
#' @importFrom Matrix printSpMatrix
#' @export
#' @noRd
print.presenceProb <- function(x,...){
  use_printSpMatrix(x$P)
  cat("  mode: ",x$mode)
}

#' Test if object if a `presenceProb` object
#'
#' @param x an object to test.
#'
#' @return logical, TRUE if the inputted object is a `presenceProb` object.
#'
#' @noRd
is.presenceProb <- function(x){
  inherits(x,"presenceProb")
}

#' Binarize from adjacency matrix
#' Internal use. Provide binary probability for each weight, taking into account the sampling effort.
#'
#' @param Adj square integers matrix of occurrences of dyads. WIP: implement method for association matrices...
#' @param total_scan integer, sampling effort. Note that 1/total_scan should be relatively small, increasingly small with increasing precision.
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. See also the weighted argument, the interpretation depends on that too. Possible values are: directed, undirected, upper, lower, max, min, plus. See details \link[igraph]{graph_from_adjacency_matrix}. Added also a vector mode.
#'
#' @details
#'
#' At the moment the probabilities are calculated to be quasi-proportional to the input weights,
#' with a factor closer to 1/total_scan as 1/total_scan become small relatively to 1 and to the weights themselves.
#'
#' The formula used in the function is to constrained output probabilities between:
#' * min_resol = 1/total_scan instead of 0 (with min_resol supposedly small),
#' * and 1-min_resol instead of 1  (with min_resol supposedly small compared to 1)
#' * cf. Keuk, unpublished yet for further details
#'
#' @return matrix of probability of presence for each dyad
#'
#' @noRd
binary.prob <- function(Adj,total_scan,
                       Adj.subfun = NULL){
  if (total_scan < max(Adj)) {stop("total_scan provided incompatible with the maximum value found in the provided adjacency matrix.")}
  bin.P <- Adj
  min_resol <- 1 / total_scan;
  prob.scaled <- (Adj[Adj.subfun(Adj)] * (1 - 2 * min_resol) / total_scan) + min_resol;
  bin.P[Adj.subfun(bin.P)] <- prob.scaled
  bin.P
}

