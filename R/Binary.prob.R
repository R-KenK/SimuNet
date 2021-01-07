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
#' @export
#'
#' @examples
#' set.seed(42)
#'
#' n<- 5;nodes<- letters[1:n];
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- sample(0:30,n*(n-1),replace = TRUE)
#' Adj
#'
#' Binary.prob(Adj,42)
Binary.prob<- function(Adj,total_scan,
                       mode = c("directed", "undirected", "max","min", "upper", "lower", "plus","vector")){
  if(total_scan<max(Adj)) {stop("total_scan provided incompatible with the maximum value found in the provided adjacency matrix.")}
  mode<- match.arg(mode)

  P<- Adj
  Adj.subfun<- switch(mode,
                      "directed" = ,
                      "undirected" = ,
                      "max" = ,
                      "min" = ,
                      "plus" = non.diagonal,
                      "upper" = upper.tri,
                      "lower" =  lower.tri,
                      "vector" = function(V){rep(TRUE,length(V))}
  )
  min_resol<- 1/total_scan;
  prob.scaled<- (Adj[Adj.subfun(Adj)]*(1-2*min_resol)/total_scan)+min_resol;
  P[Adj.subfun(P)]<- prob.scaled
  P
}
