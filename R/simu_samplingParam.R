# sampling.param generation and related functions -------------------------------

#' Wrapper for `samplingParam` objects
#'
#' @param method Character scalar, specifies if the function should return a theoretical perfect group scan, an  empirical group scan (a similarly dimensioned matrix as Adj), or a focal scan (a vector representing the given focal's row in the group scan matrix).
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
#' @param obs.prob an `obsProb` object
#' @param focal a `focal` object. Otherwise, a `focalList` object can be provided alongside a `scan.number`
#' @param scan.number Optional. Only required if inputted `focal` is a `focalList` object. An integer between in `1:total_scan`
#'
#' @return an `samplingParam` object (S3 class) containing:
#' \itemize{
#'   \item{method}{inputted `method`}
#'   \item{mode}{inputted `mode`}
#'   \item{obs.prob}{inputted `obs.prob`}
#'   \item{focal}{inputted `focal`}
#' }
#' @export
#'
#' @examples
#'
#' set.seed(42)
#'
#' n<- 5;nodes<- as.character(1:n);total_scan<- 42;
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- sample(0:total_scan,n*(n-1),replace = TRUE)
#' Adj
#'
#' obs.prob.random<- generate_obs.prob(Adj,"directed",obs.prob_fun = "random")
#' obs.prob.constant<- generate_obs.prob(Adj,"directed",obs.prob_fun = .42)
#' focal.list<- generate_focal.list(Adj,total_scan,focal.prob_fun = "even",all.sampled = TRUE)
#' focal<- generate_focal(focal.list,10)
#'
#' generate_samplingParam(method = "group",obs.prob = obs.prob.random)
#' generate_samplingParam(method = "focal",focal = focal)
#' generate_samplingParam(method = "both",obs.prob = obs.prob.constant,
#'                         focal = focal.list,scan.number = 20)
simu_samplingParam<- function(mode = c("directed","undirected","max","min","upper","lower","plus","vector"),
                              obs.prob = NULL,focal = NULL,scan.number = NULL){
  mode<- match.arg(mode);

  method<- determine_method(obs.prob = obs.prob,focal = focal)

  if(is.focalList(focal)) {
    if(is.null(scan.number)) {stop("Please provide a `scan.number` if a `focalList` object is passed through `focal`.")}
    focal<- generate_focal(focal.list = focal,scan.number = scan.number)
  }

  sampling.param<- list(
    method = method,
    mode = mode,
    obs.prob = obs.prob,
    focal = focal # also contains focal.list
  )
  class(sampling.param)<- "samplingParam"
  sampling.param
}

#' Determine the sampling parameters `method`
#' from inputted parameters
#'
#' @param obs.prob_fun either:
#' \itemize{
#'   \item{a user-defined function of (i,j,Adj) that output a probability of presence for the dyad,}
#'   \item{a single [0,1] numeric value for all dyad representing their probability of being sampled or not. (obs.prob_type will be "constant")}
#'   \item{the string "random" if each dyad should have its probability drawn from a uniform distribution between 0 and 1 (`runif(n,0,1)`).}
#' }
#'
#' @return a character scalar:
#' \itemize{
#'   \item{"user-defined function"}
#'   \item{"constant"}
#'   \item{or "random"}
#' }
#'
#' @noRd
determine_method<- function(obs.prob,focal){

  switch(method,
         "group" = if(is.null(obs.prob)) {stop("Chosen `method` requires an `obsProb` object.")},
         "focal" = if(is.null(focal)) {stop("Chosen `method` requires a `focal` or `focalList` object")},
         "both" = if(is.null(obs.prob) | is.null(focal)) {stop("Chosen `method` requires both an `obsProb` and a `focalList` objects")},
         stop("Inputted `method` not recognized.")
  )
}
