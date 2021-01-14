# obs.prob generation and related functions -------------------------------

#' Generator for `obsProb` objects
#'
#' @param Adj square integers matrix of occurrences of dyads.
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
#' @param obs.prob_fun either:
#' \itemize{
#'   \item{a user-defined function of (i,j,Adj) that output a probability of presence for the dyad}
#'   \item{a single [0,1] numeric value for all dyad representing their probability of being sampled or not (`obs.prob_type` will be `"constant"`)}
#'   \item{the string `"random"` if each dyad should have its probability drawn from a uniform distribution between 0 and 1 (`runif(n*n,0,1)`)}
#' }
#' @param Adj.subfun subsetting function of the adjacency matrix. Driven by igraph "mode" argument.
#'
#' @return an `obsProb` object (S3 class) containing:
#' \itemize{
#'   \item{`P`: a [0,1] numeric matrix of probability of observation (using the "group" scan sampling `method`) of each dyad.}
#'   \item{`Adj`: inputted `Adj`}
#'   \item{`obs.prob_type`: character scalar, either:
#'     \item{`"constant"`: if all dyad have the same probability of being sampled or not.}
#'     \item{`"random"`: if all dyad have a probability drawn from a uniform distribution between 0 and 1 (`runif(n*n,0,1)`).}
#'     \item{`"user-defined function"`: if the user inputted a function of (i,j,Adj) to calculate each dyad probability.}
#'   }
#'   *
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
#' generate_obsProb(Adj,"directed",obs.prob_fun = "random")
#' generate_obsProb(Adj,"directed",obs.prob_fun = 0.3)
#'
#' # using a user-defined function:
#' user_function.ij<- function(i,j,Adj) {i+j} # comparable to a dyad-trait-based bias
#' user_function.Adj<- function(i,j,Adj) {Adj*Adj} # comparable to a network-based bias
#'
#' generate_obsProb(Adj,"directed",obs.prob_fun = "random")
#' generate_obsProb(Adj,"directed",obs.prob_fun = user_function.ij)
#' generate_obsProb(Adj,"directed",obs.prob_fun = user_function.Adj)
#'
generate_obsProb<- function(Adj,mode,obs.prob_fun = "random",
                             Adj.subfun = NULL){
  if(is.null(Adj.subfun)){
    Adj.subfun<- determine_Adj.subfun(mode = mode)
  }
  n<- nrow(Adj);nodes_names<- rownames(Adj)

  obs.prob_type<- determine_obs.prob_type(obs.prob_fun = obs.prob_fun)

  obs.prob<- list(
    P = calculate_obs.prob(obs.prob_fun = obs.prob_fun,obs.prob_type = obs.prob_type,n = n,nodes_names = nodes_names, Adj = Adj),
    Adj = Adj,
    obs.prob_type = obs.prob_type,
    obs.prob_fun = obs.prob_fun
  )

  class(obs.prob)<- "obsProb"
  obs.prob
}

#' Print method for `obsProb` objects
#' @export
#' @noRd
print.obsProb<- function(x,...){
  print.default(x$P,...)
}

#' Test if object if a `obsProb` object
#'
#' @param x an object to test.
#'
#' @return logical, TRUE if the inputted object is a `obsProb` object.
#'
#' @noRd
is.obsProb<- function(x){
  inherits(x,"obsProb")
}

#' Determine the type of `obsProb` objects inputted or to create.
#' from the class of `obs.prob_fun` inputted.
#'
#' @param obs.prob_fun either:
#' \itemize{
#'   \item{a single [0,1] numeric value for all dyad representing their probability of being sampled or not (`obs.prob_type` will be `"constant"`)}
#'   \item{the string `"random"` if each dyad should have its probability drawn from a uniform distribution between 0 and 1 (`runif(n*n,0,1)`)}
#'   \item{a user-defined function of (i,j,Adj) that output a probability of presence for the dyad}
#' }
#'
#' @return a character scalar:
#' \itemize{
#'   \item{`"user-defined function"`}
#'   \item{`"constant"`}
#'   \item{or `"random"`}
#' }
#'
#' @noRd
determine_obs.prob_type<- function(obs.prob_fun){
  switch(class(obs.prob_fun),
         "numeric" = {
           if(obs.prob_fun<=0 & obs.prob_fun>=1){stop("incompatible numeric obs.prob_fun.")}
           "constant"
         },
         # only considering the case `obs.prob_fun = "random"`
         "character" = obs.prob_fun,
         "function" = "user-defined function",
         stop("`obs.prob_fun` class not recognized.")
  )
}

#' Calculate the dyad probabilities from provided `obs.prob_fun`
#'
#' @param obs.prob_fun either:
#' \itemize{
#'   \item{a single [0,1] numeric value for all dyad representing their probability of being sampled or not (`obs.prob_type` will be `"constant"`)}
#'   \item{the string `"random"` if each dyad should have its probability drawn from a uniform distribution between 0 and 1 (`runif(n*n,0,1)`)}
#'   \item{a user-defined function of (i,j,Adj) that output a probability of presence for the dyad}
#' }
#' @param obs.prob_type a character scalar identified with `determine_obs.prob_type` function:
#' \itemize{
#'   \item{`"user-defined function"`}
#'   \item{`"constant"`}
#'   \item{or `"random"`}
#' }
#' @param n integer, number of node in `Adj`.
#' @param nodes_names character vector or `NULL`, names of the nodes in `Adj`.
#' @param Adj square integers matrix of occurrences of dyads.
#'
#' @importFrom stats runif
#'
#' @return the observation probability matrix P (to be stored in `obs.prob$P`)
#' @noRd
calculate_obs.prob<- function(obs.prob_fun,obs.prob_type,n,nodes_names,Adj){
  P<- switch(obs.prob_type,
             "constant" = matrix(obs.prob_fun,n,n,dimnames = list(nodes_names,nodes_names)),
             # only considering the case `obs.prob_fun = "random"`
             "random" = matrix(stats::runif(n*n,0,1),n,n,dimnames = list(nodes_names,nodes_names)),  # n*n: lazy fix to be sure that all prob are drawn, but likely drawing n probabilities too many presently
             "user-defined function" = {
               dyads<- expand.grid(row = 1:n,col = 1:n)
               P<- matrix(nrow = n,ncol = n,dimnames = list(nodes_names,nodes_names),
                          data =  {
                            P<- sapply(1:nrow(dyads),
                                       function(ij) {
                                         i<- dyads[["row"]];j<- dyads[["col"]];
                                         obs.prob_fun(i,j,Adj)
                                       }
                            )
                            if(any(P<=0)|any(P>=1)){
                              P<- proportional.prob(P)
                            }
                            P
                          }
               )

             },
             stop("`obs.prob_fun` class not recognized.")
  )
  diag(P)<-0
  P
}

#' Modify a vector to be probabilities in ]0,1[
#'
#' @param P a vector of values not all in ]0,1[
#'
#' @return a vector of valid probabilities in ]0,1[
#' @noRd
proportional.prob<- function(P){
  if(any(P<0)){
    P<- P+abs(min(P)) # set the negative minimum to zero
    P<- P+min(P[P>0]) # set the minimum value (zero) to the smallest
  }
  if(any(P==0)){
    P<- P+min(P[P>0]) # set the minimum value (zero) to the smallest
  }
  if(any(P>1)){
    P<- P/(max(P)+min(P)) # set the maximal value to <1
  }
  if(any(P==1)){
    P<- P-min(P)/2 # set the maximal value to 1 - min(P)/2 and min(P) to half the previous min(P)
  }
  P
}

