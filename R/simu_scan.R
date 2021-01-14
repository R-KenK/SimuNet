# WIP: Wrapper to simulate a single scan ------------------------------------------

# #' Simulate a single scan, theoretical or empirical
# #'
# #' @param Adj square integers matrix of occurrences of dyads. Can be used to pass a `presenceProb` object too.
# #' @param total_scan integer, sampling effort. Note that 1/total_scan should be relatively small, increasingly small with increasing precision. Optional if using a `presenceProb` object.
# #' @param sampling.param Optional if a theoretical scan is needed. Otherwise `samplingParam` object that stores:
# #' \itemize{
# #'   \item{method}{a character scalar, either:
# #'     \item{"group": to use the group-scan sampling method}
# #'     \item{"focal": to use the group-scan sampling method}
# #'     \item{"both": to use both}
# #'   }
# #'   \item{obs.prob: an `obsProb` object.}
# #'   \item{focal: a `focal` object.}
# #' }
# #' @param method Character scalar, specifies if the function should return a theoretical perfect group scan, an  empirical group scan (a similarly dimensioned matrix as Adj), or a focal scan (a vector representing the given focal's row in the group scan matrix). # TO REWRITE WITH ITEMS
# #' @param obs.prob_fun Optional if a theoretical scan is needed.  either:
# #' \itemize{
# #'   \item{a user-defined function of (i,j,Adj) that output a probability of presence for the dyad,}
# #'   \item{a single [0,1] numeric value for all dyad representing their probability of being sampled or not. (obs.prob_type will be "constant")}
# #'   \item{the string "random" if each dyad should have its probability drawn from a uniform distribution between 0 and 1 (`runif(n,0,1)`).}
# #' }
# #' @param focal.prob_fun either:
# #' \itemize{
# #'   \item{a user-defined function of (n,Adj) that output a weight of being focal for each node (passed as `prob` argument to `base::sample` function)}
# #'   \item{`NULL` or `"random"`, pick focals following a uniform distribution}
# #'   \item{Special case `"even"` tries to even out the `focal.list` as much as possible before drawing randomly following a uniform distribution}
# #' }
# #' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
# #'
# #' @return
# #' @export
# #'
# #' @examples
# #' set.seed(42)
# #'
# #' n<- 5;nodes<- letters[1:n];total_scan<- 42;
# #' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
# #' Adj[non.diagonal(Adj)]<- sample(0:total_scan,n*(n-1),replace = TRUE)
# #' Adj
# #'
# #' # by default will simulate a directed theoretical scan
# #' simu_scan(Adj,total_scan)
# #'
# #' # Users can generate sampling parameters to use in simu_scan
# #' para.group.constant<- generate_samplingParam(method = "group",mode = "max",obs.prob = 0.42)
# #' simu_scan(Adj,total_scan,para.group.constant)
# #'
# #' # using a user-defined function:
# #' user_function.ij<- function(i,j,Adj) {i+j} # comparable to a dyad-trait-based bias
# #' user_function.Adj<- function(i,j,Adj) {Adj*Adj} # comparable to a network-based bias
# #'
# #' generate_obsProb(Adj,"directed",obs.prob_fun = user_function.ij)
# #' generate_obsProb(Adj,"directed",obs.prob_fun = user_function.Adj)
# simu_scan<-function(Adj = NULL,total_scan = NULL,sampling.param = NULL,method = NULL,obs.prob_fun = NULL,focal.prob_fun = NULL,
#                     mode = c("directed", "undirected", "max","min", "upper", "lower", "plus","vector")){
#   mode<- match.arg(mode)
#
#   # Check if `Adj` has been passed as a `presenceProb` object or not, retrieve other variable otherwise if needed
#   if (!is.presenceProb(Adj)) {
#     Adj.subfun<- determine_Adj.subfun(mode)
#     presence.prob<- generate_presenceProb(Adj = Adj,total_scan = total_scan,mode = mode,Adj.subfun = Adj.subfun)
#   } else {
#     presence.prob<- Adj
#   }
#
#   # use the class "scan" object generator to
#   scan<- generate_scan(presence.prob)
#
#   # output either the theoretical scan or applies sample_from_scan
#   if (is.null(sampling.param) & is.null(method) & is.null(obs.prob_fun) & is.null(focal.prob_fun)) {
#     scan
#   } else {
#     if (is.null(sampling.param)){
#       sampling.param<- generate_samplingParam(method = method,mode = mode,obs.prob = obs.prob,focal.list = focal.list) # should return an error in case of missing parameters given the chosen method
#     }
#     generate_empiScan(scan = scan,sampling.param = sampling.param)
#   }
# }
