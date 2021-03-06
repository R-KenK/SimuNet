# #' Simulate networks
# #' Simulate iterated binary group or focal scans with probabilities derived from an empirical adjacency matrix.
# #'
# #' @param Adj square integers matrix of occurrences of dyads. WIP: implement method for association matrices...
# #' @param n.simu integer, number of simulation to perform.
# #' @param total_scan integer, sampling effort. Note that 1/total_scan should be relatively small, increasingly small with increasing precision.
# #' @param method Character scalar, specifies if the function should return a theoretical perfect group scan, an  empirical group scan (a similarly dimensioned matrix as Adj), or a focal scan (a vector representing the given focal's row in the group scan matrix).
# #' @param focal.list Character vector, indicate the list of focals to consider throughout the scans. Can also be carried to `make_focal.list` as the `focal.prob_fun`argument (i.e. either "even" or a user-defined function).
# #' @param scaled logical, specifies if adjacency data should be scaled by sampling effort.
# #' @param cl Optional cluster object (cf snow package), experimentally set to put the makeCluster and stopCluster out of the bootable function. (WIP, next implementation should rethink this).
# #' @param ... additional argument to be used, to use produce a scan in a desired way.
# #' @param obs.prob either:
# #' \itemize{
# #'   \item{"a dyad observation obs.probability matrix"}{of same dimension as Adj}
# #'   \item{"a dyad observation vector"}{subsetted similarly as Adj (through the non.diagonal() function for instance)}
# #'   \item{"a general dyad observation obs.probability"}{should be in [0,1], assumed to be the case when only one value is inputed)}
# #'   \item{"a user-defined function to be passed to `make_obs.prob` as `obs.prob_fun`"}{should be a function of (i,j,Adj))}
# #' }
# #' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
# #' @param output Character scalar, specify if the function should return the list of scans, or reduce them into the bootstrapped adjacency matrix
# #' @param use.rare.opti logical: should the optimization for rare event be used? If left NULL, choice is made automatically by decide_use.rare.opti().
# #' @param cl Optional cluster object (cf snow package), left in case it's actually faster. So far benchmarks showed that with the overhead using multiple thread was actually longer...
# #'
# #' @return according to output and method: a list of iterated scans, or of adjacency matrix, with attributes to keep track of certain data
# #'
# #' @details from input adjacency matrix `Adj`, will iterate `total_scan` theoretical scans (through `iterate_fun`), alter them given an empirical sampling `method`, and sum them up (through `sum_up_fun`) into a weighted adjacency matrix. This whole process will be repeated `n.simu` times.
# #'
# #' @export
# #' @importFrom parallel detectCores
# #' @importFrom snow makeCluster
# #' @importFrom snow stopCluster
# #' @importFrom doSNOW registerDoSNOW
# #' @importFrom pbapply pblapply
# #'
# #' @examples
# #' set.seed(42)
# #'
# #' n<- 5;nodes<- as.character(1:n);
# #' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
# #' Adj[non.diagonal(Adj)]<- sample(0:42,n*(n-1),replace = TRUE)
# #' Adj
# #'
# #' obs.prob<- matrix(runif(n*n,0,1),n,n);diag(obs.prob)<- 0
# #' focal.list<- sample(nodes,42,replace = TRUE)
# #' table(focal.list)
# #'
# #' Simu_net(Adj,total_scan = 42,focal.list = focal.list,n.simu = 3,scaled = TRUE,
# #'            method = "group",use.rare.opti=FALSE,mode = "directed",obs.prob = 0.5,output = "list")
# #' Simu_net(Adj,total_scan = 42,focal.list = focal.list,n.simu = 3,scaled = FALSE,
# #'            method = "focal",use.rare.opti=FALSE,mode = "max",output = "adj")
# #' Simu_net(Adj,total_scan = 42,focal.list = "even",n.simu = 3,scaled = TRUE,
# #'            method = "focal",use.rare.opti=FALSE,mode = "plus",output = "adj")
# #' Simu_net(Adj,total_scan = 42,focal.list = function(n,Adj) 1:n*1:n,n.simu = 3,scaled = TRUE,
# #'            method = "focal",use.rare.opti=FALSE,mode = "max",output = "adj")
# #' Simu_net(Adj,total_scan = 42,obs.prob = 0.2,n.simu=3,scaled = TRUE,
# #'            method = "group",use.rare.opti=FALSE,mode = "directed",output = "list")
# #' Simu_net(Adj,total_scan = 42,obs.prob = obs.prob,n.simu=3,scaled = TRUE,
# #'            method = "both",use.rare.opti=FALSE,mode = "directed",output = "all")
# #' Simu_net(Adj,total_scan = 400,obs.prob = obs.prob,n.simu=3,scaled = TRUE,
# #'            method = "both",mode = "directed",use.rare.opti = TRUE,output = "adj")
# Simu_net<- function(Adj,total_scan,method=c("theoretical","group","focal","both"),focal.list=NULL,n.simu,...,
# scaled=FALSE,obs.prob = NULL,
#                       mode = c("directed", "undirected", "max","min", "upper", "lower", "plus","vector"),
#                       output=c("list","adjacency","all"),use.rare.opti=FALSE,cl=NULL){
#   #irrelevant bit of code, only to remove annoying note in R CMD Check...
#   # irrelevant bit of code, only to remove annoying note in R CMD Check ----
#   b<-NULL;opt.args<- list(...);if(is.null(opt.args$Adj.subfun)) {Adj.subfun<- NULL};if(is.null(opt.args$presence.prob)) {presence.prob<- NULL};
#
#   # actual algorithm ----
#   method<- match.arg(method)
#   output<- match.arg(output)
#   mode<- match.arg(mode)
#
#   # high-level check if the general method or the optimization for rare event should be used
#   Sim_fun<- ifelse(!use.rare.opti,Sim_oneNet,Sim_oneNet_rareOpti)
#
#   scan.default.args(Adj = Adj,total_scan = total_scan,method = method,obs.prob = obs.prob,mode = mode,focal.list = focal.list,...)
#
#   Bootstrap<- pbapply::pblapply(
#     1:n.simu,
#     function(sim){
#       Sim_fun(Adj = Adj,total_scan = total_scan,presence.prob = presence.prob,
#                     focal.list = focal.list,scaled = scaled,obs.prob = obs.prob,
#                     method = method,mode = mode,output = output,Adj.subfun = Adj.subfun)
#     },cl = cl
#   )
#   Bootstrap_add.attributes(Bootstrap = Bootstrap,method = method,scaled = scaled,use.rare.opti = use.rare.opti,
#                            mode = mode,output = output,total_scan = total_scan,n.simu = n.simu,presence.prob = presence.prob,obs.prob = obs.prob,focal.list = focal.list)
# }


# if(is.null(use.rare.opti)){
#   n<- nrow(Adj)
#   use.rare.opti<- decide_use.rare.opti(n = n,total_scan = total_scan,max.obs = max(Adj),alpha = 0.05)
# }
