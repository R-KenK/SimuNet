# #' Perform a scan
# #' According to edge presence probability, or a reference adjacency matrix `Adj` coupled with a sample effort `total_scan`, perform a scan in the form of a theoretical binary adjacency matrix, and simulate the empirical obtention of a group scan and/or focal scan. If the optimization for rare even is used, performs a scan knowing that at least one edge has been observed in the theoretical scan.
# #'
# #' @param Adj square integers matrix of occurences of dyads. Optional if using presence.prob. Update: now if presence prob is passed as Adj (thus all(Adj<1) is TRUE), it will be rightfully assigned to presence.prob. WIP: implement method for association matrices...
# #' @param total_scan integer, sampling effort. Note that 1/total_scan should be relatively small, increasingly small with increasing precision. Optional if using presence.prob.
# #' @param method Character scalar, specifies if the function should return a theoretical perfect group scan, an  empirical group scan (a similarly dimensioned matrix as Adj), or a focal scan (a vector representing the given focal's row in the group scan matrix).
# #' @param ... additional argument to be used, to use produce a scan in a desired way.#'
# #' \itemize{
# #'   \item{obs.prob}{either :
# #'     \item{"a dyad observation obs.probability matrix"}{of same dimension as Adj}
# #'     \item{"a dyad observation vector"}{subsetted similarly as Adj (through the non.diagonal() function for instance)}
# #'     \item{"a general dyad observation obs.probability"}{should be in [0,1], assumed to be the case when only one value is inputed)}
# #'   }
# #'   \item{focal}{Only required for method = "focal" or "both" Character scalar, indicate which focal to consider for the scan.}
# #'   \item{Adj.subfun}{subsetting function of the adjacency matrix. Driven by igraph "mode" argument}
# #'   \item{presence.prob} {square probability matrix of presence (as in Bernouilli trial) of dyads. Optional if using Adj and total_scan.}
# #'   \item{mode} {Character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.}
# #' }
# #' @param use.rare.opti logical: should the optimization for rare event be used?
# #'
# #' @return a list of the theoretical square binary matrix representing the whole group scan, and:
# #' \itemize{
# #'  \item{"method = theoretical"}{that's all,}
# #'  \item{"method = group"}{the empirical group scan square binary matrix}
# #'  \item{"method = focal"}{the empirical focal scan binary vector (and the focal identity as a "focal" attribute)}
# #'  \item{"method = both"}{or the three above}
# #'  }
# #'
# #' @export
# #' @importFrom stats rbinom
# #' @importFrom stats runif
# #'
# #' @examples
# #' set.seed(42)
# #'
# #' n<- 5;nodes<- as.character(1:n);
# #' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
# #' Adj[non.diagonal(Adj)]<- sample(0:42,n*(n-1),replace = TRUE)
# #' Adj
# #'
# #' presence.prob<- Binary.prob(Adj,50)
# #' obs.prob<- matrix(runif(n*n,0,1),n,n);diag(obs.prob)<- 0
# #'
# #' do.scan(Adj,50,method = "group",obs.prob = 0.9)
# #' do.scan(Adj,50,method = "both",obs.prob = obs.prob)
# #' do.scan(presence.prob,method = "focal")
# #' do.scan(presence.prob,method = "focal",focal = "4")
# #' do.scan(presence.prob,method = "focal",focal = 3)
# do.scan<-function(Adj=NULL,total_scan=NULL,
#                   method = c("theoretical","group","focal","both"),...,use.rare.opti=FALSE){
#   # irrelevant bit of code, only to remove annoying note in R CMD Check ----
#   opt.args<- list(...)
#   if(is.null(opt.args$obs.prob)) {obs.prob<- NULL};if(is.null(opt.args$focal)) {focal<- NULL};if(is.null(opt.args$focal.list)) {focal.list<- NULL};if(is.null(opt.args$Adj.subfun)) {Adj.subfun<- NULL};
#   if(is.null(opt.args$presence.prob)) {presence.prob<- NULL};if(is.null(opt.args$mode)) {mode<- NULL};
#   # actual algorithm ----
#   method<- match.arg(method)
#   scan.default.args(Adj,total_scan,method,...)
#
#   # set n and nodes_names according to available inputs
#   if(!is.null(Adj)){n<- nrow(Adj);nodes_names<- rownames(Adj)} else {n<- nrow(presence.prob);nodes_names<- rownames(presence.prob)}
#   # subset a presence probability vector the same way Adj is subset (depends on Adj's mode (cf. igraph))
#   presence.P<- presence.prob[Adj.subfun(presence.prob)];p<- length(presence.P)
#
#   # choose if the optimization for rare events should be performed or not (general case)
#   if(!use.rare.opti){
#     # structure the scan as a matrix filled with zeros
#     scan<- matrix(0,nrow = n,ncol = n,dimnames = list(nodes_names,nodes_names))
#     # core of the randomization: draw a (theoretical) tie or not for each (relevant, cf. triangular matrices or undirected) dyad according to its presence probability
#     # TO FIX: Are each dyad evaluated twice in the undirected (=max) case?
#     scan[Adj.subfun(scan)]<- stats::rbinom(p,1,presence.P)
#   }else{
#     # optimization for rare events: all-zeros scans were already drawn by `simulate_zeros.non.zeros()` in the parent `iterate_scans()`
#     scan<- matrix(0,nrow = n,ncol = n,dimnames = list(nodes_names,nodes_names))
#     # choose a random order between the dyads to draw
#     rand.order<- sample(1:p,p)
#     # determine the cumulative probability that a given dyad will be the first to be drawn with a 1 (tie), after adjusting to such a given conditional probability that there should be at least one tie in the scan.
#     P.cond<- cumsum(adjust.conditional.prob(presence.P[rand.order]))
#     # given the cumulative probability of each dyad (in a random order), determine which would be the first one observed with a tie (at least one will be a tie because of the cumulative distribution of presence ending with a presence probability of 1)
#     first.one<- min(which(stats::runif(1)<P.cond))
#     # set this dyad as a tie
#     scan[Adj.subfun(scan)][rand.order][first.one]<- 1
#     # TO FIX: the dyad before first.one should be zeros
#     # draw the rest of the dyads with their _original_ presence probability (now that the scan did have at least a tie)
#     scan[Adj.subfun(scan)][rand.order][-first.one]<- stats::rbinom(p-1,1,presence.P[rand.order][-first.one])
#   }
#
#   # according to the empirical method chosen, return a scan as is, after masking some dyads with observable_edges(), or after only showing the focal
#   switch(method,
#          "theoretical" = list(theoretical = scan),
#          "group" = list(theoretical = scan,
#                         group = observable_edges(Scan = scan,obs.prob = obs.prob,Adj.subfun = Adj.subfun)
#          ),
#          "focal" = list(theoretical = scan,
#                         focal = focal.scan(scan,focal)
#          ),
#          "both" = list(theoretical = scan,
#                        group = observable_edges(Scan = scan,obs.prob = obs.prob,Adj.subfun = Adj.subfun),
#                        focal = focal.scan(scan,focal)
#          )
#   )
# }

# #' Set arguments to default for do.(non.zero.)scan() function when necessary
# #'
# #' @param Adj square integers matrix of occurences of dyads. Optional if using presence.prob. WIP: implement method for association matrices...
# #' @param total_scan integer, sampling effort. Note that 1/total_scan should be relatively small, increasingly small with increasing precision. Optional if using presence.prob.
# #' @param method Character scalar, specifies if the function should return a theoretical perfect group scan, an  empirical group scan (a similarly dimensioned matrix as Adj), or a focal scan (a vector representing the given focal's row in the group scan matrix).
# #' @param ... additional argument to be used, to use produce a scan in a desired way.#'
# #' \itemize{
# #'   \item{obs.prob}{either :
# #'     \item{"a dyad observation obs.probability matrix"}{of same dimension as Adj}
# #'     \item{"a dyad observation vector"}{subsetted similarly as Adj (through the non.diagonal() function for instance)}
# #'     \item{"a general dyad observation obs.probability"}{should be in [0,1], assumed to be the case when only one value is inputed)}
# #'   }
# #'   \item{focal}{Only required for method = "focal" or "both" Character scalar, indicate which focal to consider for the scan.}
# #'   \item{Adj.subfun}{subsetting function of the adjacency matrix. Driven by igraph "mode" argument}
# #'   \item{presence.prob} {square probability matrix of presence (as in Bernouilli trial) of dyads. Optional if using Adj and total_scan.}
# #'   \item{mode} {Character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.}
# #' }
# #'
# #' @return nothing, but assign required variable for do.(non.zero.)scan() in their environment.
# #'
# #' @export
# #'
# #' @examples
# #' # Internal use.
# scan.default.args<- function(Adj,total_scan,method,...){
#   opt.args<- list(...)
#   # irrelevant bit of code, only to remove annoying note in R CMD Check ----
#   if(is.null(opt.args$obs.prob)) {obs.prob<- NULL};if(is.null(opt.args$focal)) {focal<- NULL};if(is.null(opt.args$Adj.subfun)) {Adj.subfun<- NULL};
#   if(is.null(opt.args$presence.prob)) {presence.prob<- NULL};if(is.null(opt.args$mode)) {mode<- NULL};
#   # actual algorithm ----
#
#   lapply(names(opt.args),function(name) assign(name,opt.args[[name]],parent.frame(n = 2)))
#
#   if(is.null(opt.args$obs.prob)){assign("obs.prob",1,parent.frame())}else{assign("obs.prob",opt.args$obs.prob,parent.frame(n = 1))}
#
#   if(is.null(opt.args$mode)){assign("mode","directed",parent.frame());opt.args$mode<- "directed"}else{assign("mode",opt.args$mode,parent.frame(n = 1))}
#
#   if(is.null(opt.args$Adj.subfun)){
#     assign(
#       "Adj.subfun",
#       switch(opt.args$mode,
#              "directed" = ,
#              "undirected" = ,
#              "max" = ,
#              "min" = ,
#              "plus" = non.diagonal,
#              "upper" = upper.tri,
#              "lower" =  lower.tri,
#              "vector" = function(V){rep(TRUE,length(V))}
#       ),
#       parent.frame(n = 1)
#     )
#   }else{assign("Adj.subfun",opt.args$Adj.subfun,parent.frame(n = 1))}
#
#   if(is.null(opt.args$presence.prob)){
#     if(all(Adj<1&Adj>=0)){
#       assign("presence.prob",Adj,parent.frame(n = 1)) # assuming presence prob has been passed as first argument `Adj`
#     }
#     else{
#       if(is.null(total_scan)){
#         stop("Please provide either `total_scan` and `Adj`, or `presence.prob` to run a scan.")
#       }
#       assign("presence.prob",Binary.prob(Adj,total_scan,mode = opt.args$mode),parent.frame(n = 1))
#     }
#   }else{assign("presence.prob",opt.args$presence.prob,parent.frame(n = 1))}
#
#   if(is.null(opt.args$focal)){
#     if(method!="group"){
#       if(!is.null(Adj)){
#         assign("focal",quick.sample(1:nrow(Adj),1),parent.frame(n = 1))
#       }else{
#         assign("focal",quick.sample(1:nrow(opt.args$presence.prob),1),parent.frame(n = 1))
#       }
#     }else{
#       assign("focal",NULL,parent.frame(n = 1))
#     }
#   }else{assign("focal",opt.args$focal,parent.frame(n = 1))}
#
#   if(is.null(opt.args$focal.list)&!is.null(total_scan)){
#     if(method!="group"){
#       if(!is.null(Adj)){
#         assign("focal.list","even",parent.frame(n = 1))
#       }else{
#         assign("focal.list","even",parent.frame(n = 1))
#       }
#     }else{
#       assign("focal.list",NULL,parent.frame(n = 1))
#     }
#   }else{assign("focal.list",opt.args$focal.list,parent.frame(n = 1))}
# }

# Adjacency mode tools ----------------------------------------------------

# #' Make Adjacency fit the selected mode
# #' Internal use.
# #'
# #' @param raw.scan a raw (directed by default) binary adjacency matrix, pre-theoretical or pre-empirical, to which the chosen `mode` is to be applied
# #' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
# #'
# #' @return an adjacency matrix fitting the chosen `mode` is to be applied
# #' @noRd
# apply_mode<- function(raw.scan,mode = c("directed", "undirected", "max","min", "upper", "lower", "plus","vector"),na.minimize = FALSE){
#   switch(mode,
#          "undirected" = , # as in `igraph`, consider this mode to be the same as `max`
#          "max" = {
#            min.value<- min(raw.scan,na.rm = TRUE);replace_NA_by_half.min<- function(X) {replace_NA(X = X,value = min.value/2)} # sets `NA`s to a value between zero and the min, in order to turn the directed into undirected with only one `>=` test.
#            both.na<- is.na(raw.scan) & is.na(t(raw.scan))
#            ifelse(  # `ifelse` rather than `if()` statements to return a vector or matrix after inputting a logical vector or matrix
#              test = !both.na,
#              yes = ifelse(
#                compare_with_transposed(raw.scan,comp.fun = `>=`,manage_NA.fun = replace_NA_by_half.min),
#                raw.scan,
#                t(raw.scan)
#              ),
#              no = NA # double NAs
#            )
#          },
#          "min" = {
#            min.value<- min(raw.scan,na.rm = TRUE);replace_NA_by_half.min<- function(X) {replace_NA(X,value = min.value/2)} # sets `NA`s to a value between zero and the min, in order to turn the directed into undirected with only one `<=` test.
#            both.na<- is.na(raw.scan) & is.na(t(raw.scan))
#            ifelse(
#              test = !both.na,
#              yes = ifelse(
#                compare_with_transposed(raw.scan,comp.fun = `<=`,manage_NA.fun = replace_NA_by_half.min),
#                raw.scan,
#                t(raw.scan)
#              ),
#              no = NA
#            )
#          },
#          "plus" = {  # WHAT DOES THIS MEAN FOR BINARY SCANS?
#            not.na<- !is.na(raw.scan) & !is.na(t(raw.scan))
#            ifelse(not.na,raw.scan+t(raw.scan),NA)
#          },
#          "directed" = ,
#          "upper" = ,
#          "lower" =  ,
#          "vector" = raw.scan
#   )
# }
