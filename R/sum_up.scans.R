#' Sum up and flatten list of binary scans
#' Internal use. Sum up binary adjacency matrices or focal vectors according to the method used. sum_up.scan*s*() is a wrapper for sum_up.scan() that returns summed_up scans for the differnet requested methods0
#'
#' @param scan_list list of binary adjacency matrices, binary focal vectors, or both, outputed by iterating do.scan()
#' @param scaled logical, specifies if adjacency data should be scaled by sampling effort.
#' @param use.rare.opti logical: should the optimization for rare event be used?
#' @param obs.prob either :
#' \itemize{
#'  \item{"a dyad observation probability matrix (P.obs)"}{of same dimension as Adj}
#'  \item{"a dyad observation vector"}{subsetted similarly as Adj (through the non.diagonal() function for instance)}
#'  \item{"a systematic dyad observation (P.obs constant for all i,j)"}{should be in [0,1], assumed to be the case when only one value is inputed)}
#' }
#' @param method Character scalar, specify if the function should use a whole group or a focal scan sampling method (or both).
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
#' @param focal.list Character vector, indicate the list of focals to consider throughout the scans.
#'
#' @return An integer adjacency matrix, or a list of two if method = "both". Now returns NAs for dyad scaled with no observation at all.
#' @export
#'
#' @examples
#' #Internal use for readability
sum_up.scans<- function(scan_list,scaled=FALSE,use.rare.opti=FALSE,obs.prob = NULL,focal.list = NULL,
                        method = c("theoretical","group","focal","both"),mode = c("directed", "undirected", "max","min", "upper", "lower", "plus","vector")){
  method<- match.arg(method);
  mode<- match.arg(mode);
  switch(method,
         "theoretical" = {
           list(
             theoretical = sum_up.scan(scan_list,"theoretical",mode = mode,scaled = scaled,use.rare.opti = use.rare.opti,obs.prob = obs.prob)
           )
         },
         "group" = {
           list(
             theoretical = sum_up.scan(scan_list,"theoretical",mode = mode,scaled = scaled,use.rare.opti = use.rare.opti,obs.prob = obs.prob),
             group = sum_up.scan(scan_list,"group",mode = mode,scaled = scaled,use.rare.opti = use.rare.opti,obs.prob = obs.prob)
           )
         },
         "focal" = {
           list(
             theoretical = sum_up.scan(scan_list,"theoretical",mode = mode,scaled = scaled,use.rare.opti = use.rare.opti,obs.prob = obs.prob),
             focal = sum_up.scan(scan_list,"focal",mode = mode,scaled = scaled,use.rare.opti = use.rare.opti,focal.list =  focal.list)
           )
         },
         "both" = {
           list(
             theoretical = sum_up.scan(scan_list,"theoretical",mode = mode,scaled = scaled,use.rare.opti = use.rare.opti,obs.prob = obs.prob),
             group = sum_up.scan(scan_list,"group",mode = mode,scaled = scaled,use.rare.opti = use.rare.opti,obs.prob = obs.prob),
             focal = sum_up.scan(scan_list,"focal",mode = mode,scaled = scaled,use.rare.opti = use.rare.opti,focal.list =  focal.list)
           )
         }
  )
}

#' Sum up and flatten list of binary scans
#' Internal use. Sum up binary adjacency matrices or focal vectors according to the method used. sum_up.scan*s*() is a wrapper for sum_up.scan() that returns summed_up scans for the differnet requested methods0
#'
#' @param scan_list list of binary adjacency matrices, binary focal vectors, or both, outputed by iterating do.scan()
#' @param method Character scalar, specify if the function should use a whole group or a focal scan sampling method (or both).
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
#' @param scaled logical, specifies if adjacency data should be scaled by sampling effort.
#' @param use.rare.opti logical: should the optimization for rare event be used?
#' @param obs.prob either :
#' \itemize{
#'  \item{"a dyad observation probability matrix (P.obs)"}{of same dimension as Adj}
#'  \item{"a dyad observation vector"}{subsetted similarly as Adj (through the non.diagonal() function for instance)}
#'  \item{"a systematic dyad observation (P.obs constant for all i,j)"}{should be in [0,1], assumed to be the case when only one value is inputed)}
#' }
#' @param focal.list Character vector, indicate the list of focals to consider throughout the scans.
#'
#' @return An integer adjacency matrix, or a list of two if method = "both". Now returns NAs for dyad scaled with no observation at all.
#' @export
#'
#' @examples
#' #Internal use for readability
sum_up.scan<- function(scan_list,method,mode = NULL,scaled,use.rare.opti = FALSE,obs.prob = NULL,focal.list = NULL){
 if(use.rare.opti) {n.zeros<- attr(scan_list,"n.zeros")} else {n.zeros<- 0}

  # isolate desired scan method and apply the adjacency mode
  if(mode=="min"&method=="focal"){mode<- "directed"}
  scan_list.method<- lapply(scan_list,function(scan) adjacency_mode(scan[[method]],mode = mode));
  summed_up.method<- Reduce(matrix_sum_na.rm,scan_list.method) # sum the scan list, considering NAs as zeros
  if(mode=="min"&method=="focal"){summed_up.method<- adjacency_mode(summed_up.method,mode = "min")}


  # determine sampling effort
  if(method=="theoretical"){ # shortcut for the theoretical method
    n<- nrow(summed_up.method);if(!is.null(row.names(summed_up.method))){nodes<- row.names(summed_up.method)}
    total_scan<- (length(scan_list)+n.zeros)*ifelse(mode=="plus",2,1)
    n.observed<- matrix(rep(total_scan,n*n),n,n,dimnames = list(nodes,nodes))
  }else{ # otherwise count non-NAs (considering the right mode)
    n.observed<- n.observed_edges(scan_list.method,mode = mode,diag = 1,use.rare.opti = use.rare.opti,obs.prob = obs.prob,focal.list = focal.list,n.zeros = n.zeros,method = method) # here diag = 1 because while the count of the diagonal is irrelevant, it shouldn't be 0/0.
    n.observed<- ifelse(n.observed!=0,n.observed,NA) # will return an NA instead of 0 to avoid divide by zero later
  }

  # if summed-up weighted adjacency has to be scaled by effective sampling effort
  if(scaled){summed_up.method<- summed_up.method/n.observed}

  attr(summed_up.method,"observed_edges")<- n.observed_edges(scan_list.method,mode = mode,diag = 0,use.rare.opti = use.rare.opti,obs.prob = obs.prob,focal.list = focal.list,n.zeros = n.zeros,method = method)
  summed_up.method
}

#' Get number of edge observations (for group scans with unobservable individuals)
#' quantify actual edge-wise sampling effort, considering that some weren't observable in all group scans.Internal use.
#'
#' @param scan_list list of binary group scans, with NAs when the dyad was not observable.
#' @param diag integer (mostly), value to replace the diagonal of the output matrix with. Use NULL if you consider self-loop (untested).
#' @param use.rare.opti logical: should the optimization for rare event be used?
#' @param obs.prob either :
#' \itemize{
#'  \item{"a dyad observation probability matrix (P.obs)"}{of same dimension as Adj}
#'  \item{"a dyad observation vector"}{subsetted similarly as Adj (through the non.diagonal() function for instance)}
#'  \item{"a systematic dyad observation (P.obs constant for all i,j)"}{should be in [0,1], assumed to be the case when only one value is inputed)}
#' }
#' @param n.zeros integer, the attribute outputed by `simulate_zeros.non.zeros`, representing the number of full-zero scans. Used only when use.rare.opti=TRUE
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
#' @param focal.list Character vector, indicate the list of focals to consider throughout the scans.
#' @param method Character scalar, specify if the function should use a whole group or a focal scan sampling method (or both).
#'
#' @importFrom stats rbinom
#'
#' @return a square matrix with element quantifying how many time a dyad has been sampled
#' @export
#'
#' @examples
#' #internal use.
n.observed_edges<- function(scan_list,mode,diag=NULL,use.rare.opti=FALSE,obs.prob=NULL,focal.list = NULL,n.zeros = 0,method=NULL){
  n<- nrow(scan_list[[1]]);
  if(!is.null(rownames(scan_list[[1]]))){
    nodes<- rownames(scan_list[[1]])
  }else{
    nodes<- as.character(1:n)
  }

  switch(method,
         "theoretical" = {
           n.observed<- adjacency_mode(matrix(length(scan_list)+n.zeros,n,n),mode = mode)
         },
         "focal" = {
           if(mode!="directed"){mode<- "plus"}
           focal.table<- table(focal.list);names(focal.table)<- nodes
           n.observed<- adjacency_mode(
             rbind_lapply(nodes,function(node) rep(focal.table[node],n)),
             mode = mode
           )
           colnames(n.observed)<- rownames(n.observed)<- nodes
           n.observed
         },
         "group" = {
           n.observed<- Reduce("+",
                               lapply(scan_list,
                                      function(scan) ifelse(!is.na(adjacency_mode(scan,mode = mode)),1,0) # counting part of the algorithm
                               )
           )
           if(use.rare.opti){
             if(!is.matrix(obs.prob)){obs.prob<- matrix(obs.prob,n,n)}
             n.observed_zero.scans<- rbind_lapply(1:n,function(i) rbinom(n,n.zeros,obs.prob[i,]))
             n.observed<- n.observed + adjacency_mode(n.observed_zero.scans,mode = mode)
           }
         }
  )
  if(!is.null(diag)) {diag(n.observed)<- diag} # doesn't count the diagonal by default. Left the option to count if self loops should be considered
  rownames(n.observed)<- colnames(n.observed)<- nodes
  n.observed
}
