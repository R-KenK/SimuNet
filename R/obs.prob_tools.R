# Hiding edges during group scan ------------------------------------------

#' Hide unobservable edges
#' Simulate that some dyads might reasonable not be observable during a group scan
#'
#' @param Scan a square matrix
#' @param obs.prob either :
#' \itemize{
#'  \item{"a dyad observation obs.probability matrix"}{of same dimension as Scan}
#'  \item{"a dyad observation vector"}{subsetted similarly as Scan (through the non.diagonal() function for instance)}
#'  \item{"a systematic dyad observation obs.probability"}{should be in [0,1], assumed to be the case when only one value is inputed)}
#' }
#' @param Adj.subfun subsetting function of the adjacency matrix. Driven by igraph "mode" argument
#'
#' @return a similar matrix as Scan where some non diagonal edges have a obs.probability to be NAs.
#' @export
#'
#' @examples
#' set.seed(42)
#'
#' n<- 6;nodes<- as.character(1:n);
#' Scan<- matrix(rbinom(n*n,1,0.2),n,n,dimnames = list(nodes,nodes));diag(Scan)<- 0;
#' obs.prob<- matrix(runif(n*n,0,1),n,n);diag(obs.prob)<- 0
#'
#' observable_edges(Scan,obs.prob,non.diagonal)
observable_edges<- function(Scan,obs.prob=NULL,Adj.subfun=NULL){
  observed<- Scan
  if(!is.null(obs.prob)){
    if(length(obs.prob)==1) {
      if(obs.prob<=1 & obs.prob>=0){
        obs.prob<- rep(obs.prob,length(Scan[Adj.subfun(Scan)]))
      }else{
        stop("Single observation obs.probability provided should be within [0,1]")
      }
    }else{
      if(is.matrix(obs.prob)) {
        obs.prob<- obs.prob[Adj.subfun(obs.prob)]
      }
      if(length(obs.prob)!=length(Scan[Adj.subfun(Scan)])){
        stop("Matrix or vector obs.prob dimension(s) incompatible with adjacency matrix's")
      }
    }
    missed<- rbinom(length(obs.prob),1,obs.prob)==0
    observed[Adj.subfun(observed)][missed]<- NA
  }
  observed
}

# obs.prob tools ----------------------------------------------------------

#' Produce matrix of probability of observation from user-defined function
#'
#' @param Adj square integers matrix of occurences of dyads.
#' @param obs.prob_fun either a user-defined function of (i,j,Adj) that output a probability of presence for the dyad, or a single value to indicate a constant observation probability
#' @param Adj.subfun subsetting function of the adjacency matrix. Default is non.diagonal.
#'
#' @return a matrix of probability of observation for each dyad (obs.prob)
#' @export
#'
#' @examples
#' set.seed(42)
#' n<- 6;nodes<- as.character(1:n);
#' total_scan<- 20;n.boot<- 5;
#' focal.list<- sample(nodes,total_scan,replace = TRUE)
#'
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- round(runif(n*(n-1),0,total_scan*.50))
#'
#' make_obs.prob(Adj)
#' make_obs.prob(Adj,obs.prob_fun = 0.2)
#' make_obs.prob(Adj,obs.prob_fun = function(i,j,Adj) i+j)
#' compute.EV<- function(graph,mode=NULL){
#'   if(is.matrix(graph)){
#'     graph<- igraph::graph.adjacency(graph,mode = mode,weighted = TRUE,add.colnames = TRUE)
#'   }
#'   EV<- igraph::eigen_centrality(graph, weights = igraph::E(graph)$weight,scale = FALSE)$vector
#'   if(!is.null(names(EV))){names(EV)<- igraph::vertex_attr(graph)[[1]]}
#'   EV
#' }
#' make_obs.prob(Adj,obs.prob_fun = function(i,j,Adj){EVs<- compute.EV(Adj,"directed");EVs[i]*EVs[j]})
make_obs.prob<- function(Adj,obs.prob_fun = NULL,
                         Adj.subfun = non.diagonal){
  if(is.numeric(obs.prob_fun)){
    if(length(obs.prob_fun)==1 & obs.prob_fun>0 & obs.prob_fun<1){
      return(obs.prob_fun)
    }else{
      stop("incompatible numeric obs.prob_fun.")
    }
  }

  n<- nrow(Adj);if(!is.null(rownames(Adj))){nodes<- rownames(Adj)}else{nodes<- as.character(1:n)}

  if(is.null(obs.prob_fun)){
    obs.prob<- matrix(runif(n,0,1),n,n,dimnames = list(nodes,nodes))
    diag(obs.prob)<- 0
    return(obs.prob)
  }

  dyads<- expand.grid(row = 1:n,col = 1:n)
  obs.prob<- matrix(nrow = n,ncol = n,dimnames = list(nodes,nodes),
                    data =  sapply(1:nrow(dyads),
                                   function(ij) {
                                     i<- dyads[["row"]];j<- dyads[["col"]];
                                     obs.prob_fun(i,j,Adj)
                                   }
                    )
  )
  diag(obs.prob)<- 0;P<- obs.prob[Adj.subfun(obs.prob)]
  if(any(P<=0)|any(P>=1)){
    obs.prob[Adj.subfun(obs.prob)]<- proportional.prob(P)
  }
  obs.prob
}
