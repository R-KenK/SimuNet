#' Focal scan from theoretical
#'
#' @param scan.theoretical binary matrix (the theoretical scan within do.scan())
#' @param focal either the row index or the name of the focal
#'
#' @return a focal scan as a binary matrix with matching values for the row of the focal, and NAs otherwise, and the name of the focal as a "focal" attribute.
#' @export
#'
#' @examples
#' # Internal use in do.(non.zero.)scan()
focal.scan<- function(scan.theoretical,focal){
  focal.scan<- scan.theoretical;nodes<- rownames(scan.theoretical)
  if(is.character(focal)){
    focal.index<- match(focal,nodes)
    if(is.na(focal.index)){stop("focal name not recognized.")}
    focal.name<- focal
  }else if(is.numeric(focal)){
    focal.index<- focal
    focal.name<- nodes[focal.index]
  }else{stop("focal format unrecognized")}
  focal.scan[-focal.index,]<- NA
  attr(focal.scan,"focal")<- focal.name;
  focal.scan
}

#' Produce focal.list
#'
#' @param Adj square integers matrix of occurences of dyads.
#' @param total_scan integer, sampling effort. Note that 1/total_scan should be relatively small, increasingly small with increasing precision. Optional if using presence.prob.
#' @param focal.prob_fun a user-defined function of (n,Adj) that output a weight of being focal for each node (passed as `prob` argument to sample() function). By default, pick focals following a uniform distribution. Special case "even" tries to even out the focal.list as much as possible before drawing randomly following a uniform distribution.
#' @param all.sampled logical, should all individuals be sampled before letting them be sampled according to `focal.prob_fun`? Returns an error if total_scan is smaller than the number of nodes.
#'
#' @return a vector of focals (as integers)
#' @export
#'
#' @examples
#' set.seed(42)
#' n<- 6;nodes<- as.character(1:n);
#' total_scan<- 22;n.boot<- 5;
#'
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- round(runif(n*(n-1),0,total_scan*.50))
#'
#' make_focal.list(Adj,total_scan)
#' make_focal.list(Adj,total_scan,focal.prob_fun = "even")
#' make_focal.list(Adj,total_scan,focal.prob_fun = function(n,Adj) 1:n)
#' make_focal.list(Adj,total_scan,focal.prob_fun = function(n,Adj) 1:n*1:n)
#' make_focal.list(Adj,total_scan,focal.prob_fun = function(n,Adj) exp(1:n),all.sampled = TRUE)
#' make_focal.list(Adj,total_scan,focal.prob_fun = function(n,Adj) log(1:n))
#' compute.strength<- function(graph,mode=NULL){
#'   if(is.matrix(graph)){
#'     graph<- igraph::graph.adjacency(graph,mode = mode,weighted = TRUE,add.colnames = TRUE)
#'   }
#'   stren<-igraph::strength(graph)
#'   if(!is.null(names(stren))) {names(stren)<- igraph::vertex_attr(graph)[[1]]}
#'   stren
#' }
#' make_focal.list(Adj,total_scan,focal.prob_fun = function(n,Adj) compute.strength(Adj,"directed"))
make_focal.list<- function(Adj,total_scan,
                           focal.prob_fun = "even",all.sampled = TRUE){

  n<- nrow(Adj);focal.list<- rep(NA,total_scan);

  if(is.character(focal.prob_fun)){
    if(focal.prob_fun=="even"){
      focal.list[is.na(focal.list)]<- c(rep(1:n,total_scan%/%n),sample(1:n,total_scan%%n,replace = FALSE))
      return(focal.list[sample(seq_along(focal.list))])
    }
  }

  if(all.sampled){
    if(n > total_scan){stop("total_scan is too small to sample all nodes.")}
    focal.list[sample(1:total_scan,n)]<- 1:n;total_scan<- total_scan-n;
  }

  if(is.null(focal.prob_fun)){
    focal.list[is.na(focal.list)]<- ceiling(runif(total_scan,0,n))
    return(focal.list)
  }

  P<- focal.prob_fun(n,Adj);if(any(P==0)){P<-P+min(P[P>0])}
  focal.list[is.na(focal.list)]<- sample(1:n,total_scan,replace = TRUE,prob = P)
  focal.list
}
