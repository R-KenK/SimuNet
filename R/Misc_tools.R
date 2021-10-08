#' Quick optimized equivalent to sample(x,size,replace=TRUE)
#'
#' @param x a vector
#' @param size number of elements to sample
#'
#' @importFrom stats runif
#'
#' @return sample of size "size" taken from x
#' @export
#'
#' @examples
#' quick_sample(1:20,5)
quick_sample<- function(x,size){
  x[ceiling(stats::runif(size,0,length(x)))]
}

# Make cluster object wrapper
#
# @param n.cores number of threads to use.
# @param .export vector of variable/function names to use in each clusters
#
# @return a cluster object
# @export
#
# @importFrom  snow makeCluster
# @importFrom  snow clusterExport
# @importFrom  doSNOW registerDoSNOW
#
# @examples
# # Internal use
# make_cl<- function(n.cores,.export){
#   cl<- snow::makeCluster(n.cores);snow::clusterExport(cl,list = .export);
#   cl
# }
