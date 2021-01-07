#' Row bind list of data frames
#' wrapper to one-function do.call rbind over a lapply list
#'
#' @param X a list. See details \link[base]{lapply}.
#' @param FUN a function to subset data frames (or data tables). See details \link[base]{lapply}.
#'
#' @return a row bound data frame
#' @export
#'
#' @examples
#' set.seed(42)
#'
#' X<- lapply(1:3,function(i) list(int = 42,df = data.frame(x = runif(10,0,1),y = runif(10,0,1))))
#' rbind_lapply(X,function(x) x$df)
rbind_lapply<- function(X,FUN){
  do.call(rbind,lapply(X = X,FUN = FUN))
}

#' Column bind list of data frames
#' wrapper to one-function do.call cbind over a lapply list
#'
#' @param X a list. See details \link[base]{lapply}.
#' @param FUN a function to subset data frames (or data tables). See details \link[base]{lapply}.
#'
#' @return a column bound data frame
#' @export
#'
#' @examples
#' set.seed(42)
#'
#' X<- lapply(1:3,function(i) list(int = 42,df = data.frame(x = runif(10,0,1),y = runif(10,0,1))))
#' cbind_lapply(X,function(x) x$df)
cbind_lapply<- function(X,FUN){
  do.call(cbind,lapply(X = X,FUN = FUN))
}

#' Row bind list of data frames
#' wrapper to one-function do.call rbind over a pblapply list
#'
#' @param X a list. See details \link[base]{lapply}.
#' @param FUN a function to subset data frames (or data tables). See details \link[base]{lapply}.
#' @param n.cores number or cores (threads rather) to use in the cluster.
#' @param .export character vector of the name of the variable and function to export to each worker of the cluster
#' @param cl cluster object
#'
#' @return a row bound data frame
#' @export
#'
#' @examples
#' set.seed(42)
#'
#' X<- lapply(1:3,function(i) list(int = 42,df = data.frame(x = runif(10,0,1),y = runif(10,0,1))))
#' rbind_lapply(X,function(x) x$df)
rbind_pblapply<- function(X,FUN,n.cores=NULL,.export=NULL,cl=NULL){
  if(is.null(cl)){cl<- make_cl(n.cores,.export);on.exit(snow::stopCluster(cl))}
  do.call(rbind,pbapply::pblapply(X = X,FUN = FUN,cl = cl))
}

#' Make cluster object wrapper
#'
#' @param n.cores number of threads to use.
#' @param .export vector of variable/function names to use in each clusters
#'
#' @return a cluster object
#' @export
#'
#' @importFrom  snow makeCluster
#' @importFrom  snow clusterExport
#' @importFrom  doSNOW registerDoSNOW
#'
#' @examples
#' # Internal use
make_cl<- function(n.cores,.export){
  cl<- snow::makeCluster(n.cores);snow::clusterExport(cl,list = .export);
  cl
}

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
#' quick.sample(1:20,5)
#' # microbenchmark::microbenchmark(runif={(1:20)[ceiling(runif(5,0,20))]},
#' #   quick.sample=quick.sample(1:20,5),sample=sample(1:20,5,replace = TRUE),
#' #   times = 1000,control = list("warmup"=100))
quick.sample<- function(x,size){
  x[ceiling(stats::runif(size,0,length(x)))]
}

#' Two-sample t-tests from sample statistics
#' cf. https://stats.stackexchange.com/a/30450/255116
#'
#' @param m1 numeric, sample mean in population 1.
#' @param m2 numeric, sample mean in population 2.
#' @param sd1 numeric, sample standard deviation in population 1.
#' @param sd2 numeric, sample standard deviation in population 2.
#' @param n1 integer, sample size of population 1.
#' @param n2 integer, sample size of population 1.
#' @param m0 numeric, null value for the difference in means to be tested for. Default is 0.
#' @param equal.variance logical, whether or not to assume equal variance. Default is FALSE.
#'
#' @return vector of `Difference of means`, `Std Error`, `t`, `p-value` of the t-test from summary data.
#' @export
#' @importFrom stats pt
#'
#' @examples
#' t_test_from_summary(15,18,2,1.8,20,23)
t_test_from_summary<- function(m1,m2,sd1,sd2,n1,n2,m0=0,equal.variance=FALSE){
  if(!equal.variance){
    se <- sqrt( (sd1^2/n1) + (sd2^2/n2) )
    # welch-satterthwaite df
    df <- ( (sd1^2/n1 + sd2^2/n2)^2 )/( (sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1) )
  }else{
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*sd1^2 + (n2-1)*sd2^2)/(n1+n2-2) )
    df <- n1+n2-2
  }
  t <- (m1-m2-m0)/se
  dat <- c(m1-m2, se, t, 2*stats::pt(-abs(t),df))
  names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
  dat
}

#' Modify a vector to be probabilities in ]0,1[
#'
#' @param P a vector of values not all in ]0,1[
#'
#' @return a vector of valid probabilities in ]0,1[
#' @export
#'
#' @examples
#' # internal use in make_obs.prob().
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

