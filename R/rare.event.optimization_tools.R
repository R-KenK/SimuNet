# Rare events optimization related functions ------------------------------

#' Generate new data using orthogonal polynoms of input variable
#'
#' @param n integer, number of node of the network
#' @param total_scan integer, sampling effort
#' @param max.obs integer, maximum number of observation after `total_scan` scans
#' @param algorithm character, which algorithm to produce new data for.
#'
#' @return a data.frame inputed on the same orthogonal polynomial scale as the one used by the standard and optimized models of expected time
#' @export
#' @importFrom stats predict
#'
#' @examples
#' # Internal use in opti.expected.time
new.data.poly<- function(n,total_scan,max.obs,algorithm = c("standard","optimization for rare event")){
  algorithm<- match.arg(algorithm)
  switch(algorithm,
         "standard" = {
           n.poly<- n.poly.std
           total.poly<- total.poly.std
           max.poly<- max.poly.std
         },
         "optimization for rare event" = {
           n.poly<- n.poly.opt
           total.poly<- total.poly.opt
           max.poly<- max.poly.opt
         }
  )
  cbind(data.frame(n = n,total_scan = total_scan,max.obs = max.obs),
        n = stats::predict(n.poly,n),
        total_scan = stats::predict(total.poly,total_scan),
        max.obs = stats::predict(max.poly,max.obs))
}

#' Calculating expected time using the standard scan method from glm model
#'
#' @param n number of node of the original network
#' @param total_scan sampling effort
#' @param max.obs maximum value of the weighted adjacency matrix of the original network
#' @param se.fit logical, should the standard error be computed by the predict() function?
#'
#' @return the result of the expected time model predictions (in milliseconds), i.e. a vector with values `fit`, `se.fit`, and `residual.scale`.
#' @export
#' @importFrom stats predict
#'
#' @examples
#' # Internal use in decide_use.rare.opti
standard.expected.time<- function(n,total_scan,max.obs,se.fit = TRUE){
  new.dat<- new.data.poly(n = n,total_scan = total_scan,max.obs = max.obs,algorithm = "standard")
  stats::predict(standard.model,newdata = new.dat,type = "response",se.fit = se.fit)
}

#' Calculating expected time using the optimization for rare event from glm model
#'
#' @param n number of node of the original network
#' @param total_scan sampling effort
#' @param max.obs maximum value of the weighted adjacency matrix of the original network
#' @param se.fit logical, should the standard error be computed by the predict() function?
#'
#' @return the result of the expected time model predictions (in milliseconds), i.e. a vector with values `fit`, `se.fit`, and `residual.scale`.
#' @export
#' @importFrom stats predict
#'
#' @examples
#' # Internal use in decide_use.rare.opti
opti.expected.time<- function(n,total_scan,max.obs,se.fit = TRUE){
  new.dat<- new.data.poly(n = n,total_scan = total_scan,max.obs = max.obs,algorithm = "opti")
  stats::predict(opti.model,newdata = new.dat,type = "response",se.fit = se.fit)
}

#' Decide based on expected times if the otpimization for rare event should be used
#'
#' @param n number of node of the original network, or can be an adjacency matrix
#' @param total_scan integer, sampling effort.
#' @param max.obs maximum value of the weighted adjacency matrix of the original network. Optional if n is inputted as an adjacency matrix.
#' @param alpha numerical in [0,1], type I error significance level.
#'
#' @return a logical value meaning that the optimization for rare events should be use when TRUE is returned.
#' @export
#'
#' @examples
#' decide_use.rare.opti(30,900,10)
#' decide_use.rare.opti(30,9000,10)
decide_use.rare.opti<- function(n,total_scan,max.obs=NULL,alpha=0.05){
  if(is.null(max.obs)&is.matrix(n)){max.obs<- max(n)}
  if(is.matrix(n)){n<- nrow(n)}
  expected.time.std<- standard.expected.time(n = n,total_scan = total_scan,max.obs = max.obs)
  expected.time.opt<- opti.expected.time(n = n,total_scan = total_scan,max.obs = max.obs)

  if(any(is.infinite(expected.time.std$fit),is.infinite(expected.time.opt$fit))){return(FALSE)}

  test<- t_test_from_summary(m1 = expected.time.std$fit,
                             m2 = expected.time.opt$fit,
                             sd1 = expected.time.std$se.fit,
                             sd2 = expected.time.opt$se.fit,
                             n1 = n,n2 = n,
                             m0 = 0,equal.variance = FALSE)

  if(test[["Difference of means"]]>0){ # m1-m2 <=> if optimization is faster
    test[["p-value"]]<=alpha # If difference between optimization and standard is significative. Otherwise favor the standard method
  }else{
    FALSE
  }
}


#' Simulate which scan returns an all-zeros matrix
#'
#' @param total_scan integer, sampling effort
#' @param presence.prob presence probability matrix (or vector)
#'
#' @return a list of NULL representing the non-zero scan to run with an attribute `n.zero` being the number of full-zero scans,and `non.zero.pos`
#' @export
#' @importFrom stats rbinom
#'
#' @examples
#' # Internal use
simulate_zeros.non.zeros<- function(total_scan,presence.prob){
  zero.non.zero.list<- stats::rbinom(total_scan,1,1-prod(1-presence.prob))==1
  zero.non.zero.table<- table(ifelse(zero.non.zero.list,"n.non.zeros","n.zeros"));
  scan_list<- vector(mode="list",length = zero.non.zero.table["n.non.zeros"])
  attr(scan_list,"n.zeros")<- zero.non.zero.table["n.zeros"]
  attr(scan_list,"non.zero.list")<- which(zero.non.zero.list)
  scan_list
}

#' Determine step-by-step conditional probabilities for non-zeros scans
#' Internal use. Returns the CDF of the probability that the i-th dyad (with probability presence.prob[i]) is the first to yield a 1.
#'
#' @param presence.prob presence probability matrix (or vector)
#'
#' @return a cumulative distribution function of the probability that the i-th dyad (with probability presence.prob[i]) is the first to yield a 1
#' @export
#'
#' @details Workflow is as follows: first simulate.zeros.non.zeros() determines which scans are all-zeros and which are non-zeros. Then for non zeros, at a random order each dyad is drawn in order with conditional probability that: (1) there is at least one 1 in the scan, and (2) all the previous coins were zeros. Once the first one is drawn, the rest are drawn with their regular probabilities. In details, cumulative density probability of each dyad (in a given random order) to be the first one to be a 1 is calculated, and a random draw determine which one is first, set the previous ones to zero, and draw the rest normally. cf. Supplmenentary material X.
#'
#' @examples
#' # Internal use.
adjust.conditional.prob<- function(presence.prob){
  prob.all.zeros<- 1-prod(1-presence.prob)
  previous.are.zeros<- c(1,cumprod(1-presence.prob[1:(length(presence.prob)-1)]))
  sapply(1:length(presence.prob),function(i) presence.prob[i]/prob.all.zeros*previous.are.zeros[i])
}
