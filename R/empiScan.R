# Empirical sampling related functions ------------------------------------

#' Generator for `empiScan` objects
#'
#' @param scan a `scan` object with `scan.type = "theoretical"`
#' @param sampling.param a `samplingParam` object containing:
#' \itemize{
#'   \item{method}{inputted `method`}
#'   \item{mode}{inputted `mode`}
#'   \item{obs.prob}{inputted `obs.prob`}
#'   \item{focal}{inputted `focal`}
#' }
#'
#' @return an `empiScan` object (S3 class), inheriting from `scan` object, containing:
#' \itemize{
#'   \item{raw}{a binary adjacency matrix, considered directed in the algorithm}
#'   \item{theoretical}{a binary adjacency matrix, where all ties were observed but the `mode` has been applied}
#'   \item{scan.type}{set to `"empirical"` by `sample_from_scan`}
#'   \item{method}{from inputted `sampling.param`}
#'   \item{group}{an adjacency matrix, where the observation probability `obs.prob` from `sampling.param` of each dyad has been applied. `NULL` if `method = "focal"`}
#'   \item{focal}{an adjacency matrix, where only the selected `focal` from `sampling.param` is visible. `NULL` if `method = "focal"`}
#'   \item{Adj}{`Adj` data contained in `presence.prob`}
#'   \item{total_scan}{`total_scan` data contained in `presence.prob`}
#'   \item{mode}{`mode` data contained in `presence.prob`}
#'   \item{weighted}{logical, at this stage can only be `TRUE` if `mode = plus` (some edges can become `2`)}
#'   \item{Adj.subfun}{`Adj.subfun` data contained in `presence.prob`}
#'   \item{presence.prob}{`presence.prob$P` (only the probability matrix) data contained in `presence.prob`}
#' }
#' @noRd
generate_empiScan<- function(scan,sampling.param){
  scan$scan.type<- "empirical"
  scan$method<- sampling.param$method
  scan$sampling.param<- sampling.param$sampling.param
  if(!is.null(sampling.param$obs.prob)) {
    scan$group<- sample_from_scan(scan = scan,sampling.param = sampling.param,method = "group")
  }
  if(!is.null(sampling.param$focal)) {
    scan$focal<- sample_from_scan(scan = scan,sampling.param = sampling.param,method = "focal")
  }
  class(scan)<- c("empiScan","scan")
  scan
}

#' Print method for `scan` objects
#' @export
#' @noRd
print.empiScan<- function(x,...){
  cat("Scan type: theoretical\n\n")
  print.default(x$theoretical,...)
  if(!is.null(x$group)) {
    cat("\n\nScan type: group scan\n\n")
    print.default(x$group,...)
  }
  if(!is.null(x$focal)) {
    cat("\n\nScan type: focal scan\n\n")
    print.default(x$focal,...)
  }
}

#' Test if object if a `empiScan` object
#'
#' @param scan an object to test.
#'
#' @return logical, TRUE if the inputted object is a `empiScan` object.
#'
#' @noRd
is.empiScan<- function(scan){
  inherits(scan,"empiScan")
}

#' Plot method for `scan` objects
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph plot.igraph
#' @export
#' @noRd
plot.empiScan<- function(x,...){ # Need a way to make it print two in case of method = "both"
  x<- igraph::graph_from_adjacency_matrix(x$theoretical,mode = x$mode,weighted = x$weighted)
  igraph::plot.igraph(x,...)
}

#' Empirically sample from a theoretical scan
#'
#'
#' @param scan a `scan` object with `scan.type = "theoretical"`
#' @param sampling.param a `samplingParam` object containing:
#' \itemize{
#'   \item{method}{inputted `method`}
#'   \item{mode}{inputted `mode`}
#'   \item{obs.prob}{inputted `obs.prob`}
#'   \item{focal}{inputted `focal`}
#' }
#'
#' @return
#' @noRd
sample_from_scan<- function(scan,sampling.param,method){
  # according to the empirical method chosen, applies some "empirical" missed observation via the correct sampling method
  switch(method,
         "group" = group_sample(scan = scan,obs.prob = sampling.param$obs.prob),
         "focal" = focal_sample(scan = scan,focal = sampling.param$focal)
  )
}

#' Perform a group scan sampling
#' Internal use. Presently, the sampling is performed on the _theoretical_ scan
#'
#' @param scan a `scan` object with `scan.type = "theoretical"`
#' @param obs.prob an `obsProb` object
#'
#' @importFrom stats rbinom
#'
#' @return a binary adjacency matrix with potentially some dyads not observed (turned into `NA`)
#' @noRd
group_sample<- function(scan,obs.prob){
  observed<- scan$theoretical  # set the sampled scan to be like the `raw` one (i.e. this is a binary _directed_ adjacency matrix)
  obs.P<- obs.prob$P[scan$Adj.subfun(obs.prob$P)]  # subset obs.prob like the adjacency matrix (i.e. triangular matrix) into a vector of observation probabilities
  missed<- stats::rbinom(length(obs.P),1,obs.P)==0  # draw the missed observation
  observed[scan$Adj.subfun(observed)][missed]<- NA  # set the missed observation to `NA`
  observed
}

#' Perform a group scan sampling
#' Internal use. Presently, the sampling is performed on the _raw_ scan
#'
#' @param scan a `scan` object with `scan.type = "theoretical"`
#' @param obs.prob an `obsProb` object
#'
#' @importFrom stats rbinom
#'
#' @return a binary adjacency matrix with potentially some dyads not observed (turned into `NA`)
#' @noRd
group_sample.old<- function(scan,obs.prob){
  observed<- scan$raw  # set the sampled scan to be like the `raw` one (i.e. this is a binary _directed_ adjacency matrix)
  obs.P<- obs.prob$P[scan$Adj.subfun(obs.prob$P)]  # subset obs.prob like the adjacency matrix (i.e. triangular matrix) into a vector of observation probabilities
  missed<- stats::rbinom(length(obs.P),1,obs.P)==0  # draw the missed observation
  observed[scan$Adj.subfun(observed)][missed]<- NA  # set the missed observation to `NA`
  apply_mode(observed,mode = scan$mode)
}

#' Perform a focal scan sampling
#' Internal use.
#'
#' @param scan a `scan` object with `scan.type = "theoretical"`
#' @param focal an `focal` object
#'
#' @return a binary adjacency matrix with dyads not in the `focal$focal` row and column turned into `NA`
#' @noRd
focal_sample<- function(scan,focal){
  observed<- scan$theoretical  # set the sampled scan to be like the `raw` one (i.e. this is a binary _directed_ adjacency matrix)
  observed[-focal$focal,-focal$focal]<- NA # set other rows and columns than those of the `focal` to `NA`
  observed
}

#' Perform a focal scan sampling
#' Internal use.
#'
#' @param scan a `scan` object with `scan.type = "theoretical"`
#' @param focal an `focal` object
#'
#' @return a binary adjacency matrix with dyads not in the `focal$focal` row and column turned into `NA`
#' @noRd
focal_sample.old<- function(scan,focal){
  observed<- scan$raw  # set the sampled scan to be like the `raw` one (i.e. this is a binary _directed_ adjacency matrix)
  observed[-focal$focal,-focal$focal]<- NA # set other rows and columns than those of the `focal` to `NA`
  apply_mode(observed,mode = scan$mode)
}

