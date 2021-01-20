# Single scan related functions -------------------------------------------

#' Generator for `scan` objects
#' Internal use. The user should rather rely on `simu_scan` as a wrapper for the
#' different steps needed to perform the scan from inputted data.
#'
#' @param presence.prob a `presenceProb` object
#' @param scans.to.do Optional. Only required if inputted `focal` is a `focalList` object. Either:
#'  \itemize{
#'   \item{an integer vector included in `1:total_scan` of the scans to perform}
#'   \item{the special case `"all"` (default) sets `scans.to.do` to `1:total_scan` and set the simulation to perform all the scans}
#' }
#'
#' @return a `scan` object (S3 class) containing:
#' \itemize{
#'   \item{`raw.scan`: a list of raw binary adjacency matrix shaped like the
#'   `Adj` contained in `presence.prob`, considered directed in the
#'   algorithm}
#'   \item{`theoretical.scan`: a list of binary adjacency matrix, where all ties were
#'   observed but the `mode` has been applied}
#'   \item{`scan.type`: character scalar. `generate_scan` sets it to
#'   "theoretical", `sample_from_scan` will set it to "empirical" and append the
#'   empirical matrix}
#'   \item{`Adj`: `Adj` data contained in `presence.prob`}
#'   \item{`total_scan`: `total_scan` data contained in `presence.prob`}
#'   \item{`scans.to.do`: inputted `scans.to.do`}
#'   \item{`mode`: `mode` data contained in `presence.prob`}
#'   \item{`weighted`: logical, at this stage can only be `TRUE` if `mode =
#'   plus` (some edges can become `2`)}
#'   \item{`Adj.subfun`: `Adj.subfun` data contained in `presence.prob`}
#'   \item{`presence.prob`: `presence.prob$P` (only the probability matrix) data
#'   contained in `presence.prob`}
#' }
#'
#' @noRd
generate_scan<- function(presence.prob,scans.to.do){
  if(length(scans.to.do) == 1) {if(scans.to.do == "all") {scans.to.do <- 1:presence.prob$total_scan}}
  raw.scan<- draw_raw.scan(presence.prob = presence.prob,scans.to.do = scans.to.do)
  scan<- list(
    raw.scan = raw.scan,
    theoretical.scan = apply_mode(raw.scan = raw.scan,mode = presence.prob$mode),
    scan.type = "theoretical",
    Adj = presence.prob$Adj,
    total_scan = presence.prob$total_scan,
    scans.to.do = scans.to.do,
    mode = presence.prob$mode,
    weighted = ifelse(presence.prob$mode == "plus",FALSE,TRUE), # only case for which an edge can be > 1
    Adj.subfun = presence.prob$Adj.subfun,
    presence.prob = presence.prob$P # in here it is not a presenceProb object anymore, to avoid storing redundant variables
  )
  class(scan)<- "scan"
  scan
}

#' Print method for `scan` objects
#' @export
#' @noRd
print.scan<- function(x,...){
  cat("Scan type: ",x$scan.type,", mode: ",x$mode,"\n\n",sep = "")
  n <- length(x$theoretical.scan)
  if(length(x$theoretical.scan) >= 10) {
    print.default(x$theoretical.scan[1:3],...)
    cat("... (",n-5," more scans)\n\n\n",sep = "")
    cat("[[",n-1,"]]\n",sep = "")
    print.default(x$theoretical.scan[[(n-1)]],...)
    cat("[[",n,"]]\n",sep = "")
    print.default(x$theoretical.scan[[n]],...)
  } else {
    print.default(x$theoretical.scan,...)
  }
}

#' Plot method for `scan` objects
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph plot.igraph
#' @export
#' @noRd
plot.scan<- function(x,...){
  x<- igraph::graph_from_adjacency_matrix(x$theoretical,mode = x$mode,weighted = x$weighted)
  igraph::plot.igraph(x,...,main = "theoretical")
}

#' Draw raw binary scan(s) from presence.prob
#' Internal use in `generate_scan`. Draw the raw scan as an always directed binary adjacency matrix. This way, keeps track of how the theoretical scan has been obtained, before applying `apply_mode` to this raw scan
#'
#' @param presence.prob a `presenceProb` object
#' @param scans.to.do either:
#'  \itemize{
#'   \item{an integer vector included in `1:total_scan` of the scans to perform}
#'   \item{the special case `"all"` (default) sets `scans.to.do` to `1:total_scan` and set the simulation to perform all the scans}
#' }
#'
#' @return a list of raw binary adjacency matrix shaped like the `Adj` contained in `presence.prob`
#'
#' @noRd
draw_raw.scan <- function(presence.prob,scans.to.do){
  n<- nrow(presence.prob$Adj);nodes_names<- rownames(presence.prob$Adj)
  presence.P.vec<- presence.prob$P[presence.prob$Adj.subfun(presence.prob$P)]; p<- length(presence.P.vec)  # subset a presence probability vector the same way Adj is subset (depends on `Adj`'s mode (cf. igraph)) <- NOW THIS IS JUST A UPPER/LOWER.TRI VS NON.DIAGONAL THING...
  raw.scan <-
    rep(
      list( # required for rep to output a list
        matrix(0,nrow = n,ncol = n,dimnames = list(nodes_names,nodes_names))  # structure the scan as a matrix filled with zeros
      ),
      length(scans.to.do)
    )
  lapply(
    raw.scan,
    function(s) {
      s[presence.prob$Adj.subfun(s)]<- stats::rbinom(p,1,presence.P.vec)  # core of the randomization: draw a (raw.scan) tie or not for each (relevant, cf. triangular matrices or undirected) dyad according to its presence probability
      s
    }
  )
}
