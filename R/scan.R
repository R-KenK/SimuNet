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
#'   \item{`raw.scan.list`: a list of raw binary adjacency matrix shaped like the
#'   `Adj` contained in `presence.prob`, considered directed in the
#'   algorithm}
#'   \item{`theoretical.scan.list`: a list of binary adjacency matrix, where all ties were
#'   observed but the `mode` has been applied}
#'   \item{`scan.type`: character scalar. `generate_scan` sets it to
#'   "theoretical", `sample_from_scan` will set it to "empirical" and append the
#'   empirical matrix}
#'   \item{`Adj`: `Adj` data contained in `presence.prob`}
#'   \item{`total_scan`: `total_scan` data contained in `presence.prob`}
#'   \item{`scans.to.do`: inputted `scans.to.do`}
#'   \item{`mode`: `mode` data contained in `presence.prob`}
#'   \item{`Adj.subfun`: `Adj.subfun` data contained in `presence.prob`}
#'   \item{`presence.prob`: `presence.prob$P` (only the probability matrix) data
#'   contained in `presence.prob`}
#' }
#'
#' @noRd
generate_scan<- function(presence.prob,scans.to.do){
  if(length(scans.to.do) == 1) {if(scans.to.do == "all") {scans.to.do <- 1:presence.prob$total_scan}}
  raw.scan.list <- draw_raw.scan.list(presence.prob = presence.prob,scans.to.do = scans.to.do)
  theoretical.scan.list <- apply_mode(raw.scan.list = raw.scan.list,mode = presence.prob$mode)
  scan<- list(
    raw.scan.list = raw.scan.list,
    theoretical.scan.list = theoretical.scan.list,
    scan.type = "theoretical",
    Adj = presence.prob$Adj,
    total_scan = presence.prob$total_scan,
    scans.to.do = scans.to.do,
    mode = presence.prob$mode,
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
  scans.to.do <- explicit_scans.to.do(x)
  n <- length(scans.to.do)
  # print the general simulation info
  if (n >= 15) {
    scans.to.do <- paste(do.call(paste,as.list(scans.to.do)[1:5]),
                         "...",
                         do.call(paste,as.list(scans.to.do)[(n-4):n])
    )
  } else {
    scans.to.do <- scans.to.do
  }
  cat("\nScan(s) performed: ",scans.to.do,sep=" ")
  cat("\nScan type: ",x$scan.type,", mode: ",x$mode,"\n\n",sep = "")
  shorten_list.to.print(x$theoretical.scan.list)
}

#' Summary method for `scan` objects
#' @export
#' @noRd
summary.scan <- function(object,...) {
  scans.to.do <- explicit_scans.to.do(object)
  mode <- object$mode
  theoretical.sum <- sum_scan.list(object$theoretical.scan.list)
  theoretical.sampled <- sum_scan.sampled(object,method = "theoretical")
  theoretical.scaled <- theoretical.sum/ifelse(theoretical.sampled != 0,theoretical.sampled,1)

  scan.summary <- list(
    theoretical.sum = theoretical.sum,
    theoretical.sampled = theoretical.sampled,
    theoretical.scaled = theoretical.scaled,
    scans.to.do = scans.to.do,
    mode = mode
    # more things can be added here# more things can be added here
  )
  class(scan.summary) <- "summary.scan"
  scan.summary
  # more info to be displayed
}

#' Print method for `summary.scan` objects
#' @export
#' @noRd
print.summary.scan<- function(x,scaled = FALSE,...){
  cat("Theoretical weighted adjacency matrix:\n")
  if (scaled) {
    to.print <- x$theoretical.scaled
  } else {
    to.print <- x$theoretical.sum
  }
  print.default(to.print,...)
  cat(paste0("\nobtained after summing ", length(x$scans.to.do), " binary scans (mode = \"", x$mode,"\")", "\n\n"))
}

#' Plot method for `scan` objects
#' @export
#' @noRd
plot.scan<- function(x,...){
  x.summary <- summary.scan(x)
  plot.summary.scan(x.summary,...)
}

#' Plot method for `summary.scan` objects
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph plot.igraph
#' @export
#' @noRd
plot.summary.scan<- function(x,scaled = FALSE,vertex.size = NULL,vertex.size.mul = nrow(x$theoretical.sum)/2,vertex.size.min = 3,vertex.size.fun = compute.strength,...){
  if (scaled) {theoretical.adj <- x$theoretical.scaled} else {theoretical.adj <- x$theoretical.sum}
  total_scan <- max(x$theoretical.sampled,na.rm = TRUE)
  max.adj <- max(theoretical.adj,na.rm = TRUE)
  x.G <- igraph::graph_from_adjacency_matrix(theoretical.adj,mode = x$mode,weighted = TRUE)
  if (is.null(vertex.size)) {
    vertex.size <- vertex.size.fun(x.G,mode = x$mode)
    vertex.size <- (vertex.size.mul * (vertex.size / max(vertex.size))) + vertex.size.min
  }
  igraph::plot.igraph(x.G,...,main = "Scan type: theoretical",
                      sub = paste0("obtained after summing ", total_scan, " binary scans (mode = \"", x$mode,"\")"),
                      vertex.size = vertex.size
                      )
}

#' Compute node degree from graph or adjacency matrix
#'
#' @param graph an igraph object (or an adjacency matrix)
#' @param mode optional, only if `graph` is an adjacency matrix. Othewise character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
#'
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph degree
#' @importFrom igraph vertex_attr

#' @return a vector of degree values for each node
#' @noRd
compute.deg<- function(graph,mode=NULL){
  if(is.matrix(graph)){graph<- igraph::graph_from_adjacency_matrix(graph,mode = mode,weighted = TRUE,add.colnames = TRUE)}
  deg<- igraph::degree(graph)
  if(!is.null(names(deg))) {names(deg)<- igraph::vertex_attr(graph)[[1]]} # dirty: does not actually test if the order of the vertex centrality is the same as the name, but I suspect igraph does that by default...
  deg
}

#' Compute node strength from graph or adjacency matrix
#'
#' @param graph an igraph object (or an adjacency matrix)
#' @param mode optional, only if `graph` is an adjacency matrix. Otherwise character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
#'
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph strength
#' @importFrom igraph vertex_attr
#'
#' @return a vector of strength values for each node
#' @noRd
compute.strength<- function(graph,mode = NULL){
  if (is.matrix(graph)) {graph <- igraph::graph_from_adjacency_matrix(graph,mode = mode,weighted = TRUE,add.colnames = TRUE)}
  stren <- igraph::strength(graph)
  if (!is.null(names(stren))) {names(stren) <- igraph::vertex_attr(graph)[[1]]} # dirty: does not actually test if the order of the vertex centrality is the same as the name, but I suspect igraph does that by default...
  stren
}

#' Compute eigenvector centrality values from graph or adjacency matrix
#'
#' @param graph an igraph object (or an adjacency matrix)
#' @param mode optional, only if `graph` is an adjacency matrix. Othewise character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.

#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph strength
#' @importFrom igraph vertex_attr
#' @importFrom igraph E
#'
#' @return a vector of eigenvector centrality values for each node
#' @noRd
compute.EV<- function(graph,mode = NULL){
  if (is.matrix(graph)) {graph <- igraph::graph_from_adjacency_matrix(graph,mode = mode,weighted = TRUE,add.colnames = TRUE)}
  EV <- igraph::eigen_centrality(graph, weights = igraph::E(graph)$weight,scale = FALSE)$vector
  if (!is.null(names(EV))) {names(EV) <- igraph::vertex_attr(graph)[[1]]} # dirty: does not actually test if the order of the vertex centrality is the same as the name, but I suspect igraph does that by default...
  EV
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
draw_raw.scan.list <- function(presence.prob,scans.to.do){
  n<- nrow(presence.prob$Adj);nodes_names<- rownames(presence.prob$Adj)
  presence.P.vec<- presence.prob$P[presence.prob$Adj.subfun(presence.prob$P)]; p<- length(presence.P.vec)  # subset a presence probability vector the same way Adj is subset (depends on `Adj`'s mode (cf. igraph)) <- NOW THIS IS JUST A UPPER/LOWER.TRI VS NON.DIAGONAL THING...
  raw.scan.list <-
    rep(
      list( # required for rep to output a list
        matrix(0,nrow = n,ncol = n,dimnames = list(nodes_names,nodes_names))  # structure the scan as a matrix filled with zeros
      ),
      length(scans.to.do)
    )
  lapply(
    raw.scan.list,
    function(s) {
      s[presence.prob$Adj.subfun(s)]<- stats::rbinom(p,1,presence.P.vec)  # core of the randomization: draw a (raw.scan.list) tie or not for each (relevant, cf. triangular matrices or undirected) dyad according to its presence probability
      s
    }
  )
}
