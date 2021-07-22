#' Calculate a social network metric from simulation
#'
#' @param scan a scan or empiScan object
#' @param method character vector, among "theoretical","group", and "focal", on which to calculate
#'   the network metric
#' @param on.what character scalar, among "scan.list","sum","scaled", describing the type of
#'   matrix/network on which the function should be applied
#' @param SNm.fun function of (matrix,mode) (for the moment) that returns a named numeric vector of
#'   the metric(s) of interest
#' @param output character scalar, among "vector" and "list" (for the moment), to know if the vector
#'   of metric should be concatenated into a vector, or outputted as a list of vector
#'
#' @return a named vector or list of vectors of the metric of interest
#' @export
#'
#' @examples
#' set.seed(42)
#'
#' n<- 5;nodes<- letters[1:n];total_scan<- 42;
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- sample(0:total_scan,n*(n-1),replace = TRUE)
#' Adj
#'
#' theo.scans <- simu_scan(Adj,total_scan,mode = "min",scans.to.do = 1:3)
#' calculate_SNm(theo.scans,method = "theoretical","sum",compute.strength)
#'
#' para.group.constant<- simu_samplingParam(Adj,total_scan,mode =
#'                                          "min",group.scan_param = 0.42)
#' group.constant.scans <- simu_scan(sampling.param = para.group.constant)
#' calculate_SNm(group.constant.scans,
#'   method = c("theoretical","group"),"scaled",
#'   compute.strength,output = "list"
#'  )
calculate_SNm <- function(scan,method = c("theoretical","group","focal"),
                          on.what = c("scan.list","sum","scaled"),SNm.fun,output = c("vector","list")) {
  # method <- match.arg(method)
  on.what <- match.arg(on.what)
  output <- match.arg(output)
  mode <- scan$mode

  arg.name <- paste0(method,".",on.what)
  scan.arg <- switch(on.what,
                     "sum" = ,
                     "scaled" = {
                       summary(scan)
                     },
                     "scan.list" = {
                       scan
                     }
  )
  l <- compute.SNm.list(scan.arg,arg.name,mode,SNm.fun)
  format_output(output,l)
}
# early recycling from validation script
# extract_SNm <- function(Adj.list,SNm_fun,transform.to.igraph = FALSE,mode = NULL,weighted = TRUE,...,Xapply = lapply) {
#   if (!is.list(Adj.list)) {Adj.list <- list(Adj.list)}
#   if (transform.to.igraph) {
#     Adj.list <- lapply(Adj.list,igraph::graph.adjacency,mode = mode,weighted = weighted)
#   }
#   Xapply(Adj.list,SNm_fun,...)
# }
#
#
# # wrapper to extract several metrics passed as a vector of functions ------
#
#
# extract_SNm.vec <- Vectorize(FUN = extract_SNm,
#                              vectorize.args = c("SNm_fun","transform.to.igraph"),
#                              SIMPLIFY = FALSE)

#' TO WRITE
#'
#' @param scan.arg TO WRITE
#' @param arg.name TO WRITE
#' @param mode TO WRITE
#' @param SNm.fun TO WRITE
#'
#' @return TO WRITE
#' @noRd
compute.SNm.list <- function(scan.arg,arg.name,mode,SNm.fun) {
  if (inherits(scan.arg,"scan")) {
    lapply(
      arg.name,
      function(a) {
        lapply(
          scan.arg[[a]],
          function(s) {
            SNm.fun(s,mode)
          }
        )
      }
    )
  } else {
    lapply(
      arg.name,
      function(a) {
        SNm.fun(scan.arg[[a]],mode)
      }
    )
  }
}

#' TO WRITE
#'
#' @param output TO WRITE
#' @param scan.arg TO WRITE
#' @param arg.name TO WRITE
#' @param mode TO WRITE
#'
#' @return TO WRITE
#' @noRd
format_output <- function(output,l) {
  switch(output,
         "vector" = do.call(c,l),
         "list" = l
  )
}


# additional SN metrics ---------------------------------------------------

#' TO WRITE
#'
#' @param Adj TO WRITE
#' @param mode TO WRITE
#'
#' @return TO WRITE
#' @importFrom DirectedClustering ClustF
#'
#' @noRd
GlobalCC <- function(Adj,mode) {
  if (inherits(Adj,"igraph")) {
    Adj <- igraph::get.adjacency(Adj,type = "both",sparse = FALSE,attr = "weight")
  }
  switch(mode,
         "upper" = ,
         "lower" = {
           Adj <- Adj + t(Adj)
           type <- "undirected"
         },
         "directed" = type <- "directed",
         type <- "undirected"
  )
  DirectedClustering::ClustF(Adj,type = type)[["GlobalCC"]]
}

#' TO WRITE
#'
#' @param G TO WRITE
#' @param mode TO WRITE
#' @param centrality.fun TO WRITE
#'
#' @return TO WRITE
#' @noRd
sort_central <- function(G,mode,centrality.fun) {
  sort(centrality.fun(G,mode),decreasing = TRUE)
}

#' TO WRITE
#'
#' @param G TO WRITE
#' @param mode TO WRITE
#' @param centrality.fun TO WRITE
#'
#' @return TO WRITE
#' @noRd
get_central.node <- function(G,mode,centrality.fun) {
  names(sort_central(G,mode,centrality.fun)[1])
}

#' TO WRITE
#'
#' @param G TO WRITE
#' @param mode TO WRITE
#' @param centrality.fun TO WRITE
#' @param central.node TO WRITE
#'
#' @return TO WRITE
#' @noRd
remove_central.node <- function(G,mode = NULL,centrality.fun = NULL,central.node = NULL) {
  central.node <- get_central.node(G = G,mode = mode,centrality.fun = centrality.fun)
  to.keep <- row.names(G) != central.node
  G[to.keep,to.keep]
}
