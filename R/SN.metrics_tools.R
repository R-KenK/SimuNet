#' Calculate a social network metric from simulation
#'
#' @param scan a scan or empiScan object
#' @param method character vector, among "theoretical","group", and "focal", on which to calculate the network metric
#' @param on.what character scalar, among "scan.list","sum","scaled", describing the type of matrix/network on which the function should be applied
#' @param SNm.fun function of (matrix,mode) (for the moment) that returns a named numeric vector of the metric(s) of interest
#' @param output character scalar, among "vector" and "list" (for the moment), to know if the vector of metric should be concatenated into a vector, or outputted as a list of vector
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
#' calculate_SNm(group.constant.scans,method = c("theoretical","group"),"scaled",compute.strength,output = "list")
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

