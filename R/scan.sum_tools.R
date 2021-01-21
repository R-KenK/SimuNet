# Sum up scan routine transposition to new framework ----------------------

#' Sum up `scan.list`s into weighted adjacency matrices
#'
#' @param scan.list
#'
#' @return
#' @noRd
sum_scan.list<- function(scan.list){
  Reduce(matrix_sum_na.rm,scan.list) # sum the scan list, considering NAs as zeros
}

#' TO WRITE
#'
#' @param scan TO WRITE
#' @param method TO WRITE
#'
#' @return TO WRITE
#' @noRd
sum_scan.sampled<- function(scan,method = c("theoretical","group","focal")) {
  method <- match.arg(method)
  switch(method,
         "theoretical" = {
           theoretical.sampled <- scan$Adj
           non.diagonal(theoretical.sampled) <- length(explicit_scan.to.do(scan))
           theoretical.sampled
         },
         "group" = count_NA(scan$group.scan.list),
         "focal" = count_NA(scan$focal.scan.list)
  )
}

#' TO WRITE
#'
#' @param scan TO WRITE
#'
#' @return TO WRITE
#' @noRd
explicit_scan.to.do <- function(scan) {
  if (length(scan$scans.to.do) == 1) {
    if (scan$scans.to.do == "all") {
      return(1:scan$total_scan)
    }
  }
  scan$scans.to.do
}
