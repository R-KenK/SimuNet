# Sum up scan routine transposition to new framework ----------------------

#' Sum up `scan.list` into weighted adjacency matrix
#'
#' @param scan.list a list of binary adjacency matrices, where an unobserved dyad (whether it is theoretically 0 or 1) is `NA`
#'
#' @return a weighted adjacency matrix representing the total number of occurrences where each dyad has been observed associated (and where `NA`s were considered as zero)
#' @noRd
sum_scan.list<- function(scan.list){
  Reduce(matrix_sum_na.rm,scan.list) # sum the scan list, considering NAs as zeros
}

#' Sum up `scan.list` into weighted sampling effort matrix
#'
#' @param scan a `scan` object (thus containing among others `Adj`, `scans.to.do`, and potentially either (or both) `group.scan.list` and `focal.scan.list`)
#' @param method character scalar, either `"theoretical"`, `"group"`, or `"focal"`
#'
#' @return an integer matrix representing the total number of occurrences where each dyad has been observed - whether it being `0` or `1` (and where `NA`s were considered as zero)
#' @noRd
sum_scan.sampled<- function(scan,method = c("theoretical","group","focal")) {
  method <- match.arg(method)
  X.sampled <- switch(method,
         "theoretical" = {
           theoretical.sampled <- scan$Adj
           non.diagonal(theoretical.sampled) <- length(explicit_scans.to.do(scan))
           theoretical.sampled
         },
         "group" = count_nonNA(scan$group.scan.list,scan$Adj.subfun),
         "focal" = count_nonNA(scan$focal.scan.list,scan$Adj.subfun)
  )
  X.sampled[!scan$Adj.subfun(X.sampled)] <- 0L
  X.sampled
}

#' Explicit the vector of scans to do
#' Internal use. Explicit the case where `scans.to.do = "all"` to be `1:total_scan`
#'
#' @param scan a `scan` object (or an object with a `scans.to.do` S3 class attribute)
#'
#' @return the integer vector `scans.to.do` contained in `scan`, or its explicit form `1:total_scan` in the case where `scans.to.do = "all"`
#' @noRd
explicit_scans.to.do <- function(scan) {
  if (length(scan$scans.to.do) == 1) {
    if (scan$scans.to.do == "all") {
      return(1:scan$total_scan)
    }
  }
  scan$scans.to.do
}
