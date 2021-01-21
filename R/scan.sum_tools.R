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

# sum_sampled<- function(scan,method){
#   switch(method,
#          "theoretical" = {
#            theoretical.sampled <- scan$theoretical.sum
#            if (length(scan$scans.to.do) == 1) {if (scan$scans.to.do == "all") {scan$scans.to.do <- 1:scan$total_scan}}
#            theoretical.sampled[scan$Adj.subfun(theoretical.sampled)] <- length(scan$scans.to.do)
#            theoretical.sampled
#          },
#          ""
#   )
# }
