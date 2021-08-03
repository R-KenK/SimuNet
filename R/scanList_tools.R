#' Generator for `scanList` objects
#' Internal use. The user should rather rely on `simunet` as a wrapper for the
#' different steps needed to perform the scan from inputted data.
#'
#' @param edge.Prob a `edgeProb` object
#' @param n.scans TO WRIte
#'
#' @return a `scanList` object (S3 class) containing:
#' \itemize{
#'   \item{`raw.scan.list`: a list of raw binary adjacency matrix shaped like the
#'   `Adj` contained in `edge.Prob`, considered directed in the
#'   algorithm}
#'   \item{`theoretical.scan.list`: a list of binary adjacency matrix, where all ties were
#'   observed but the `mode` has been applied}
#'   \item{`scanList.type`: character scalar. `generate_scan` sets it to
#'   "theoretical", `sample_from_scan` will set it to "empirical" and append the
#'   empirical matrix}
#'   \item{`Adj`: `Adj` data contained in `edge.Prob`}
#'   \item{`total_scan`: `total_scan` data contained in `edge.Prob`}
#'   \item{`n.scans`: inputted `n.scans`}
#'   \item{`mode`: `mode` data contained in `edge.Prob`}
#'   \item{`Adj.subfun`: `Adj.subfun` data contained in `edge.Prob`}
#'   \item{`edge.Prob`: `edge.Prob$P` (only the probability matrix) data}
#'   \item{`use.snPackMat`: logical, if scans should be `snPackMat` objects or
#'   regular matrices}
#' }
#'
#' @noRd
generate_scanList <- function(edge.Prob,n.scans){
  raw.scanList <- draw_raw_scanList(edge.Prob = edge.Prob,n.scans = n.scans)
  scanList <- apply_mode(raw.scanList = raw.scanList,mode = edge.Prob$mode)
  attr(scanList,"attrs") <-
    list(
      scanList.type = "theoretical",
      raw.scanList = raw.scanList,
      Adj = edge.Prob$Adj,
      samp.effort = edge.Prob$samp.effort,
      n.scans = n.scans,
      mode = edge.Prob$mode,
      Adj.subfun = edge.Prob$Adj.subfun,
      edge.Prob = edge.Prob$P # in here it is not a edgeProb object anymore, to avoid storing redundant variables,
    )
  class(scanList) <- c("theoretical","scanList")
  scanList
}

#'  TO WRITE
#'
#' @param scan.list TO WRITE
#' @param exp.design TO WRITE
#'
#' @return TO WRITE
#' @noRd
generate_empiscanList <- function(scan.list,exp.design) {
  empiscanList <- exp.design$FUN.seq(scan.list)
  attrs(empiscanList,"scanList.type") <- "empirical"
  attrs(empiscanList,"theoretical.scanList") <- without_attrs(scan.list)
  class(empiscanList)<- c("empirical","scanList")
  empiscanList
}

#'  TO WRITE
#'
#' @param edge.Prob TO WRITE
#' @param n.scans TO WRITE
#'
#' @return TO WRITE
#' @noRd
draw_raw_scanList <- function(edge.Prob,n.scans) {
  sL <- vapply(
    1:n.scans,
    \(s) {
      stats::rbinom(edge.Prob$P,1L,edge.Prob$P)
    },edge.Prob$Adj
  )
  class(sL) <- "scanList"
  sL
}

# scanList tools ------------------------------------------------------------------------------

#'  TO WRITE
#'
#' @param scanList TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @examples
#' # TO WRITE
get_attrs <- function(scanList) {
  attr(scanList,"attrs")
}

#'  TO WRITE
#'
#' @param scanList TO WRITE
#'
#' @return TO WRITE
#' @noRd
without_attrs <- function(scanList) {
  attr(scanList,"attrs") <- NULL
  scanList
}

#'  TO WRITE
#'
#' @param scanList TO WRITE
#' @param a TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @examples
#' # TO WRITE
attrs <- function(scanList,a = NULL) {
  if (is.null(a)) return(get_attrs(scanList))
  get_attrs(scanList)[[a]]
}

#'  TO WRITE
#'
#' @param x TO WRITE
#' @param which TO WRITE
#' @param value TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @examples
#' # TO WRITE
`attrs<-` <- function(x,which,value) {
  new <- get_attrs(x)
  new[[which]] <- value
  attr(x,"attrs") <- new
  x
}

copy_attrs_to <- function(from,to) {
  if (!inherits(to,"scanList")) {class(to) <- class(from)}
  attr(to,"attrs") <- attrs(from)
  to
}

#' TO WRITE
#'
#' @param sL TO WRITE
#' @param .f TO WRITE
#' @param ... TO WRITE
#' @param USE.NAMES TO WRITE
#'
#' @return TO WRITE
#' @keywords internal
#'
#' @examples
#' # TO WRITE
sLapply <- function(sL,.f,...,USE.NAMES = TRUE) {
  vapply(
    X = 1:(dim(sL)[3]),
    FUN = function(x) .f(sL[,,x]),
    FUN.VALUE = sL[,,1],
    ... = ...,
    USE.NAMES = USE.NAMES
  )
}

#' Print method for `scanList` objects
#' @export
#' @noRd
print.scanList <- function(x,...) {
  print.default(without_attrs(x))
  cat("\n\nHidden attributes:",names(get_attrs(x)))
}

#' transpose method for `scanList` objects
#' @export
#' @noRd
t.scanList <- function(x) {
  aperm(x,c(2,1,3))
}


#' rbind method for `scanList` objects
#' @importFrom abind abind
#' @export
#' @noRd
rbind.scanList <- function(...,deparse.level = 1) {
  do.call(abind::abind,list(...,along = 3))
}
