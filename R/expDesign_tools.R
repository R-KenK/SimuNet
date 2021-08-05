# expDesign object functions ------------------------------------------------------------------

#' TO WRITE
#'
#' @param scan.list TO WRITE
#' @param exp.design TO WRITE
#' @param ... TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @examples
#' # TO WRITE
perform_exp <- function(scan.list,exp.design = NULL,...){
  if (!inherits(scan.list,"scanList")) {stop("scan.list inputted is not a scanList object.")}
  if (is.null(exp.design)) {
    return(scan.list)
  } else if (!inherits(exp.design,"expDesign")) {stop("exp.design inputted is not a expDesign object.")}
  if (missing(...)) generate_empiscanList(scan.list,exp.design)
  else {
    expD.list <- list(exp.design,...)
    sL.list <- lapply(expD.list,\(expD) generate_empiscanList(scan.list = scan.list,exp.design = expD))
    class(sL.list) <- "sLlist"
    sL.list
  }
}

#' TO WRITE
#'
#' @param ... TO WRITE
#' @param .dir TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @importFrom purrr compose
#'
#' @examples
#' # TO WRITE
design_exp <- function(...,.dir = c("forward", "backward")) {
  .dir <- match.arg(.dir)
  FUN.seq <- purrr::compose(... = ...,.dir = .dir)
  generate_expDesign(FUN.seq = FUN.seq)
}

#' TO WRITE
#'
#' @param FUN.seq TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @examples
#' # TO WRITE
generate_expDesign <- function(FUN.seq) {
  expD <-
    list(
      FUN.seq = FUN.seq
    )
  class(expD) <- "expDesign"
  expD
}

# scanList convenience and compatibility functions ----
#'  TO WRITE
#'
#' @param array.3D TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @examples
#' # TO WRITE
array2matList <- function(array.3D) {
  vapply(1:dim(array.3D)[3],\(s) array.3D[,,s],FUN.VALUE = array.3D[,,1])
}
