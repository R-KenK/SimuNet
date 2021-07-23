# expDesign object functions ------------------------------------------------------------------

#' TO WRITE
#'
#' @param scan.list TO WRITE
#' @param exp.design TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @examples
#' # TO WRITE
perform_exp <- function(scan.list,exp.design = NULL){
  if (!inherits(scan.list,"scanList")) {stop("scan.list inputted is not a scanList object.")}
  if (is.null(exp.design)) {
    return(scan.list)
  } else if(!inherits(exp.design,"expDesign")) {stop("exp.design inputted is not a expDesign object.")}
  generate_empiscanList(scan.list,exp.design)
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


