# expDesign object functions ------------------------------------------------------------------

#' Title
#'
#' @param scan_list
#' @param exp.design
#'
#' @return
#' @export
#'
#' @examples
perform_exp <- function(scan.list,exp.design = NULL){
  if (!inherits(scan.list,"scanList")) {stop("scan.list inputted is not a scanList object.")}
  if (is.null(exp.design)) {
    return(scan.list)
  } else if(!inherits(exp.design,"expDesign")) {stop("exp.design inputted is not a expDesign object.")}
  generate_empiscanList(scan.list,exp.design)
}

#' Title
#'
#' @param ...
#' @param .dir
#'
#' @return
#' @export
#'
#' @examples
design_exp <- function(...,.dir = c("forward", "backward")) {
  .dir <- match.arg(.dir)
  FUN.seq <- purrr::compose(... = ...,.dir = .dir)
  generate_expDesign(FUN.seq = FUN.seq)
}

#' Title
#'
#' @param FUN.seq
#'
#' @return
#' @export
#'
#' @examples
generate_expDesign <- function(FUN.seq) {
  expD <-
    list(
      FUN.seq = FUN.seq
    )
  class(expD) <- "expDesign"
  expD
}


