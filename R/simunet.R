#' TO WRITE
#'
#' @param Adj TO WRITE
#' @param samp.effort TO WRITE
#' @param mode TO WRITE
#' @param n.scans TO WRITE
#' @param exp.design TO WRITE
#' @param ... TO WRITE
#' @param edge.Prob TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @examples
#' # TO WRITE
simunet <- function(Adj = NULL,
                    samp.effort = NULL,
                    mode = c("directed","undirected","max","min","upper","lower","plus","vector"),
                    n.scans = NULL,
                    exp.design = NULL,
                    ...,
                    edge.Prob = NULL
) {
  mode <- match.arg(mode)

  scan.list <-
    determine_edgeProb(Adj = Adj,
                       mode = mode,
                       samp.effort = samp.effort,
                       edge.Prob = edge.Prob) |>
    generate_scanList(n.scans = n.scans)

  perform_exp(scan.list = scan.list,exp.design = exp.design,...)
}

