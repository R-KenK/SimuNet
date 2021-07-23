#' TO WRITE
#'
#' @param Adj TO WRITE
#' @param samp.effort TO WRITE
#' @param n.scans TO WRITE
#' @param exp.design TO WRITE
#' @param edge.Prob TO WRITE
#' @param mode TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @examples
#' # TO WRITE
simunet <-
  function(Adj = NULL,samp.effort = NULL,n.scans = NULL,exp.design = NULL,edge.Prob = NULL,
           mode = c("directed","undirected","max","min","upper","lower","plus","vector")
  ){
    edge.Prob <- determine_edgeProb(Adj = Adj,samp.effort = samp.effort,edge.Prob = edge.Prob)

    scan.list <- generate_scanList(edge.Prob = edge.Prob,n.scans = n.scans)

    perform_exp(scan.list = scan.list,exp.design = exp.design)
  }

