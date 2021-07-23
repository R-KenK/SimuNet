#' Title
#'
#' @param Adj
#' @param samp.effort
#' @param n.scans
#' @param exp.design
#' @param edge.Prob
#' @param mode
#'
#' @return
#' @export
#'
#' @examples
simunet <-
  function(Adj = NULL,samp.effort = NULL,n.scans = NULL,exp.design = NULL,edge.Prob = NULL,
           mode = c("directed","undirected","max","min","upper","lower","plus","vector")
  ){
    edge.Prob <- determine_edgeProb(Adj = Adj,samp.effort = samp.effort,edge.Prob = edge.Prob)

    scanList <- generate_scanList(edge.Prob = edge.Prob,n.scans = n.scans)

    perform_exp(scan_list = scanList,exp.design = exp.design)
  }

