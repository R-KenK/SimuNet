# graphml file import -----------------------------------------------------

#' TO WRITE
#'
#' @param path TO WRITE
#' @param output TO WRITE
#'
#' @return TO WRITE
#' @importFrom igraph read_graph
#' @importFrom igraph get.adjacency
#' @noRd
import_from_graphml <- function(path,
                                output = c("graph","adjacency"),
                                type = c("both", "upper", "lower")) {
  output <- match.arg(output)
  G <- igraph::read_graph(path,format = "graphml")

  switch(output,
         "graph" = G,
         "adjacency" = {
           Adj <- igraph::get.adjacency(G,type = type,attr = "weight",sparse = FALSE)
           add_Adj_names(G,Adj)
         }
  )
}

#' Add names to Adj from graph G
#'
#' @param G an igraph network object
#' @param Adj an adjacency matrix
#'
#' @return an Adjacency matrix which col and row names have been filled with the names contained in G
#' @importFrom igraph vertex_attr
#' @noRd
add_Adj_names <- function(G, Adj) {
  if(!is.null(igraph::vertex_attr(G,"name")) | !is.null(igraph::vertex_attr(G,"id"))) {
    if(is.null(igraph::vertex_attr(G,"name"))) {
      rownames(Adj) <- igraph::vertex_attr(G,"id")
    } else {
      rownames(Adj) <- igraph::vertex_attr(G,"name")
    }
    colnames(Adj) <- rownames(Adj)
  } else {
    rownames(Adj) <- as.character(1:nrow(Adj));colnames(Adj)<- as.character(1:ncol(Adj))
  }
  Adj
}