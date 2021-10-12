# Adjacency mode-related functions -----------------------------------------

#' Determine the right adjacency matrix subsetting function according to the chosen igraph mode
#'
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. Default
#'   here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus.
#'   Added vector too. See details
#'   [`graph_from_adjacency_matrix()`][igraph::graph_from_adjacency_matrix()].
#'
#' @return a subsetting function among `non.diagonal`, `upper.tri`, `lower.tri` or `function(V)
#'   {rep(TRUE,length(V))}`
#'
#' @keywords internal
#'
#' @export
determine_Adj.subfun <- function(mode){
  switch(mode,
         "directed" = ,
         "undirected" = ,
         "max" = ,
         "min" = ,
         "plus" = non.diagonal,
         "upper" = upper.tri,
         "lower" =  lower.tri,
         "vector" = function(V) {rep(TRUE,length(V))} # kept for now but I don't know if that is useful
  )
}

#' Make Adjacency fit the selected mode
#' Internal use.
#'
#' @param raw.scanList a list of (directed by default) binary adjacency matrix, pre-theoretical, to
#'   which the chosen `mode` is to be applied
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. Default
#'   here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus.
#'   Added vector too. See details
#'   [`graph_from_adjacency_matrix()`][igraph::graph_from_adjacency_matrix()].
#'
#' @return an list of adjacency matrix fitting the chosen `mode` is to be applied
#'
#' @export
#'
#' @keywords internal
apply_mode <-
  function(raw.scanList,
           mode = c("directed", "undirected", "max","min", "upper", "lower", "plus","vector")){
    switch(mode,
           "undirected" = , # as in `igraph`, consider this mode to be the same as `maraw.scanList`
           "maraw.scanList" = ifelse(raw.scanList >= t(raw.scanList),raw.scanList,t(raw.scanList)),
           "min" = ifelse(raw.scanList <= t(raw.scanList),raw.scanList,t(raw.scanList)),
           "plus" = {  # WHAT DOES THIS MEAN FOR BINARY SCANS?
             ifelse(!is.na(raw.scanList) & !is.na(t(raw.scanList)),raw.scanList,t(raw.scanList))
           },
           "directed" = ,
           "upper" = ,
           "lower" =  ,
           "vector" = raw.scanList
    )
  }

