# Adjacency mode-related functions -----------------------------------------

#' Determine the right adjacency matrix subsetting function according to the chosen igraph mode
#'
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. See also the weighted argument, the interpretation depends on that too. Possible values are: directed, undirected, upper, lower, max, min, plus. See details \link[igraph]{graph_from_adjacency_matrix}. Added also a vector mode.
#'
#' @return a subsetting function among `non.diagonal`, `upper.tri`, `lower.tri` or `function(V) {rep(TRUE,length(V))}`
  #' @noRd
determine_Adj.subfun<- function(mode){
  switch(mode,
         "directed" = ,
         "undirected" = ,
         "max" = ,
         "min" = ,
         "plus" = non.diagonal,
         "upper" = upper.tri,
         "lower" =  lower.tri,
         "vector" = function(V) {rep(TRUE,length(V))}  # kept for now but I don't know if that is useful
  )
}

#' Make Adjacency fit the selected mode
#' Internal use.
#'
#' @param raw.scan.list a list of (directed by default) binary adjacency matrix, pre-theoretical, to which the chosen `mode` is to be applied
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
#'
#' @return an list of adjacency matrix fitting the chosen `mode` is to be applied
#' @noRd
apply_mode<- function(raw.scan.list,mode = c("directed", "undirected", "max","min", "upper", "lower", "plus","vector")){
  switch(mode,
         "undirected" = , # as in `igraph`, consider this mode to be the same as `max`
         "max" = lapply(
           raw.scan.list,
           function(scan) {
             ifelse(scan >= t(scan),scan,t(scan))
           }
         ),
         "min" = lapply(
           raw.scan.list,
           function(scan) {
             ifelse(scan <= t(scan),scan,t(scan))
           }
         ),
         "plus" = {  # WHAT DOES THIS MEAN FOR BINARY SCANS?
           lapply(
             raw.scan.list,
             function(scan) {
               not.na <- !is.na(scan) & !is.na(t(scan))
               ifelse(not.na,scan+t(scan),NA)
             }
           )
         },
         "directed" = ,
         "upper" = ,
         "lower" =  ,
         "vector" = raw.scan.list
  )
}

