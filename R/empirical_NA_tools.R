#' Replace `NA`s with known values when they can be determined
#' Internal use. When the `"min"` or `"max"` algorithm is chosen to convert directed networks into undirected one, some `NA`s can be deduced from their transposed values (cf. below).
#'
#' @param empirical.scan.list a list of binary adjacency matrix potentially containing `NA`s, in which some `NA`s' value can actually be resolved depending on the chosen mode: `scan[i,j] = NA & scan[j,i] = 0 => scan[i,j] = 0` when `mode = "min"`, and `scan[i,j] = NA & scan[j,i] = 1 => scan[i,j] = 1` when `mode = "max"`
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
#'
#' @return a list of binary adjacency matrix where resolvable `NA`s are changed depending on the chosen mode: `scan[i,j] = NA & scan[j,i] = 0 => scan[i,j] = 0` when `mode = "min"`, and `scan[i,j] = NA & scan[j,i] = 1 => scan[i,j] = 1` when `mode = "max"`
#' @noRd
resolve_NA <- function(empirical.scan.list,mode = c("directed", "undirected", "max","min", "upper", "lower", "plus","vector")){
  switch(mode,
         "undirected" = , # as in `igraph`, consider this mode to be the same as `max`
         "max" = lapply(
           empirical.scan.list,
           function(scan) {
             ifelse(
               is.na(scan) | is.na(t(scan)),
               ifelse(
                 scan == 1 | t(scan) == 1,
                 1L,
                 NA),
               scan
             )
           }
         ),
         "min" = lapply(
           empirical.scan.list,
           function(scan) {
             ifelse(
               is.na(scan) | is.na(t(scan)),
               ifelse(
                 scan == 0 | t(scan) == 0,
                 0L,
                 NA),
               scan
             )
           }
         ),
         "plus" = {  # WHAT DOES THIS MEAN FOR BINARY SCANS?
           lapply(
             empirical.scan.list,
             function(scan) {
               ifelse(
                 is.na(scan) | is.na(t(scan)),
                 NA,
                 scan
               )
             }
           )
         },
         "directed" = ,
         "upper" = ,
         "lower" =  ,
         "vector" = empirical.scan.list
  )
}
#' Count observed edges (non-`NA`s) for each edge in list of scans
#' Internal use. Used to determine the sampling effort across all scans performed
#'
#' @param scan.list a list of binary adjacency matrices, where an unobserved dyad (whether it is theoretically 0 or 1) is `NA`
#'
#' @return an integer matrix representing the sampling effort for each dyad
#' @noRd
count_nonNA <- function(scan.list,Adj.subfun) {
  scan.sampled <- Reduce("+",
                         lapply(
                           scan.list,
                           function(scan) {
                             ifelse(!is.na(scan),1L,0L) # counting part of the algorithm
                           }
                         )
  )
  scan.sampled[Adj.subfun(scan.sampled)] <- 0L
  scan.sampled
}

#' Replace NAs by zeros in vectors/matrices
#'
#' @param X a vector or matrix
#'
#' @return similarly dimensioned vector or matrix with zeros instead of NAs
#' @noRd
zero_NA<- function(X){
  ifelse(!is.na(X),X,0L)
}

#' Replace NAs by a set value in vectors/matrices
#'
#' @param X a vector or matrix
#' @param value a numeric scalar
#'
#' @return similarly dimensioned vector or matrix with `value` instead of NAs
#' @noRd
replace_NA<- function(X,value){
  ifelse(!is.na(X),X,value)
}

#' Compare elements of a matrix with its transposed
#' wrapper using `zero_NA`
#'
#' @param X a numeric matrix
#' @param comp.fun a function to compare X with t(X), in this order. Default is superior or equal
#' @param manage_NA.fun a function to replace `NA`s in `X` with a value. Default is `zero_NA` to replace them with zeros (recycled code)
#'
#' @return a logical matrix
#' @noRd
compare_with_transposed<- function(X,comp.fun = `>=`,manage_NA.fun = zero_NA){
  comp.fun(manage_NA.fun(X),manage_NA.fun(t(X)))
}
