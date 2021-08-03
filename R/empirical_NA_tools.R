#' Replace `NA`s with known values when they can be determined
#' Internal use. When the `"min"` or `"max"` algorithm is chosen to convert directed networks into
#' undirected one, some `NA`s can be deduced from their transposed values (cf. below).
#'
#' @param scan.list a list of binary adjacency matrix potentially containing `NA`s, in which some
#'   `NA`s' value can actually be resolved depending on the chosen mode:
#'    \itemize{
#'      \item{`scan[i,j] = NA & scan[j,i] = 0 => scan[i,j] = 0` when `mode = "min"`}
#'      \item{`scan[i,j] = NA & scan[j,i] = 1 => scan[i,j] = 1` when `mode = "max"`}
#'    }
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. Default
#'   here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus.
#'   Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
#'
#' @return a list of binary adjacency matrix where resolvable `NA`s are changed depending on the
#'   chosen mode:
#'    \itemize{
#'      \item{`scan[i,j] = NA & scan[j,i] = 0 => scan[i,j] = 0` when `mode = "min"`}
#'      \item{`scan[i,j] = NA & scan[j,i] = 1 => scan[i,j] = 1` when `mode = "max"`}
#'    }
#'
#' @keywords internal
#' @export
resolve_NA <- function(scan.list,mode = c("directed", "undirected", "max","min", "upper", "lower", "plus","vector")){
  if (length(mode) > 1L) mode <- attrs(scan.list,"mode")
  else mode <- match.arg(mode)
  switch(mode,
         "undirected" = , # as in `igraph`, consider this mode to be the same as `max`
         "max" = scan.list |>
           {\(x) ifelse(is.na(x) | is.na(t(x)),ifelse(x == 1 | t(x) == 1,1L,NA),x)}(),
         "min" = scan.list |>
           {\(x) ifelse(is.na(x) | is.na(t(x)),ifelse(x == 0 | t(x) == 0,0L,NA),x)}(),
         "plus" = scan.list |>
           {\(x) ifelse(is.na(x) | is.na(t(x)),as.integer(NA),x)}(),  # WHAT DOES THIS MEAN FOR BINARY SCANS?
         "directed" = ,
         "upper" = ,
         "lower" =  ,
         "vector" = scan.list
  )
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
