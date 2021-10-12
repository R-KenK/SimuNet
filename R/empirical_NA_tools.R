#' Replace `NA`s with known values when they can be determined
#' Internal use. When the `"min"` or `"max"` algorithm is chosen to convert directed networks into
#' undirected one, some `NA`s can be deduced from their transposed values (cf. below).
#'
#' @param scan.list a `scanList` object. See objects returned by [`simunet()`][simunet()]
#' @param mode character scalar, See [`simunet()`][simunet()]
#'
#'
#' @return a list of binary adjacency matrix where resolvable `NA`s are changed depending on the
#'   chosen mode:
#'    * `scan[i,j] = NA & scan[j,i] = 0 => scan[i,j] = 0` when `mode = "min"`
#'    * `scan[i,j] = NA & scan[j,i] = 1 => scan[i,j] = 1` when `mode = "max"`
#'
#' @keywords internal
#' @export
resolve_NA <- function(scan.list,mode = c("directed", "undirected", "max","min", "upper", "lower", "plus","vector")){
  if (length(mode) > 1L) mode <- attrs(scan.list,"mode")
  else mode <- match.arg(mode)
  switch(mode,
         "undirected" = , # as in `igraph`, consider this mode to be the same as `max`
         "max" = ifelse(is.na(scan.list) | is.na(t(scan.list)),
                        ifelse(scan.list == 1 | t(scan.list) == 1,1L,NA),scan.list),
         "min" = ifelse(is.na(scan.list) | is.na(t(scan.list)),
                        ifelse(scan.list == 0 | t(scan.list) == 0,0L,NA),scan.list),
         "plus" = ifelse(is.na(scan.list) | is.na(t(scan.list)),
                         as.integer(NA),scan.list),  # WHAT DOES THIS MEAN FOR BINARY SCANS?
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
