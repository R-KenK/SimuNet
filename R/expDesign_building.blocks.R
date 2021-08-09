# expDesign "building-block" functions --------------------------------------------------------
. <- NULL

# scanList manipulations ----

## sum ----

#'  TO WRITE
#'
#' @param scan.list TO WRITE
#' @param ... TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @examples
#' # TO WRITE
sum_scans <- function(scan.list,...) {
  UseMethod("sum_scans")
}


#'  TO WRITE
#' @export
#' @noRd
sum_scans.scanList <- function(scan.list,which = c("auto","theoretical","raw"),...) {
  which <- match.arg(which)
  sf <- attrs(scan.list,"Adj.subfun")
  sL.ori <- scan.list
  switch(which,
         "auto" = {},
         "theoretical" = scan.list <- attrs(scan.list,"theoretical.scanList"),
         "raw" = scan.list <- attrs(scan.list,"raw.scanList.type")
  )
  summed <- rowSums(scan.list,na.rm = TRUE,dims = 2L)
  summed <- copy_attrs_to(sL.ori,summed)
  attrs(summed,"summed.scanList") <- without_attrs(sL.ori)
  attrs(summed,"sampled") <- scan.list |> count_nonNA()
  class(summed) <- c("sum",class(scan.list))
  summed
}

#'  TO WRITE
#' @export
#' @noRd
sum_scans.sLlist <- function(scan.list,which = c("auto","theoretical","raw"),...) {
  sLlapply(scan.list,sum_scans,which = which)
}

#'  TO WRITE
#' @export
#' @noRd
sum_scans.empirical <- function(scan.list,which = c("auto","theoretical","raw"),...) {
  NextMethod()
}

## scale ----

#'  TO WRITE
#'
#' @param summed TO WRITE
#' @param ... TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @examples
#' # TO WRITE
scale_scans <- function(summed,...) {
  UseMethod("scale_scans")
}

#'  TO WRITE
#' @export
#' @noRd
scale_scans.sum <- function(summed,...) {
  sf <- attrs(summed,"Adj.subfun")
  sampled <- attrs(summed,"sampled")
  scaled <- summed
  scaled[sf(scaled)] <- summed[sf(summed)] / sampled[sf(sampled)]
  ifelse(!is.infinite(scaled),sampled,NA)
  scaled <- copy_attrs_to(summed,scaled)
  class(scaled) <- c("scaled",class(summed))
  scaled
}

#'  TO WRITE
#' @export
#' @noRd
scale_scans.scanList <- function(summed,...) {
  summed |>
    sum_scans() |>
    scale_scans()
}

#'  TO WRITE
#' @export
#' @noRd
scale_scans.empirical <- function(summed,...) {
  NextMethod()
}

#'  TO WRITE
#' @export
#' @noRd
scale_scans.sLlist <- function(summed,...) {
  sLlapply(summed,scale_scans)
}

## count NAs ----
#'  TO WRITE
#'
#' @param scan.list TO WRITE
#' @param ... TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @examples
#' # TO WRITE
count_NA <- function(scan.list,...) {
  UseMethod("count_NA")
}

#'  TO WRITE
#' @export
#' @noRd
count_NA.scanList <- function(scan.list,empirical = TRUE,...) {
  sf <- attrs(scan.list,"Adj.subfun")
  scan.sampled <- scan.list |> is.na() |> ifelse(1L,0L) %>% copy_attrs_to(from = scan.list)
  scan.sampled <- scan.sampled |> rowSums(na.rm = TRUE,dims = 2L)
  scan.sampled[!sf(scan.sampled)] <- 0L
  scan.sampled
}

#'  TO WRITE
#' @export
#' @noRd
count_NA.empirical <- function(scan.list,empirical = TRUE,...) {
  NextMethod()
}

#'  TO WRITE
#' @export
#' @noRd
count_NA.sLlist <- function(scan.list,empirical = TRUE,...) {
  sLlapply(scan.list,count_NA,empirical = empirical)
}

## count non NAs ----

#' Count observed edges (non-`NA`s) for each edge in list of scans
#' Internal use. Used to determine the sampling effort across all scans performed
#'
#' @param scan.list a list of binary adjacency matrices, where an unobserved dyad (whether it is theoretically 0 or 1) is `NA`
#' @param ... TO WRITE
#'
#' @return an integer matrix representing the sampling effort for each dyad
#'
#' @export
#'
#' @examples
#' # TO WRITE
count_nonNA <- function(scan.list,...) {
  UseMethod("count_nonNA")
}

#'  TO WRITE
#' @export
#' @noRd
count_nonNA.scanList <- function(scan.list,...) {
  sf <- attrs(scan.list,"Adj.subfun")
  scan.sampled <- scan.list |> is.na() |> ifelse(0L,1L) %>% copy_attrs_to(from = scan.list)
  scan.sampled <- scan.sampled |> rowSums(na.rm = TRUE,dims = 2L)
  scan.sampled[!sf(scan.sampled)] <- 0L
  scan.sampled
}

#'  TO WRITE
#' @export
#' @noRd
count_nonNA.empirical <- function(scan.list,...) {
  NextMethod()
}

#'  TO WRITE
#' @export
#' @noRd
count_nonNA.sLlist <- function(scan.list,...) {
  sLlapply(scan.list,count_nonNA)
}

## add scans ----

#'  TO WRITE
#'
#' @param scan.list TO WRITE
#' @param ... TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @examples
#' # TO WRITE
add_scans <- function(scan.list,...) {
  UseMethod("add_scans")
}

#'  TO WRITE
#' @export
#' @noRd
add_scans.scanList <- function(scan.list,new.scans,...) {
  edge.Prob <- reconstruct_edgeProb(scan.list)
  n.scans <- attrs(scan.list,"n.scans")

  new.scan.list <- generate_scanList(edge.Prob = edge.Prob,n.scans = new.scans)

  new.scan.list <- rbind(scan.list,new.scan.list)
  new.scan.list <- copy_attrs_to(scan.list,new.scan.list)

  total.scans <- attrs(new.scan.list,"n.scans") + new.scans
  attr(total.scans,"scans.performed") <-
    c(new.scans,
      if (is.null(attr(n.scans,"scans.performed"))) n.scans
      else attr(n.scans,"scans.performed")
    )
  attrs(new.scan.list,"n.scans") <- total.scans
  new.scan.list
}

#'  TO WRITE
#' @export
#' @noRd
add_scans.empirical <- function(scan.list,new.scans,...) {
  NextMethod()
}

#'  TO WRITE
#' @export
#' @noRd
add_scans.sLlist <- function(scan.list,new.scans,...) {
  sLlapply(scan.list,add_scans,new.scans = new.scans)
}

## remove peripheral individual ----

#'  TO WRITE
#'
#' @param scan.list TO WRITE
#' @param ... TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @examples
#' # TO WRITE
remove_mostPeripheral <- function(scan.list,...) {
  UseMethod("remove_mostPeripheral")
}

#'  TO WRITE
#' @export
#' @noRd
remove_mostPeripheral.scanList <- function(scan.list,...) {
  mode <- attrs(scan.list,"mode")
  directed <- switch(mode,"directed" = TRUE,FALSE)
  scan.list |>
    sum_scans() |>
    igraph::graph.adjacency(weighted = TRUE) |>
    igraph::eigen_centrality(directed = directed) %>%
    .$vector |>
    which.min() %>%
    {scan.list[-c(.),-c(.),]} |>
    copy_attrs_to(from = scan.list)
}

#'  TO WRITE
#' @export
#' @noRd
remove_mostPeripheral.sLlist <- function(scan.list,...) {
 sLlapply(remove_mostPeripheral,scan.list)
}

