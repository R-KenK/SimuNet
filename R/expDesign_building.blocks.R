# expDesign "building-block" functions --------------------------------------------------------
. <- NULL

# scanList manipulations ----
#'  TO WRITE
#'
#' @param scan.list TO WRITE
#' @param which TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @examples
#' # TO WRITE
sum_scans <- function(scan.list,which = c("auto","theoretical","raw")) {
  which <- match.arg(which)
  sL.ori <- scan.list
  switch(which,
         "auto" = {},
         "theoretical" = scan.list <- attrs(scan.list,"theoretical.scanList"),
         "raw" = scan.list <- attrs(scan.list,"raw.scanList.type")
  )
  summed <- rowSums(scan.list,na.rm = TRUE,dims = 2L)
  summed <- copy_attrs_to(sL.ori,summed)
  attrs(summed,"summed.scanList") <- without_attrs(sL.ori)
  class(summed) <- c("sum",class(scan.list))
  summed
}

#'  TO WRITE
#'
#' @param summed TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @examples
#' # TO WRITE
scale_scans <- function(summed) {
  scaled <- summed / attrs(summed,"Adj")
  scaled <-
    scaled |> is.nan() |> ifelse(0,scaled)
  attr(scaled,"attrs") <- get_attrs(summed)
  class(scaled) <- c("scaled",class(summed))
  scaled
}

#'  TO WRITE
#'
#' @param scan.list TO WRITE
#' @param empirical TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @examples
#' # TO WRITE
count_na <- function(scan.list,empirical = TRUE) {
  scan.list |> is.na() |> ifelse(1L,0L) |> sum_scans()
}

#'  TO WRITE
#'
#' @param scan.list TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @examples
#' # TO WRITE
remove_mostPeripheral <- function(scan.list) {
  scan.list |>
    sum_scans() |>
    igraph::graph.adjacency(weighted = TRUE) |>
    igraph::eigen_centrality(directed = TRUE) %>%
    .$vector |>
    which.min() %>%
    {scan.list[-c(.),-c(.),]}
}
