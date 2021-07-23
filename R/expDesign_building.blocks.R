# expDesign "building-block" functions --------------------------------------------------------
#' Title
#'
#' @param scan.list
#' @param obs.P
#'
#' @return
#' @export
#'
#' @examples
group_samp.list <- function(scan.list,obs.P){
  vapply(
    1:dim(scan.list)[3],
    function(s) {
      scan.list[,,s] %>% {ifelse(rbinom(.,1L,obs.P) == 1,.,NA)}
    },scan.list[,,1]
  )
}

#' Title
#'
#' @param scan.list
#'
#' @return
#' @export
#'
#' @examples
sum_scans <- function(scan.list,which = c("auto","theoretical","raw")) {
  which <- match.arg(which)
  sL.ori <- scan.list
  switch(which,
         "auto" = {},
         "theoretical" = scan.list <- attrs(scan.list,"theoretical.scanList"),
         "raw" = scan.list <- attrs(scan.list,"raw.scanList.type")
  )
  summed <- rowSums(scan.list,na.rm = TRUE,dims = 2L)
  attr(summed,"attrs") <- get_attrs(sL.ori)
  class(summed) <- c("sum",class(scan.list))
  summed
}

#' Title
#'
#' @param scan.list
#'
#' @return
#' @export
#'
#' @examples
scale_scans <- function(summed) {
  scaled <- summed / attrs(summed,"Adj")
  scaled <-
    scaled |> is.nan() |> ifelse(0,scaled)
  attr(scaled,"attrs") <- get_attrs(summed)
  class(scaled) <- c("scaled",class(summed))
  scaled
}

#' Title
#'
#' @param scan.list
#'
#' @return
#' @export
#'
#' @examples
count_na <- function(scan.list,empirical = TRUE) {
  scan.list |> is.na() |> ifelse(1L,0L) |> sum_scans()
}

#' Title
#'
#' @param scan.list
#'
#' @return
#' @export
#'
#' @examples
remove_mostPeripheral <- function(scan.list) {
  scan.list |>
    sum_scans() |>
    igraph::graph.adjacency(weighted = TRUE) |>
    igraph::eigen_centrality(directed = TRUE) %>%
    .$vector |>
    which.min() %>%
    {scan.list[-c(.),-c(.),]}
}

#' Title
#'
#' @param array.3D
#'
#' @return
#' @export
#'
#' @examples
array2matList <- function(array.3D) {
  vapply(1:dim(array.3D)[3],\(s) array.3D[,,s],FUN.VALUE = array.3D[,,1])
}

#' Title
#'
#' @param method
#' @param obsProb
#'
#' @return
#' @export
#'
#' @examples
customize_sampling <- function(method = c("group","focal"),obsProb) {
  method <- match.arg(method)
  switch(method,
         "group" = \(scan.list) group_samp.list(scan.list,obs.P = obsProb),
         "focal" = stop("not implemented")
  )
}
