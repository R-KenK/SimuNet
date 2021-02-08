# Empirical sampling related functions ------------------------------------

#' Generator for `empiScan` objects
#'
#' @param scan a `scan` object with `scan.type = "theoretical"`
#' @param sampling.param a `samplingParam` object containing:
#' \itemize{
#'   \item{method}{inputted `method`}
#'   \item{mode}{inputted `mode`}
#'   \item{scans.to.do: inputted `scans.to.do`}
#'   \item{obs.prob}{inputted `obs.prob`}
#'   \item{focal}{inputted `focal`}
#' }
#'
#' @return an `empiScan` object (S3 class), inheriting from `scan` object,
#'   containing:
#' \itemize{
#'   \item{`raw.scan.list`: a list of binary adjacency matrix, considered directed in
#'   the algorithm}
#'   \item{`theoretical.scan.list`: a list of binary adjacency matrix, where all ties
#'   were observed but the `mode` has been applied}
#'   \item{`scan.type`: set to `"empirical"` by `sample_from_scan`}
#'   \item{`method`: from inputted `sampling.param`}
#'   \item{`scans.to.do`: from inputted `sampling.param`}
#'   \item{`group.scan.list`: a list of adjacency matrix, where the observation
#'   probability `obs.prob` from `sampling.param` of each dyad has been applied.
#'   `NULL` if `method = "focal"`}
#'   \item{`focal.scan.list`: a list of adjacency matrix, where only the selected
#'   `focal` from `sampling.param` is visible. `NULL` if `method = "group"`}
#'   \item{`Adj`: `Adj` data contained in `presence.prob`}
#'   \item{`total_scan`: `total_scan` data contained in `presence.prob`}
#'   \item{`mode`: `mode` data contained in `presence.prob`}
#'   \item{`Adj.subfun`: `Adj.subfun` data contained in `presence.prob`}
#'   \item{`presence.prob`: `presence.prob$P` (only the probability matrix) data
#'   contained in `presence.prob`}
#' }
#' @noRd
generate_empiScan <- function(scan, sampling.param) {
  scan$scan.type <- "empirical"
  method <- sampling.param$method
  mode <- scan$mode
  scans.to.do <- sampling.param$scans.to.do
  if (!is.null(sampling.param$obs.prob)) {
    group.scan.list <-
      sample_from_scan(scan = scan,
                       sampling.param = sampling.param,
                       method = "group")
    group.scan.na.resolved <- resolve_NA(empirical.scan.list = group.scan.list,mode = mode)
    group.scan.sum <- sum_scan.list(group.scan.na.resolved)
    obs.prob <- sampling.param$obs.prob
  } else {
    group.scan.sum <- group.scan.list <- obs.prob <- NULL
  }
  if (!is.null(sampling.param$focal)) {
    focal.scan.list <-
      sample_from_scan(scan = scan,
                       sampling.param = sampling.param,
                       method = "focal")
    focal.scan.na.resolved <- resolve_NA(empirical.scan.list = focal.scan.list,mode = mode)
    focal.scan.sum <- sum_scan.list(focal.scan.na.resolved)
    focal <- sampling.param$focal
  } else {
    focal.scan.sum <- focal.scan.list <- focal <- NULL
  }
  scan <- list(
    raw.scan.list = scan$raw.scan.list,
    theoretical.scan.list = scan$theoretical.scan.list,
    group.scan.list = group.scan.list,
    obs.prob = obs.prob,
    focal.scan.list = focal.scan.list,
    focal = focal,
    scan.type = "empirical",
    Adj = scan$Adj,
    total_scan = scan$total_scan,
    scans.to.do = scan$scans.to.do,
    mode = mode,
    Adj.subfun = scan$Adj.subfun,
    presence.prob = scan$presence.prob # in here it is not a presenceProb object anymore, to avoid storing redundant variables
  )
  class(scan) <- c("empiScan", "scan")
  scan
}

#' Print method for `scan` objects
#' @export
#' @noRd
print.empiScan <- function(x, ...) {
  scans.to.do <- explicit_scans.to.do(x)
  n <- length(scans.to.do)
  # print the general simulation infos
  if (n >= 15) {
    scans.to.do <- paste(do.call(paste,as.list(scans.to.do)[1:5]),
                         "...",
                         do.call(paste,as.list(scans.to.do)[(n-4):n])
    )
  } else {
    scans.to.do <- scans.to.do
  }
  cat("\nScan(s) performed: ",scans.to.do,sep=" ")
  cat("\n\nScan type: theoretical, mode: ", x$mode, "\n\n", sep = "")

  # display the theoretical scans
  shorten_list.to.print(x$theoretical.scan.list)

  # display the optional group scans
  if (!is.null(x$group.scan.list)) {
    cat(
      "\n\nScan type: group scan, mode: ",
      x$mode,
      "\nobs.prob type: ",
      x$obs.prob$obs.prob_type,
      "\n\n",
      sep = ""
    )
    shorten_list.to.print(x$group.scan.list)
  }

  # display the optional focal scans
  if (!is.null(x$focal.scan.list)) {
    cat(
      "\n\nScan type: focal scan, mode: ",
      x$mode,
      "\nfocal.prob type: ",
      x$focal$focal.list$focal.prob_type,
      "\n\n",
      sep = ""
    )
    shorten_list.to.print(x$focal.scan.list)
  }
}

#' Summary method for `empiScan` objects
#' @export
#' @noRd
summary.empiScan <- function(object,...) {
  scans.to.do <- explicit_scans.to.do(object)
  theoretical.sum <- sum_scan.list(object$theoretical.scan.list)
  theoretical.sampled <- sum_scan.sampled(object,method = "theoretical")
  theoretical.scaled <- theoretical.sum/ifelse(theoretical.sampled != 0,theoretical.sampled,1)
  mode <- object$mode
  scans.to.do <- object$scans.to.do
  if (!is.null(object$group.scan.list)) {
    group.scan.na.resolved <- resolve_NA(empirical.scan.list = object$group.scan.list,mode = mode)
    group.sum <- sum_scan.list(group.scan.na.resolved)
    group.sampled <- sum_scan.sampled(object,method = "group")
    group.scaled <- group.sum/ifelse(group.sampled != 0,group.sampled,1)
    obs.prob_type <- object$obs.prob$obs.prob_type
  } else {
    obs.prob_type <- group.scaled <- group.sum <- group.sampled <- group.scan.list <- obs.prob <- NULL
  }
  if (!is.null(object$focal.scan.list)) {
    focal.scan.na.resolved <- resolve_NA(empirical.scan.list = object$focal.scan.list,mode = mode)
    focal.sum <- sum_scan.list(focal.scan.na.resolved)
    focal.sampled <- sum_scan.sampled(object,method = "focal")
    focal.scaled <- focal.sum/ifelse(focal.sampled != 0,focal.sampled,1)
    focal.prob_type <- object$focal$focal.list$focal.prob_type
  } else {
    focal.prob_type <- focal.scaled <- focal.sum <- focal.sampled <- focal.scan.list <- focal <- NULL
  }
  scan.summary <- list(
    theoretical.sum = theoretical.sum,
    theoretical.sampled = theoretical.sampled,
    theoretical.scaled = theoretical.scaled,
    group.sum = group.sum,
    group.sampled = group.sampled,
    group.scaled = group.scaled,
    obs.prob_type = obs.prob_type,
    focal.sum = focal.sum,
    focal.sampled = focal.sampled,
    focal.scaled = focal.scaled,
    focal.prob_type = focal.prob_type,
    scans.to.do = scans.to.do,
    mode = mode#,
    # more things can be added here
  )
  class(scan.summary) <- c("summary.empiScan","summary.scan")
  scan.summary
}

#' Print method for `summary.empiScan` objects
#' @importFrom Matrix Matrix
#' @importFrom Matrix printSpMatrix
#' @export
#' @noRd
print.summary.empiScan<- function(x,...){
  print.summary.scan(x,...)
  if (!is.null(x$group.sum)) {
    cat("Group-scan sampling method weighted adjacency matrix:\n")
    use_printSpMatrix(x$group.sum)
    cat(paste0("\nobtained after the following per-edge sampling matrix:", "\n\n"))
    use_printSpMatrix(x$group.sampled)
  }
  if (!is.null(x$focal.sum)) {
    cat("Focal-scan sampling method weighted adjacency matrix:\n")
    use_printSpMatrix(x$focal.sum)
    cat(paste0("\nobtained after the following per-edge sampling matrix:", "\n\n"))
    use_printSpMatrix(x$focal.sampled)
  }
}

#' Plot method for `empiScan` objects
#' @export
#' @noRd
plot.empiScan<- function(x,method = c("both","theoretical","group","focal"),...){
  method <- match.arg(method)
  x.summary <- summary.empiScan(x)
  plot.summary.empiScan(x.summary,method = method,...)
}

#' Plot method for `scan` objects
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph plot.igraph
#' @importFrom igraph layout_with_fr
#' @importFrom graphics layout
#' @importFrom graphics par
#' @export
#' @noRd
plot.summary.empiScan <- function(x,
                                 method = c("both","theoretical","group","focal"),
                                 scaled = FALSE,
                                 vertex.size = NULL,
                                 vertex.size.mul = nrow(x$theoretical.sum)/2,
                                 vertex.size.min = 3,
                                 vertex.size.fun = compute.strength,
                                 layout = NULL,
                                 ...){
  method <- match.arg(method)

  switch(method,
         "theoretical" = {
           plot.summary.scan(x,
                             scaled = scaled,
                             vertex.size = vertex.size,
                             vertex.size.mul = vertex.size.mul,
                             vertex.size.min = vertex.size.min,
                             vertex.size.fun = vertex.size.fun,
                             layout = layout,
                             ...)
         },
         "group" = ,
         "focal" = {
           plot_empirical(x = x,
                          method = method,scaled = scaled,
                          vertex.size = vertex.size,
                          vertex.size.mul = vertex.size.mul,
                          vertex.size.min = vertex.size.min,
                          vertex.size.fun = vertex.size.fun,
                          layout = layout,
                          ...)
         },
         "both" = {
           # graphics::par(mar=c(5,1,2,1))
           graphics::layout(matrix(c(0,1,1,0,2,2,3,3), 2, 4, byrow = TRUE))
           if (scaled) {theoretical.adj <- x$theoretical.scaled} else {theoretical.adj <- x$theoretical.sum}
           if (is.null(layout) | is.function(layout)) {
             x.G <- igraph::graph_from_adjacency_matrix(theoretical.adj,mode = x$mode,weighted = TRUE)
             layout <- layout(x.G)
           }
           plot.summary.scan(x,
                             scaled = scaled,
                             vertex.size = vertex.size,
                             vertex.size.mul = vertex.size.mul,
                             vertex.size.min = vertex.size.min,
                             vertex.size.fun = vertex.size.fun,
                             layout = layout,
                             ...)
           # par(mar=c(1,2,1,2))
           plot_empirical(x = x,
                          method = "group",scaled = scaled,
                          vertex.size = vertex.size,
                          vertex.size.mul = vertex.size.mul,
                          vertex.size.min = vertex.size.min,
                          vertex.size.fun = vertex.size.fun,
                          layout = layout,
                          ...)
           plot_empirical(x = x,
                          method = "focal",scaled = scaled,
                          vertex.size = vertex.size,
                          vertex.size.mul = vertex.size.mul,
                          vertex.size.min = vertex.size.min,
                          vertex.size.fun = vertex.size.fun,
                          layout = layout,
                          ...)
           # par(mar = c(5,4,4,2))
           graphics::layout(matrix(c(1,1,1), 1, 1, byrow = TRUE))
           graphics::par(mar=c(5,4,4,2))
         }
  )
}

#' Plot method for `scan` objects
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph plot.igraph
#' @noRd
plot_empirical<- function(x,
                          method = c("group","focal"),
                          scaled = FALSE,
                          vertex.size = NULL,
                          vertex.size.mul = nrow(x$theoretical.sum)/2,
                          vertex.size.min = 3,
                          vertex.size.fun = compute.strength,
                          ...){
  method <- match.arg(method)
  if (scaled) {empirical.adj <- x[[paste0(method,".scaled")]]} else {empirical.adj <- x[[paste0(method,".sum")]]}
  total_scan <- max(x$theoretical.sampled,na.rm = TRUE)
  max.adj <- max(empirical.adj,na.rm = TRUE)
  x.G <- igraph::graph_from_adjacency_matrix(empirical.adj,mode = x$mode,weighted = TRUE)
  if (is.null(vertex.size)) {
    vertex.size <- vertex.size.fun(x.G,mode = x$mode)
    vertex.size <- (vertex.size.mul * (vertex.size / max(vertex.size))) + vertex.size.min
  }
  switch(method,
         "group" = {
           X.prob_type <- x$obs.prob_type
           X.sentence <- " type of edge observation probability "
         },
         "focal" = {
           X.prob_type <- x$focal.prob_type
           X.sentence <- " type of focal sampling probability "
         }
  )

  igraph::plot.igraph(x.G,...,main = paste0("Scan type: ",method," scan sampling"),
                      sub = paste0("relying on a ",X.prob_type,X.sentence,"(mode = \"", x$mode,"\")"),
                      vertex.size = vertex.size
  )
}


#' Test if object if a `empiScan` object
#'
#' @param scan an object to test.
#'
#' @return logical, `TRUE` if the inputted object is a `empiScan` object.
#'
#' @noRd
is.empiScan <- function(scan) {
  inherits(scan, "empiScan")
}

#' Empirically sample from a theoretical scan
#'
#'
#' @param scan a `scan` object with `scan.type = "theoretical"`
#' @param sampling.param a `samplingParam` object containing:
#' \itemize{
#'   \item{`method`: inputted `method`}
#'   \item{`mode`: inputted `mode`}
#'   \item{`scans.to.do`: inputted `scans.to.do`}
#'   \item{`obs.prob`: inputted `obs.prob`}
#'   \item{`focal`: inputted `focal`}
#' }
#'
#' @return either:
#' \itemize{
#'   \item{a list of binary adjacency matrix with potentially some dyads not observed
#'   (turned into `NA`)}
#'   \item{a list of binary adjacency matrix with dyads not in the `focal$focal` row and
#'   column turned into `NA`}
#' }
#'
#' @noRd
sample_from_scan <- function(scan, sampling.param, method) {
  # according to the empirical method chosen, applies some "empirical" missed
  # observation via the correct sampling method
  switch(
    method,
    "group" = group_sample(scan = scan, obs.prob = sampling.param$obs.prob),
    "focal" = focal_sample(scan = scan, focal = sampling.param$focal)
  )
}

#' Perform a group scan sampling
#' Internal use. Presently, the sampling is performed on the _theoretical_ scan
#'
#' @param scan a `scan` object with `scan.type = "theoretical"`
#' @param obs.prob an `obsProb` object
#'
#' @importFrom stats rbinom
#'
#' @return a list of binary adjacency matrix with potentially some dyads not observed
#'   (turned into `NA`)
#' @noRd
group_sample <- function(scan, obs.prob) {
  # set the sampled scan to be like the `raw.scan.list` one (i.e. this is a binary
  # _directed_ adjacency matrix)
  observed <-
    scan$theoretical.scan.list
  # subset obs.prob like the adjacency matrix (i.e. triangular matrix) into a
  # vector of observation probabilities
  obs.P <-
    obs.prob$P[scan$Adj.subfun(obs.prob$P)]
  lapply(
    observed,
    function(s) {
      # draw the missed observation
      missed <-
        stats::rbinom(length(obs.P), 1, obs.P) == 0
      # set the missed observation to `NA`
      s <- unpack_snPackMat(s)
      s[scan$Adj.subfun(s)][missed] <-
        NA
      s # standard
      # Matrix::pack(as.matrix(s)) # Matrix.packed
      generate_snPackMat(M = s,Adj.subfun = scan$Adj.subfun,mode = scan$mode) # Matrix.packed
    }
  )
}

#' Perform a focal scan sampling
#' Internal use.
#'
#' @param scan a `scan` object with `scan.type = "theoretical"`
#' @param focal an `focal` object
#'
#' @return a list of binary adjacency matrix with dyads not in the `focal$focal` row and
#'   column turned into `NA`
#' @noRd
focal_sample <- function(scan, focal) {
  # set the sampled scan to be like the `raw.scan.list` one (i.e. this is a binary
  # _directed_ adjacency matrix)
  observed <-
    scan$theoretical.scan.list
  lapply(
    seq_along(observed),
    function(s) {
      obs <- observed[[s]]
      obs <- unpack_snPackMat(obs)

      foc <- focal$focal[s]
      # set other rows and columns than those of the `focal` to `NA`
      obs[-foc,-foc] <- NA
      obs[!scan$Adj.subfun(obs)] <- 0L
      obs # standard
      # Matrix::pack(as.matrix(obs)) # Matrix.packed
      generate_snPackMat(M = obs,Adj.subfun = scan$Adj.subfun,mode = scan$mode) # Matrix.packed
    }
  )
}
