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
    theoretical.scan.sum = scan$theoretical.scan.sum,
    group.scan.list = group.scan.list,
    group.scan.sum = group.scan.sum,
    obs.prob = obs.prob,
    focal.scan.list = focal.scan.list,
    focal.scan.sum = focal.scan.sum,
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
  if (length(x$scans.to.do) == 1) {
    if (x$scans.to.do == "all") {
      scans.to.do <- 1:x$total_scan
    } else {
      scans.to.do <- x$scans.to.do
    }
  } else {
    scans.to.do <- x$scans.to.do
  }
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
  cat("\nScans performed: ",scans.to.do,sep=" ")
  cat("\n\nScan type: theoretical, mode: ", x$mode, "\n\n", sep = "")

  # display the theoretical scans
  if (n >= 10) {
    print.default(x$theoretical.scan.list[1:3],...)
    cat("... (",n-5," more scans)\n\n\n",sep = "")
    cat("[[",n-1,"]]\n",sep = "")
    print.default(x$theoretical.scan.list[[(n-1)]],...)
    cat("[[",n,"]]\n",sep = "")
    print.default(x$theoretical.scan.list[[n]],...)
  } else {
    print.default(x$theoretical.scan.list,...)
  }

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
    if (n >= 10) {
      print.default(x$group.scan.list[1:3],...)
      cat("... (",n-5," more scans)\n\n\n",sep = "")
      cat("[[",n-1,"]]\n",sep = "")
      print.default(x$group.scan.list[[(n-1)]],...)
      cat("[[",n,"]]\n",sep = "")
      print.default(x$group.scan.list[[n]],...)
    } else {
      print.default(x$group.scan.list, ...)
    }
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
    if (n >= 10) {
      print.default(x$focal.scan.list[1:3],...)
      cat("... (",n-5," more scans)\n\n\n",sep = "")
      cat("[[",n-1,"]]\n",sep = "")
      print.default(x$focal.scan.list[[(n-1)]],...)
      cat("[[",n,"]]\n",sep = "")
      print.default(x$focal.scan.list[[n]],...)
    } else {
      print.default(x$focal.scan.list, ...)
    }
  }
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

#' Plot method for `scan` objects
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph plot.igraph
#' @export
#' @noRd
plot.empiScan <-
  function(x, ...) {
    # Need a way to make it print two in case of method = "both"
    x <-
      igraph::graph_from_adjacency_matrix(x$theoretical.scan.list,
                                          mode = x$mode,
                                          weighted = x$weighted)
    igraph::plot.igraph(x, ...)
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
      s[scan$Adj.subfun(s)][missed] <-
        NA
      s
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
      foc <- focal$focal[s]
      # set other rows and columns than those of the `focal` to `NA`
      obs[-foc,-foc] <-
        NA
      obs
    }
  )
}
