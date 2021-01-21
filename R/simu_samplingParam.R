# sampling.param generation and related functions -------------------------------

#' Define the sampling parameters to use in the network simulations
#' To use before calling `simu_scan`. Wrapper for `samplingParam` objects.
#'
#' @param Adj square integers matrix of occurrences of dyads. Empirical
#'   adjacency matrix for the simulation to inspire its internal probabilities
#'   from.
#' @param total_scan integer, sampling effort. Note that 1/total_scan should be
#'   relatively small, increasingly small with increasing precision.
#' @param mode Character scalar, specifies what type of igraph network `mode`
#'   should be used to convert the supplied matrix. Possible values are:
#' \itemize{
#'   \item{`"directed"` (the default): for non-symmetrical adjacency matrix where
#'   `Adj[i,j]` doesn't have the same meaning as `Adj[j,i]`}
#'   \item{`"undirected"`: same as `"max"`}
#'   \item{`"upper"`: undirected matrix focusing only on the upper triangle of
#'   `Adj` (relying on `upper.tri`). Either `"upper"` or `"lower"` could be
#'   favor if only one of `Adj[i,j]` and `Adj[j,i]` should be randomized}
#'   \item{`"lower"`: undirected matrix focusing only on the lower triangle of
#'   `Adj` (relying on `lower.tri`)}
#'   \item{`"max"`: from a `"directed"` randomization process (both `Adj[i,j]`
#'   and `Adj[j,i]` will be drawn at each scan), `max(Adj[i,j],Adj[j,i])` will
#'   be kept for both}
#'   \item{`"min"`: from a `"directed"` randomization process (both `Adj[i,j]`
#'   and `Adj[j,i]` will be drawn at each scan), `min(Adj[i,j],Adj[j,i])` will
#'   be kept for both}
#'   \item{`"plus."`:  from a `"directed"` randomization process (both
#'   `Adj[i,j]` and `Adj[j,i]` will be drawn at each scan), `Adj[i,j] +
#'   Adj[j,i]` will be kept for both}
#'   \item{`"vector"`: experimental. To consider adjacency matrices as flat
#'   vectors to be randomized. Relevance unexplored yet.}
#'   \item{See details \link[igraph]{graph_from_adjacency_matrix}}
#' }
#' @param group.scan_param internally the `obs.prob_fun` variable to be passed
#'   to `generate_obsProb`. Either:
#' \itemize{
#'   \item{a single [0,1] numeric value for all dyad representing their
#'   probability of being sampled or not (`obs.prob_type` will be `"constant"`)}
#'   \item{the string `"random"` if each dyad should have its probability drawn
#'   from a uniform distribution between 0 and 1 (`runif(n*n,0,1)`)}
#'   \item{a user-defined function of (i,j,Adj) that output a probability of
#'   presence for the dyad}
#' }
#' @param focal.scan_param internally the `focal.prob_fun` variable to be passed
#'   to `generate_focalList`. Either:
#' \itemize{
#'   \item{the special case `"even"` tries to even out the `focal.list` as much
#'   as possible before drawing randomly following a uniform distribution}
#'   \item{`"random"`: pick focals following a uniform distribution}
#'   \item{a user-defined function of (n,Adj) that output a weight of being
#'   focal for each node (passed as the `prob` argument to `base::sample`
#'   function)}
#' }
#' @param scans.to.do either:
#'  \itemize{
#'   \item{an integer vector included in `1:total_scan` of the scans to perform}
#'   \item{the special case `"all"` (default) sets `scans.to.do` to
#'   `1:total_scan` and set the simulation to perform all the scans}
#' }
#' @param all.sampled logical, should all individuals be sampled at least once
#'   by design? Default is `TRUE`. Internally passed to `generate_focalList`.
#'   After drawing each node once, the rest of the `focal.list` is  sampled
#'   according to `focal.prob_fun`. Returns an error if `total_scan` is smaller
#'   than the number of nodes.
#'
#' @return a `samplingParam` object (S3 class) containing:
#' \itemize{
#'   \item{method: inputted `method`}
#'   \item{mode: inputted `mode`}
#'   \item{obs.prob: inputted `obs.prob`}
#'   \item{focal: inputted `focal`}
#' }
#' @export
#'
#' @examples
#'
#' set.seed(42)
#'
#' n<- 5;nodes<- as.character(1:n);total_scan<- 42;
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- sample(0:total_scan,n*(n-1),replace = TRUE)
#' Adj
#'
#' simu_samplingParam(Adj,total_scan,focal.scan_param = "even")
#' simu_samplingParam(Adj,total_scan,mode = "max",group.scan_param = 0.42)
#' simu_samplingParam(Adj,total_scan,mode = "min",
#'                    group.scan_param = 0.42,
#'                    focal.scan_param = "random",scans.to.do = 1:4)
simu_samplingParam <-
  function(Adj,
           total_scan,
           mode = c("directed",
                    "undirected",
                    "max",
                    "min",
                    "upper",
                    "lower",
                    "plus",
                    "vector"),
           group.scan_param = NULL,
           focal.scan_param = NULL,
           all.sampled = TRUE,
           scans.to.do = "all") {
    mode <- match.arg(mode)

    method <-
      determine_method(group.scan_param = group.scan_param, focal.scan_param = focal.scan_param)

    switch(
      method,
      "group" = {
        obs.prob <-
          generate_obsProb(Adj = Adj,
                           total_scan = total_scan,
                           mode = mode,
                           obs.prob_fun = group.scan_param)
        focal <- NULL
      },
      "focal" = {
        focal.list <-
          generate_focalList(
            Adj = Adj,
            total_scan = total_scan,
            focal.prob_fun = focal.scan_param,
            all.sampled = all.sampled
          )
        focal <-
          generate_focal(focal.list = focal.list, scans.to.do = scans.to.do)
        obs.prob <- NULL
      },
      "both" = {
        obs.prob <-
          generate_obsProb(Adj = Adj,
                           total_scan = total_scan,
                           mode = mode,
                           obs.prob_fun = group.scan_param)
        focal.list <-
          generate_focalList(
            Adj = Adj,
            total_scan = total_scan,
            focal.prob_fun = focal.scan_param,
            all.sampled = all.sampled
          )
        focal <-
          generate_focal(focal.list = focal.list, scans.to.do = scans.to.do)
      },
      stop("Inputted `method` not recognized.")
    )
    generate_samplingParam(
      method = method,
      mode = mode,
      obs.prob = obs.prob,
      focal = focal,
      scans.to.do = scans.to.do
    )
  }

#' Determine the sampling parameters `method`
#' from inputted parameters
#'
#' @param group.scan_param internally the `obs.prob_fun` variable to be passed
#'   to `generate_obsProb`. Either:
#' \itemize{
#'   \item{a single [0,1] numeric value for all dyad representing their
#'   probability of being sampled or not (`obs.prob_type` will be `"constant"`)}
#'   \item{the string `"random"` if each dyad should have its probability drawn
#'   from a uniform distribution between 0 and 1 (`runif(n*n,0,1)`)}
#'   \item{a user-defined function of (i,j,Adj) that output a probability of
#'   presence for the dyad}
#' }
#' @param focal.scan_param internally the `focal.prob_fun` variable to be passed
#'   to `generate_focalList`. Either:
#' \itemize{
#'   \item{the special case `"even"` tries to even out the `focal.list` as much
#'   as possible before drawing randomly following a uniform distribution}
#'   \item{`"random"`: pick focals following a uniform distribution}
#'   \item{a user-defined function of (n,Adj) that output a weight of being
#'   focal for each node (passed as the `prob` argument to `base::sample`
#'   function)}
#' }
#'
#' @return a character scalar:
#' \itemize{
#'   \item{`"group"`}
#'   \item{`"focal"`}
#'   \item{or `"both"`}
#' }
#'
#' @noRd
determine_method <- function(group.scan_param, focal.scan_param) {
  is.group <-
    !is.null(group.scan_param)
  is.focal <- !is.null(focal.scan_param)
  if (is.group & is.focal) {
    return("both")
  } else {
    if (is.group) {
      return("group")
    }
    if (is.focal) {
      return("focal")
    }
    stop("Please input enough parameters to determine which method to use.")
  }
}
