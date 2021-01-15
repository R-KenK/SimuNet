# WIP: Wrapper to simulate a single scan ------------------------------------------

#' Simulate a single scan, theoretical or empirical
#'
#' @param Adj square integers matrix of occurrences of dyads. Empirical
#'   adjacency matrix for the simulation to inspire its internal probabilities
#'   from.
#' @param total_scan integer, sampling effort. Note that 1/total_scan should be
#'   relatively small, increasingly small with increasing precision. Optional if
#'   using a `presenceProb` object.
#' @param mode Character scalar, specifies what type of igraph network `mode`
#'   should be used to convert the supplied matrix. Ignored if `sampling.param` is
#'   provided. Possible values are:
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
#' @param sampling.param Optional if a theoretical scan is needed. Otherwise
#'   a `samplingParam` object (S3 class) containing:
#' \itemize{
#'   \item{`method`: a character scalar:
#'     \item{`"group"`: an empirical group-scan is performed, in which some
#'     dyads can be missed (`NA` in the empirical adjacency matrix)}
#'     \item{`"focal"`: an empirical focal-scan is performed, in which only the
#'     row and column of the chosen focal is recorded (others are `NA` in the
#'     empirical adjacency matrix)}
#'     \item{or `"both"`: both methods are done in parallel from the same
#'     theoretical scan}
#'   }
#'   \item{mode: inputted `mode`}
#'   \item{obs.prob: inputted `obs.prob`}
#'   \item{focal: inputted `focal`}
#' }
#'
#' @return
#' @export
#'
#' @examples
#' set.seed(42)
#'
#' n<- 5;nodes<- letters[1:n];total_scan<- 42;
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- sample(0:total_scan,n*(n-1),replace = TRUE)
#' Adj
#'
#' # by default will simulate a directed theoretical scan
#' simu_scan(Adj,total_scan)
#'
#' # Users can generate sampling parameters through `simu_samplingParam` to use
#' in `simu_scan`
#' para.group.constant<- simu_samplingParam(Adj,total_scan,mode =
#'                                            "max",group.scan_param = 0.42)
#' simu_scan(Adj,total_scan,para.group.constant)
#'
#' # Users can also define functions to use trait- or network- based sampling
#' # biases (cf. ?simu_samplingParam)
#' obs.prob.trait.bias.ij<- function(i,j,Adj) {i+j} # comparable to a dyad-trait-based bias
#' obs.prob.net.bias.Adj<- function(i,j,Adj) {Adj*Adj} # comparable to a network-based bias
#'
#' simu_scan(Adj,total_scan,obs.prob.trait.bias.ij)
#'
#' generate_obsProb(Adj,"directed",obs.prob_fun = user_function.ij)
#' generate_obsProb(Adj,"directed",obs.prob_fun = user_function.Adj)
simu_scan <-
  function(Adj = NULL,
           total_scan = NULL,
           mode = c("directed","undirected","max","min","upper","lower","plus","vector"),
           sampling.param = NULL) {
    if (!is.null(sampling.param)) {
      if (!is.samplingParam(sampling.param)) {
        stop(
          "Please provide a valid `samplingParam` object. You can use simu_samplingParam() for instance"
        )
      }
      mode <- sampling.param$mode
    } else {
      mode <- match.arg(mode)
    }

    # Check if `Adj` has been passed as a `presenceProb` object or not, retrieve
    # other variable otherwise if needed
    if (!is.presenceProb(Adj)) {
      Adj.subfun <- determine_Adj.subfun(mode)
      presence.prob <-
        generate_presenceProb(
          Adj = Adj,
          total_scan = total_scan,
          mode = mode,
          Adj.subfun = Adj.subfun
        )
    } else {
      presence.prob <- Adj
    }

    # use the class "scan" object generator to
    scan <- generate_scan(presence.prob)

    # output either the theoretical scan or applies sample_from_scan
    if (is.null(sampling.param)) {
      scan
    } else {
      if (is.null(sampling.param)) {
        sampling.param <-
          generate_samplingParam(
            method = method,
            mode = mode,
            obs.prob = obs.prob,
            focal.list = focal.list
          ) # should return an error in case of missing parameters given the chosen method
      }
      generate_empiScan(scan = scan, sampling.param = sampling.param)
    }
  }
