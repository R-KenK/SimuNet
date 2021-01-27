# WIP: Wrapper to simulate a single scan ------------------------------------------

#' Simulate a single scan, theoretical or empirical
#'
#' @param Adj square integers matrix of occurrences of dyads. Empirical
#'   adjacency matrix for the simulation to inspire its internal probabilities
#'   from.
#' @param total_scan integer, sampling effort. Note that 1/total_scan should be
#'   relatively small, increasingly small with increasing precision. Optional if
#'   using a `presenceProb` object.
#' @param scans.to.do Optional. Only required if no `sampling.param` inputted (in which case it is ignored). Either:
#'  \itemize{
#'   \item{an integer vector included in `1:total_scan` of the scans to perform}
#'   \item{the special case `"all"` (default) sets `scans.to.do` to `1:total_scan` and set the simulation to perform all the scans}
#' }
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
#' @return depending if `sampling.param` are provided, either:
#'  \itemize{
#'    \item{a `scan` object (S3 class) containing:
#'      \item{`raw.scan.list`: a list of raw binary adjacency matrix shaped like the
#'      `Adj` contained in `presence.prob`, considered directed in the
#'      algorithm}
#'      \item{`theoretical.scan.list`: a list of binary adjacency matrix, where all
#'      ties were observed but the `mode` has been applied}
#'      \item{`scan.type`: character scalar. `generate_scan` sets it to
#'      "theoretical", `sample_from_scan` will set it to "empirical" and append
#'      the empirical matrix}
#'      \item{`Adj`: `Adj` data contained in `presence.prob`}
#'      \item{`total_scan`: `total_scan` data contained in `presence.prob`}
#'      \item{`scans.to.do`: inputted `scans.to.do`}
#'      \item{`mode`: `mode` data contained in `presence.prob`}
#'      \item{`weighted`: logical, at this stage can only be `TRUE` if `mode =
#'      plus` (some edges can become `2`)}
#'      \item{`Adj.subfun`: `Adj.subfun` data contained in `presence.prob`}
#'      \item{`presence.prob`: `presence.prob$P` (only the probability matrix)
#'      data contained in `presence.prob`}
#'    }
#'    \item{an `empiScan` object (S3 class) containing:
#'      \item{`raw.scan.list`: a list of binary adjacency matrix, considered directed
#'      in the algorithm}
#'      \item{`theoretical.scan.list`: a list of binary adjacency matrix, where all
#'      ties were observed but the `mode` has been applied}
#'      \item{`scan.type`: set to `"empirical"` by `sample_from_scan`}
#'      \item{`method`: from inputted `sampling.param`}
#'      \item{`scans.to.do`: from inputted `sampling.param`}
#'      \item{`group.scan.list`: a list of adjacency matrix, where the observation
#'      probability `obs.prob` from `sampling.param` of each dyad has been
#'      applied. `NULL` if `method = "focal"`}
#'      \item{`focal.scan.list`: a list of adjacency matrix, where only the selected
#'      `focal` from `sampling.param` is visible. `NULL` if `method = "group"`}
#'      \item{`Adj`: `Adj` data contained in `presence.prob`}
#'      \item{`total_scan`: `total_scan` data contained in `presence.prob`}
#'      \item{`mode`: `mode` data contained in `presence.prob`}
#'      \item{`weighted`: logical, at this stage can only be `TRUE` if `mode =
#'      plus` (some edges can become `2`)}
#'      \item{`Adj.subfun`: `Adj.subfun` data contained in `presence.prob`}
#'      \item{`presence.prob`: `presence.prob$P` (only the probability matrix)
#'      data contained in `presence.prob`}
#'    }
#'  }
#'
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
#' # by default will simulate a single directed theoretical scan
#' simu_scan(Adj,total_scan)
#'
#' # but other mode can be used, as well as more scans:
#' theo.scans <- simu_scan(Adj,total_scan,mode = "min",scans.to.do = 1:3)
#' theo.scans
#'
#' # `summary()` methods can be used to sum up scans contained in a `scan` object
#' # (and can be accessed in the `$theoretical.sum` of the object returned by
#' # `summary()`)
#' summary(theo.scans)
#'
#'
#' # Users can generate sampling parameters through `simu_samplingParam` to use in
#' # `simu_scan`, in which case `Adj` and `total_scan` arguments are optional.
#' para.group.constant<- simu_samplingParam(Adj,total_scan,mode =
#'                                          "min",group.scan_param = 0.42)
#' group.constant.scans <- simu_scan(sampling.param = para.group.constant)
#' group.constant.scans
#'
#' # `summary()` methods can be used to sum up scans contained in an `empiScan`
#' # object (and can be accessed in the `$group.sum` or `$focal.sum` of the object
#' # returned by `summary()`)
#' summary(group.constant.scans)
#'
#' # Users can also define functions to use trait- or network- based sampling
#' # biases for group-scan sampling (cf. `?simu_samplingParam`)
#' obs.prob.trait.bias_fun<- function(i,j,Adj) {i+j} # comparable to a dyad-trait-based bias
#' para.group.trait.bias<- simu_samplingParam(Adj,total_scan,mode ="directed",
#'                                            group.scan_param = obs.prob.trait.bias_fun,
#'                                            scans.to.do = 1:15)
#' para.group.net.bias<- simu_samplingParam(Adj,total_scan,mode =
#'                                          "max",group.scan_param = function(i,j,Adj) {Adj*Adj})
#'
#' simu_scan(sampling.param = para.group.trait.bias)
#' simu_scan(sampling.param = para.group.net.bias)
#'
#' # or for biases regarding which focals to draw for focal-scan sampling (cf.
#' # `?simu_samplingParam`)
#' focal.trait.bias_fun<- function(n,Adj) {1:n} # comparable to a dyad-trait-based bias
#' para.focal.trait.bias<- simu_samplingParam(Adj,
#'                                            total_scan,mode = "directed",
#'                                            focal.scan_param = focal.trait.bias_fun,
#'                                            scans.to.do = 1:10)
#' para.focal.net.bias<- simu_samplingParam(Adj,total_scan,mode = "max",
#'                                          focal.scan_param = function(n,Adj) {colSums(Adj*Adj)},
#'                                          scans.to.do = 20)
#'
#' simu_scan(Adj,total_scan,sampling.param = para.focal.trait.bias)
#' simu_scan(Adj,total_scan,sampling.param = para.focal.net.bias)
#'
#' # Users can also perform both group-scan and focal-scan sampling methods at the
#' # same time by inputting both a `group.scan_param` and a `focal.scan_param` when
#' # creating their sampling parameter object
#' para.both.no.bias <- simu_samplingParam(Adj,
#'                                         total_scan,mode = "min",
#'                                         group.scan_param = 0.9,
#'                                         focal.scan_param = "even",
#'                                         scans.to.do = "all")
#' simu.both <- simu_scan(sampling.param = para.both.no.bias)
#' simu.both
#' summary(simu.both)
simu_scan <-
  function(Adj = NULL,
           total_scan = NULL,scans.to.do = NULL,
           mode = c("directed","undirected","max","min","upper","lower","plus","vector"),
           sampling.param = NULL) {
    if (!is.null(sampling.param)) {
      if (is.samplingParam(sampling.param)) {
        if (!is.null(sampling.param$obs.prob)) {
          Adj <- sampling.param$obs.pro$Adj
          total_scan <- sampling.param$obs.pro$total_scan
        }
        if (!is.null(sampling.param$focal)) {
          Adj <- sampling.param$focal$focal.list$Adj
          total_scan <- sampling.param$focal$focal.list$total_scan
        }
        method <- sampling.param$method
        mode <- sampling.param$mode
        obs.prob <- sampling.param$obs.prob
        focal <- sampling.param$focal
        if(!is.null(scans.to.do)) {
          if(scans.to.do != sampling.param$scans.to.do) {
            warning("Inputted `scans.to.do` different from the one stored in `sampling.param`.\nThe latter has been used.")
          }
        }
        scans.to.do <- sampling.param$scans.to.do
      } else {
        stop(
          "Please provide a valid `samplingParam` object. You can use simu_samplingParam() to create one."
        )
      }
    } else {
      mode <- match.arg(mode)
      if (is.null(scans.to.do)) {scans.to.do <- 1} # in case nothing is provided at all, still only performs a single scan
    }

    # Check if `Adj` has been passed as a `presenceProb` object or not, retrieve
    # other variable otherwise if needed
    if (!is.presenceProb(Adj)) {
      Adj.subfun <- determine_Adj.subfun(mode)
      mode(Adj) <- "integer"  # helps storing derived objects or elements as integer matrix to save on memory allocation
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
    scan <- generate_scan(presence.prob = presence.prob,scans.to.do = scans.to.do)

    # output either the theoretical scan or applies sample_from_scan
    if (is.null(sampling.param)) {
      scan
    } else {
      generate_empiScan(scan = scan, sampling.param = sampling.param)
    }
  }
