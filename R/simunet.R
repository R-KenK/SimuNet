#' Perform social network simulations
#'
#' @param Adj integer matrix, the reference adjacency matrix to base edge probabilities on. Users
#'   can either:
#'   * import their own adjacency matrix
#'   * or rely on the [`import_from_asnr()`][import_from_asnr()] function to interact with networks
#'   from the [Animal Social Network Repository](https://github.com/bansallab/asnr)
#' @param samp.effort integer scalar, the sampling effort, or number of scans, that led to obtaining
#'   of Adj
#' @param mode character scalar, specifies what igraph network `mode` should be used to convert the
#'   supplied matrix. Ignored if `sampling.param` is provided. Possible values are:
#'   * `"directed"` (the default): for non-symmetrical adjacency matrix where `Adj[i,j]` doesn't
#'   have the same meaning as `Adj[j,i]`
#'   * `"undirected"`: same as `"max"`
#'   * `"upper"`: undirected matrix focusing only on the upper triangle of `Adj` (relying on
#'   `upper.tri`). Either `"upper"` or `"lower"` could be favor if only one of `Adj[i,j]` and
#'   `Adj[j,i]` should be randomized
#'   * `"lower"`: undirected matrix focusing only on the lower triangle of `Adj` (relying on
#'   `lower.tri`)
#'   * `"max"`: from a `"directed"` randomization process (both `Adj[i,j]` and `Adj[j,i]` will be
#'   drawn at each scan), `max(Adj[i,j],Adj[j,i])` will be kept for both
#'   * `"min"`: from a `"directed"` randomization process (both `Adj[i,j]` and `Adj[j,i]` will be
#'   drawn at each scan), `min(Adj[i,j],Adj[j,i])` will be kept for both
#'   * `"plus."`:  from a `"directed"` randomization process (both `Adj[i,j]` and `Adj[j,i]` will be
#'   drawn at each scan), `Adj[i,j] + Adj[j,i]` will be kept for both
#'   * `"vector"`: experimental. To consider adjacency matrices as flat vectors to be randomized.
#'   Relevance unexplored yet.
#'   * See details in the relevant [`igraph`][igraph::graph_from_adjacency_matrix()] package
#'   documentation.
#' @param n.scans integer scalar, number of scans to generate in the simulation
#' @param exp.design `expDesign` object (cf. [`design_exp()`][design_exp()] function, or
#'   `?design_exp`), that consists in a sequence of experimental manipulations (functions) to
#'   perform on the (theoretical) scanList to be simulated to obtain a empirical scanList.
#' @param ... additional arguments to be passed to the function. Specifically, this is used at the
#'   moment to pass more than one `expDesign` object to run multiple experiments on a given
#'   theoretical `scanList`. This cause the returned object to be a list of empirical `scanList`,
#'   i.e. a `sLlist` object.
#' @param edge.Prob optional. An `edgeProb` object (cf.
#'   [`generate_edgeProb()`][generate_edgeProb()]) that consists in the edge presence probability
#'   matrix at each scan. The probability matrix is drawn from a beta distribution determined via
#'   Bayesian inference, from `Adj` and `samp.effort`.
#'
#'   `edgeProb` object are actually lists that
#'   contain the following components:
#'   * `P`: the edge presence probability matrix
#'   * `Adj`: the inputted adjacency matrix
#'   * `samp.effort`: the inputted sampling effort
#'   * `mode`: the inputted igraph `mode`
#'   * `Adj.subfun`: a matrix function, determined from the igraph `mode` (cf. ), that return a
#'   logical matrix with `TRUE` values only for matrix cells relevant to the igraph `mode.` e.g.
#'   only the upper triangle for `mode = "upper"`
#'
#' @return a `scanList` object, primarily a 3 dimensional array representing the (binary) adjacency
#'   matrices (coded within the first two dimensions of the 3D-array) obtained at each simulated
#'   scan (coded as the 3rd dimension of the 3D-array), and a list of attributes.
#'
#'   `scanList` objects have this common structure:
#'   * the 3D-array where the first 2 dimensions are the adjacency matrices (with the node names
#'   from `Adj`) and the 3rd dimension is the simulated scan number
#'   * an attribute named `attrs`: a list of objects - see `attrs` as a flat list of attributes -
#'   that are recorded throughout the simulation and subsequent experimental manipulations.
#'   We provide an equivalent to r base's `attr()` function - [`attrs()`][attrs()] - to retrieve
#'   scanList objects' named attributes contained in their `attrs`
#'
#' @export
#'
#'
#' @examples
#' set.seed(42)
#' n <- 5L
#' samp.effort <- 241L
#'
#' # Adjacency matrix import
#' ## random directed adjacency matrix
#' Adj <- sample(1:samp.effort,n * n) |>
#'   matrix(nrow = 5,dimnames = list(letters[1:n],letters[1:n]))
#' diag(Adj) <- 0L
#' Adj
#'
#' ## manual lower adjacency matrix
#' Adj <- c(0, 0, 0, 0,0,
#'          1, 0, 0, 0,0,
#'          2, 3, 0, 0,0,
#'          4, 5, 6, 0,0,
#'          7, 8, 9,10,0) |>
#'   matrix(nrow = 5,byrow = TRUE,dimnames = list(as.character(1:n),as.character(1:n)))
#' Adj
#'
#' ## upper adjacency matrix imported from ASNR (https://github.com/bansallab/asnr)
#' Adj <- import_from_asnr(class = "Mammalia",
#'                         species = "kangaroo_proximity_weighted",
#'                         output = "adjacency",type = "upper")
#' Adj
#'
#' # social network simulations
#' ## theoretical scans
#' sL <- simunet(Adj = Adj,samp.effort = samp.effort,mode = "upper",n.scans = 120L)
#' sL
#' sL |> sum_scans()
#'
#' ## group-scan sampling
#' ### Designing the experiment: setting a constant probability of not observing edges
#' group.scan <- design_exp(customize_sampling(method = "group",sampling = 0.8))
#'
#' ### simulation can be directly run through the simunet() function
#' simunet(Adj = Adj,samp.effort = samp.effort,mode = "upper",n.scans = 120L,
#'         exp.design = group.scan)
#' ### or the experiment can be applied to a theoretical scanList object
#' group.sL <- perform_exp(sL,group.scan)
#' group.sL |> count_nonNA()
#' group.sL |> sum_scans()
#'
#' ## add more scans, perform even focal sampling, then remove the overall most peripheral node
#' foc.peri_removed <- design_exp(function(x) add_scans(x,200),
#'                                customize_sampling(method = "focal",sampling = "even"),
#'                                remove_mostPeripheral
#' )
#' ### or the experiment can be applied to a theoretical scanList object
#' foc.peri_removed.sL <- perform_exp(sL,foc.peri_removed)
#' foc.peri_removed.sL |> count_nonNA()
#' foc.peri_removed.sL |> sum_scans()
simunet <- function(Adj = NULL,
                    samp.effort = NULL,
                    mode = c("directed","undirected","max","min","upper","lower","plus","vector"),
                    n.scans = NULL,
                    exp.design = NULL,
                    ...,
                    edge.Prob = NULL
) {
  mode <- match.arg(mode)

  scan.list <-
    determine_edgeProb(Adj = Adj,
                       mode = mode,
                       samp.effort = samp.effort,
                       edge.Prob = edge.Prob) |>
    generate_scanList(n.scans = n.scans)

  perform_exp(scan.list = scan.list,exp.design = exp.design,...)
}

