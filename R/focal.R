#' Generator for `focal` objects
#'
#' @param focal.list a `focalList` object
#' @param scans.to.do either:
#'  \itemize{
#'   \item{an integer vector included in `1:total_scan` of the scans to perform}
#'   \item{the special case `"all"` (default) sets `scans.to.do` to `1:total_scan` and set the simulation to perform all the scans}
#' }
#'
#' @return an `focal` object (S3 class) containing:
#' \itemize{
#'   \item{`focal`: named integer vector representing the _index_ of the node(s) to sample at scan number(s) `scans.to.do`}
#'   \item{`focal.list`: inputted `focal.list` object}
#'   \item{`scans.to.do`:  inputted `scans.to.do`. If `"all"` was inputted, is set to `1:total_scan`}
#' }
#' @export
#'
#' @examples
#'
#' set.seed(42)
#'
#' n<- 5;nodes<- letters[1:n];total_scan<- 42;
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- sample(0:total_scan,n*(n-1),replace = TRUE)
#' Adj
#'
#' focal.list<- generate_focalList(Adj,total_scan,focal.prob_fun = "even")
#' generate_focal(focal.list,5)
#' generate_focal(focal.list,3:10)
#' generate_focal(focal.list,"all")
#' # below returns an error if `scans.to.do` is incompatible with `total_scan`
#' # generate_focal(focal.list,50)
#' # generate_focal(focal.list,40:43)
#'
generate_focal<- function(focal.list,
                          scans.to.do){
  if (length(scans.to.do) == 1) {
    if (scans.to.do == "all") {
      scans.to.do<- 1:focal.list$total_scan
    }
  }
  if (any(scans.to.do > focal.list$total_scan)) {stop("Inputted `scans.to.do` involves number(s) higher than `total_scan`.")}

  focal<- list(
    focal = focal.list$focals[scans.to.do],
    focal.list = focal.list,
    scans.to.do = scans.to.do
  )
  class(focal)<- "focal"
  focal
}

#' Print method for `focal` objects
#' @export
#' @noRd
print.focal<- function(x,...){
  cat(
    paste0(
      ifelse(!is.null(names(x$focal)),"\n  node name: ",""),names(x$focal), # won't display anything if nodes don't have names
      "\n node index: ",x$focal,
      "\nscan number: ",x$scans.to.do," out of ",x$focal.list$total_scan,"\n"
    ),sep = ""
  )
}

#' Test if object if a `focal` object
#'
#' @param x an object to test.
#'
#' @return logical, `TRUE` if the inputted object is a `focal` object.
#'
#' @noRd
is.focal<- function(x){
  inherits(x,"focal")
}

