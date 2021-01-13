#' Generator for `focal` objects
#'
#' @param focal.list a `focalList` object
#' @param scan.number an integer between in `1:total_scan`
#'
#' @return an `focal` object (S3 class) containing:
#' \itemize{
#'   \item{focal}{named integer scalar representing the _index_ of the node to sample at scan number `scan.number`}
#'   \item{focal.list}{inputted `focal.list` object}
#'   \item{scan.number}{inputted `scan.number`}
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
#' focal.list<- generate_focal.list(Adj,total_scan,focal.prob_fun = "even")
#' generate_focal(focal.list,5)
#' # generate_focal(focal.list,50) # returns an error if `scan.number` is incompatible with `total_scan`
#'
generate_focal<- function(focal.list,
                          scan.number){
  if(scan.number > focal.list$total_scan) {stop("Inputted `scan.number` higher than `total_scan`.")}

  focal<- list(
    focal = focal.list$focals[scan.number],
    focal.list = focal.list,
    scan.number = scan.number
  )
  class(focal)<- "focal"
  focal
}

#' Print method for `focal` objects
#' @export
print.focal<- function(focal,...){
  cat(
    paste0(
      ifelse(!is.null(names(focal$focal)),paste0("  node name: ",names(focal$focal)),""),
      "\n node index: ",focal$focal,
      "\nscan number: ",focal$scan.number," out of ",focal$focal.list$total_scan
    )
  )
}

#' Test if object if a `focal` object
#'
#' @param focal an object to test.
#'
#' @return logical, TRUE if the inputted object is a `focal` object.
#'
#' @noRd
is.focal<- function(focal){
  inherits(focal,"focal")
}

