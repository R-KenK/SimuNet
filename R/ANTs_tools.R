
# managing the import/export ----------------------------------------------

#' Convert a binary matrix into an edge-list
#'
#' @param mat a binary matrix
#' @param method character scalar, the scan.list to transform into an edge-list. Either "theoretical","group", or "focal"
#'
#' @return a data frame with columns i and j representing the association between node i and j (more or less considered always directed at present)
#' @noRd
convert_mat_to_df <- function(mat,method,l) {
  edge.present <- mat == 1
  indices <- which(edge.present,arr.ind = TRUE)
  data.frame(i = row.names(mat)[indices[,1]],
             j = colnames(mat)[indices[,2]],
             n.scan = rep(l,nrow(indices))
  )
}

#' Convert a scan object into a data frame
#' Helper function to ensure compatibility with the ANTs package
#'
#' @param scan a scan object
#' @param method character scalar, the scan.list to transform into an edge-list. Either "theoretical","group", or "focal"
#'
#' @return a data frame with columns i and j representing the association between node i and j (more or less considered always directed at present)
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
#' simu <- simu_scan(Adj,total_scan,"all")
#' convert_scan_to_df(simu)
convert_scan_to_df <- function(scan,method = c("theoretical","group","focal")) {
  if (!inherits(scan,"empiScan")) {
    if (length(method) == 1 && method != "theoretical") {
      warning("Requested method not available for non empirical scans")
    }
    method <- "theoretical"
  } else {
    method <- match.arg(method)
  }
  scan.list.name <- paste0(method,".scan.list")

  if (method == "focal") {
    scan <- convert_to_focal_line(scan)
  }

  rbind_lapply(
    seq_along(scan[[scan.list.name]]),
    function(l) {
      s <- scan[[scan.list.name]][[l]]
      convert_mat_to_df(s,method,l)
    }
  )
}

convert_to_focal_line <- function(scan) {
  scan$focal.scan.list <-
    lapply(
      seq_along(scan$focal.scan.list),
      function(s) {
        m <- scan$focal.scan.list[[s]]
        foc <- scan$focal$focal[s]
        m <- m + t(m)
        m[-foc,] <- NA
        m
      }
    )
  scan
}
