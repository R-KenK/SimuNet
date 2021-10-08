
#####################################################
# WIP for potential interaction with the ANTs package
#####################################################

# managing the import/export ----------------------------------------------

#' Convert a binary matrix into an edge-list
#'
#' @param mat a binary matrix
#' @param method character scalar, the scan.list to transform into an edge-list. Either
#'   "theoretical","group", or "focal"
#'
#' @return a data frame with columns i and j representing the association between node i and j (more
#'   or less considered always directed at present)
#' @noRd
convert_mat_to_df <- function(mat,method,l) {
  edge.present <- mat == 1
  indices <- which(edge.present,arr.ind = TRUE)
  data.frame(i = row.names(mat)[indices[,1]],
             j = colnames(mat)[indices[,2]],
             n.scan = rep(l,nrow(indices))
  )
}
