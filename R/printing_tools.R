# TO RECYCLE:
# #' Wrapper to shorten vectors to print with ellipsis
# #'
# #' @param v a vector to print
# #' @param threshold the length of x above which the vector should be shortened with an ellipsis before being printed
# #' @param before number of elements to display before the ellipsis
# #' @param after number of elements to display after the ellipsis
# #'
# #' @return either v or a shortened version of v to display
# #' @noRd
# shorten_vec.to.print <- function(v,threshold = 15,before = 5,after = 5) {
#   n <- length(v)
#   if (n >= threshold) {
#     paste(do.call(paste,as.list(v)[1:before]),
#           "...",
#           do.call(paste,as.list(v)[(n-after+1):n])
#     )
#   } else {
#     v
#   }
# }

# TO RECYCLE:
# #' Print list element in custom format
# #'
# #' @param l a list
# #' @param i the index of the element to print
# #'
# #' @importFrom Matrix Matrix
# #' @importFrom Matrix printSpMatrix
# #'
# #' @return
# #' @noRd
# print_list_element <- function(l,i) {
#   cat("[[",i,"]]\n",sep = "")
#   # l[[i]] <- Matrix::Matrix(l[[i]],sparse = TRUE)
#   use_printSpMatrix(unpack_snPackMat(l[[i]]))
#   cat("\n")
#   invisible(l)
# }

# TO RECYCLE:
# #' Wrapper for using Matrix package's print method
# #'
# #' @param M a regular matrix
# #' @param ... additional argument for `printSpMatrix`
# #' @param sparse logical. Should a sparse matrix be printed?
# #' @param digits integer. Number of digits to display
# #' @param note.dropping.colnames logical. Should the column names be dropped?
# #' @param align character string. among "right", "left" and "center"
# #'
# #' @return print the matrix M using the Matrix package's print method
# #' @noRd
# use_printSpMatrix <- function(M,...,sparse = TRUE, digits = 3, note.dropping.colnames = FALSE,align = "right") {
#   if (is.snPackMat(M)) {M <- unpack_snPackMat(M)}
#   M <- Matrix::Matrix(M,sparse = sparse)
#   Matrix::printSpMatrix(M,digits = digits,note.dropping.colnames = note.dropping.colnames,align = align,...)
# }

# TO RECYCLE:
# #' Wrapper to shorten vectors to print with ellipsis
# #'
# #' @param l a list to print
# #' @param threshold the length of x above which the vector should be shortened with an ellipsis before being printed
# #' @param before number of elements to display before the ellipsis
# #' @param after number of elements to display after the ellipsis
# #'
# #' @return either v or a shortened version of v to display
# #' @noRd
# shorten_list.to.print <- function(l,threshold = 10,before = 3,after = 2) {
#   n <- length(l)
#   if (n >= threshold) {
#     lapply(1:before,function(i) print_list_element(l,i))
#     cat("... (",n-before," more scans)\n\n\n",sep = "")
#     lapply((n-after+1):n,function(i) print_list_element(l,i))
#   } else {
#     lapply(1:n,function(i) print_list_element(l,i))
#   }
# }
