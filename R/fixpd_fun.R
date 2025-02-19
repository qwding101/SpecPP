#' Fix the positive definite issue for a matrix
#' @param mat A complex-valued matrix.

fixpd_fun = function(mat){
  min.eigen = min(eigen(mat, only.values = TRUE)$values)
  if (min.eigen <= 0){
    mat = mat + (1e-3 - min.eigen)*diag(1, nrow(mat))
  }
  return(mat)
}
