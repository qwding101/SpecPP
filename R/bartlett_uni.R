#' Scaled Battlett kernel
#'
#' @description
#' This function calculates the scaled kernel \eqn{K_b(v) = b^{-1}K(v/b)} where
#' \deqn{K(v) = \max\{1 - |v|\}.}
#'
#' @param v Input value.
#' @param b Bandwidth.
bartlett_uni = function(v, b){
  return(ifelse(abs(v/b) > 1, 0, (1-abs(v/b))/b))
}
