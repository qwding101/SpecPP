#' Scaled Battlett kernel (univariate)
#'
#' @param v Input value.
#' @param b Bandwidth.
bartlett_uni = function(v, b){
  return(ifelse(abs(v/b) > 1, 0, (1-abs(v/b))/b))
}
