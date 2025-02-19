#' Univariate data taper
#'
#' @param x Cartesian coordinate of events. This argument is vectorized.
#' @param a A value in unit interval.

taper = function(x, a){
  # x is vectorized, while a is a value
  if (a == 0){
    result = rep(1, length(x))
  }else{
    result = rep(NA, length(x))
    idx = x >= -0.5 & x <= a - 0.5
    result[idx] = ((x[idx]+0.5)/a) - ((1/(2*pi))*sin(2*pi*(x[idx]+0.5)/a))
    result[x > a-0.5 & x < 0.5-a] = 1
    idx = x >= 0.5-a & x <= 0.5
    result[idx] = ((0.5-x[idx])/a) - ((1/(2*pi))*sin(2*pi*(0.5-x[idx])/a))
  }
  return(result)
}
