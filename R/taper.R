#' Univariate data taper
#'
#' @description
#' This function implements
#' \deqn{
#' h(x) =
#' \begin{cases}
#' (x+\frac{1}{2})/a - \frac{1}{2\pi}\sin\left(2\pi(x+\frac{1}{2})/a\right), &
#' -\frac{1}{2} \leq x \leq -\frac{1}{2}+a \\
#' 1, & -\frac{1}{2}+a \leq x \leq \frac{1}{2}-a \\
#' h(-x), & \frac{1}{2}-a \leq x \leq \frac{1}{2}
#' \end{cases}
#' ,}
#' where \eqn{a\in (0,1/2)}. If \eqn{a = 0}, set \eqn{h(x) = 1}.
#'
#' @param x Cartesian coordinate of events. This argument is vectorized.
#' @param a A value within \eqn{(0,1/2)}. If `a = 0`, then taper is not applied.

taper = function(x, a){
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
