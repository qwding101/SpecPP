#' Compute \eqn{\hat{H}^{(n)}_{h\hat{\lambda}^{(i)},1}} in the centered DFT
#'
#' @description
#' This function computes \deqn{\hat{H}^{(n)}_{h\hat{\lambda}^{(i)},1}(\boldsymbol\omega) =
#' \int_{D_n}h(\boldsymbol{x}/\boldsymbol{A})\hat{\lambda}^{(i)}(\boldsymbol{x})
#' \exp(-i\boldsymbol{x}^\intercal\boldsymbol\omega)d\boldsymbol{x}}
#' for the \eqn{i}th point process, \eqn{i\in\{1,2,\ldots,m\}}.
#'
#' @param w1,w2 A numeric value or vector of frequency values at horizontal and
#' vertical directions, respectively.
#' @param a Taper coefficient, a value within unit interval. If `a = 0`, then
#'  taper is not applied, i.e., \eqn{h(\boldsymbol{x}/\boldsymbol{A}) = 1}.
#' @param taper Data taper function \eqn{h}.
#' @param A1,A2 Side lengths of the observation window.
#' @param inten.fitted Fitted intensity function of individual point pattern,
#' \eqn{\hat{\lambda}^{(i)}(\cdot)}.
#'
#' @returns A value.

H.h.lambda.1 = Vectorize(function(w2, w1, a, taper, A1, A2, inten.fitted){
  surface = spatstat.geom::as.im(X = function(x,y,a,A1,A2,w1,w2) taper(x,a)*taper(y,a)*inten.fitted(x*A1,y*A2)*exp(-1i*(A1*x*w1+A2*y*w2)),
                                 W = spatstat.geom::owin(xrange = c(-.5,.5), yrange = c(-.5,.5)),
                                 a = a, A1 = A1, A2 = A2, w1 = w1, w2 = w2)
  return(spatstat.geom::integral.im(surface)) # This part is the computation bottleneck, which takes most of the time
}, vectorize.args = c("w2","w1") # Vectorize the argument w1 & w2 for `outer()` to use
)
