#' Compute H.h.lambda.1 in the centered DFT
#'
#' @param w1,w2 A numeric value or vector of frequency values at horizontal and
#' vertical directions, respective.
#' @param a Taper coefficient, a value within unit interval. If `a = 1`, then no
#'  data taper is used.
#' @param taper Data taper function.
#' @param A1,A2 Side lengths of the observation window.
#' @param inten.fitted Intensity function of individual point pattern.

H.h.lambda.1 = Vectorize(function(w2, w1, a, taper, A1, A2, inten.fitted){
  surface = spatstat.geom::as.im(X = function(x,y,a,A1,A2,w1,w2) taper(x,a)*taper(y,a)*inten.fitted(x*A1,y*A2)*exp(-1i*(A1*x*w1+A2*y*w2)),
                                 W = spatstat.geom::owin(xrange = c(-.5,.5), yrange = c(-.5,.5)),
                                 a = a, A1 = A1, A2 = A2, w1 = w1, w2 = w2)
  return(spatstat.geom::integral.im(surface)) # This part is the computation bottleneck, which takes most of the time
}, vectorize.args = c("w2","w1") # Vectorize the argument w1 & w2 for `outer()` to use
)
