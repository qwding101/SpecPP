#' Smooth an element from a matrix by averaging its neighbors
#'
#' @description
#' This function applies kernel smoothing on a matrix, \eqn{\hat{I}(\cdot)}
#' under our context, as follows:
#' \deqn{\hat{F}_{\boldsymbol{b}}(\boldsymbol\omega) = \frac{\sum_{\boldsymbol{k}
#' \in\mathbb{Z}^2}K_{\boldsymbol{b}}(\boldsymbol\omega - \boldsymbol{x_{k,\Omega}})
#' \hat{I}(\boldsymbol{x_{k,\Omega}})}{\sum_{\boldsymbol{k}\in\mathbb{Z}^2}K_{
#' \boldsymbol{b}}(\boldsymbol\omega - \boldsymbol{x_{k,\Omega}})}}
#'
#' @param w Input frequency \eqn{\boldsymbol\omega = (\omega_1,\omega_2)^\intercal}.
#' @param period.mat Matrix. The naive spectral estimate, i.e., periodogram \eqn{
#' \hat{I}(\cdot)}, could be complex-valued.
#' @param w.k1,w.k2 Vectors containing the whole frequencies at horizontal and
#' vertical directions.
#' @param b1,b2 Numeric. Bandwidth vector \eqn{\boldsymbol b = (b_1,b_2)^\intercal}.
#' @param loo Logical. If `TRUE`, conduct leave-one-out kernel smoothing (the
#' center `w` will be excluded when averaging). Otherwise, keep the center `w`
#' for averaging.
#' @param kernel_uni Univariate kernel function \eqn{K}.
#'
#' @return A matrix, \eqn{\hat{F}_{\boldsymbol{b}}(\boldsymbol\omega)}.

smoother = function(w, period.mat, w.k1, w.k2, b1=1, b2=b1,
                    loo = FALSE, kernel_uni){

  kernel_prod = function(x1, x2, omega1, omega2, b1, b2=b1){
    return(kernel_uni(omega1 - x1, b1)*kernel_uni(omega2 - x2, b2))
  }

  # dx = w.k1[2] - w.k1[1]
  dy = w.k2[2] - w.k2[1]
  # w1.pt2smooth.theory = 2*floor(b1/dx) + 1 # Number of points to smooth for x-direction in theory
  # w2.pt2smooth.theory = 2*floor(b2/dy) + 1 # Number of points to smooth for y-direction in theory

  # Narrow the range of discretized values to iterate,
  # which is related to the kernel bandwidth b1 and b2
  idx.w1 = w.k1 >= w[1] - b1 & w.k1 <= w[1] + b1
  idx.w2 = w.k2 >= w[2] - b2 & w.k2 <= w[2] + b2
  w1.range = w.k1[idx.w1]
  w2.range = w.k2[idx.w2]
  w1.pt2smooth.actual = length(w1.range) # Number of points to smooth for x-direction in reality
  w2.pt2smooth.actual = length(w2.range) # Number of points to smooth for y-direction in reality

  # Detect small bandwidth problem: there is no any neighbor to smooth for at least one direction
  if (w1.pt2smooth.actual == 1 || w2.pt2smooth.actual == 1){
    if (loo){
      stop(paste0("Either the minimal bandwidth is too small or the frequency grid is too coarse. ",
                  "Given current frequency grid, the bandwidth should be larger than ", round(dy,4), "."),
           call. = FALSE)
    }else{
      stop(paste0("Either the bandwidth is too small or the frequency grid is too coarse. ",
                  "Given current frequency grid, the bandwidth should be larger than ", round(dy,4), "."),
           call. = FALSE)
    }
  }

  # Given a (w[c], w[r]) pair, calculate the kernel product matrix
  # for the range of discretized values x1.range and x2.range.
  kernel.product = t(outer(X = w1.range, Y = w2.range, kernel_prod,
                           omega1 = w[1], omega2 = w[2], b1 = b1, b2 = b2))

  if (loo){ # For LOOCV, remove the kernel density for central point
    kernel.product[which(w2.range == w[2]), which(w1.range == w[1])] = 0
  }

  # Volume = height of endpoint * subrectangle area, where area = dx*dy,
  # which is cancelled out in below ratio
  numerator = sum((kernel.product*period.mat[idx.w2, idx.w1]))
  denominator = sum(kernel.product)

  return(ifelse(denominator == 0, 0, numerator/denominator))
}
