#' Frequency grid generator
#'
#' @description
#' This function generates the frequency grid which is used in `periodogram()`
#' and [periodogram_smooth()].
#'
#' @param A1,A2 Side lengths of the observational window.
#' @param ext.factor A positive value indicating the extension factor of frequency.
#' If `NULL`, no extension is conducted.
#' @param return.comb Logical. If `TURE`, also return the `data.frame`.
#' @param endpt A positive value indicating the scale factor of the endpoint frequency.
#'
#' @returns A list of two frequency vectors. If `return.comb = TRUE`, also return a
#' `data.frame` containing all combinations of frequencies.
#' @details
#' The frequency grid \eqn{(\omega_1,\omega_2)} is constructed
#' as follows. For \eqn{i=1,2}, \eqn{\omega_i = 1.5\pi k/A_i}, where \eqn{k\in\{
#' -A_i,-A_i+1,\ldots,A_i\}} with side length \eqn{A_i}. The original frequency
#' ranges from \eqn{-1.5\pi} to \eqn{1.5\pi}.
#'
#' The extended frequency ranges from \eqn{-1.5\pi\times\text{ext.factor}} to
#' \eqn{-1.5\pi\times\text{ext.factor}}. This grid is used in [periodogram_smooth()]
#' to correct the edge effect by setting `correct = TRUE`, where `ext.factor`
#' argument in [generate_freq()] is set to `2`.
#' @examples
#' generate_freq(A1 = 8, A2 = 5, return.comb = TRUE)
#' @export
generate_freq = function(A1, A2, ext.factor = NULL, return.comb = FALSE, endpt = 1.5){

  freq = list()
  # Generate frequency grid
  stopifnot("The value of side lengths should be larger than zero." = A1*A2 > 0)
  A1 = round(A1); A2 = round(A2)
  k1 = -A1:A1; k2 = -A2:A2
  omega1 = endpt*pi*k1/A1; omega2 = endpt*pi*k2/A2

  # Extend frequency gird to handle boundary effect when conducting kernel smoothing
  # E.g., original frequency: -1.5π ~ 1.5π
  #    => extended frequency: -1.5π*ext.factor ~ 1.5π*ext.factor
  if (!is.null(ext.factor)){
    dx = omega1[2] - omega1[1]
    dy = omega2[2] - omega2[1]
    ext.x = seq(max(omega1), max(omega1)*ext.factor, by = dx)[-1]
    ext.y = seq(max(omega2), max(omega2)*ext.factor, by = dy)[-1]
    omega1 = c(sort(-ext.x), omega1, ext.x)
    omega2 = c(sort(-ext.y), omega2, ext.y)
  }
  freq[[1]] = omega1
  freq[[2]] = omega2
  names(freq) = c("omega1", "omega2")

  if (return.comb){
    # Below code is basically `as.matrix(expand.grid(list(omega1, omega2, NA)))` but much faster
    freq[[3]] = cbind(omega1 = rep(omega1, length(omega2)),
                      omega2 = rep(omega2, each = length(omega1)),
                      Density = NA)
    names(freq)[3] = "omega.comb"
  }

  return(freq)
}

