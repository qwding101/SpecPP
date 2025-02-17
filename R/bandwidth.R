#' Bandwidth selection procedure for kernel spectral estimator
#'
#' @param ppp A point pattern of class `"ppp"`.
#' @param inten.formula A [`formula`] syntax in character format specifying the log-liner model for the
#'  intensity function, which is passed to [`ppm`]. The default is constant intensity `inten.formula = "~1"`.
#' @param data.covariate Optional. The values of spatial covariates passed to the `data` argument in [`ppm`].
#' @param a Taper coefficient, a value within unit interval. If `a = 1`, then no data taper is used.
#' @param band.range Numeric vector. Search space for the optimal bandwidth.
#' @param correct Logical. If `TRUE` (default), conduct edge correction when computing the kernel spectral estimator.
#' @param A1,A2 Optional. Side lengths of the observation window.
#' @param endpt A positive value indicating the scale factor of the endpoint frequency.
#' @param equal Logical. Whether to use the same bandwidth for both x and y axis. The default is `TRUE`.
#' @param kern Univariate scaled kernel function, e.g., Barlett kernel (default).
#'
#' @return
#' A list.
#' @examples
#' library(spatstat)
#' lam <- function(x, y, m) {(x^2 + y) * ifelse(m == "A", 2, 1)}
#' spp <- rmpoispp(lambda = lam, win = square(5), types = c("A","B"))
#' b <- seq(0.95, 1.8, 0.1) # You may use finer grid to search optimal bandwidth
#'
#' cv <- select_band(spp, inten.formula = "~ x + y", band.range = b)
#' cv$`Optimal bandwidth` # Print the optimal bandwidth
#' plot(cv$Result[1,], cv$Result[2,], type = "b", pch = 20, las = 1,
#'      xlab = "Bandwidth", ylab = "Whittle likelihood")
#' abline(v = cv$`Optimal bandwidth`, col = "blue")
#' @importFrom foreach %dopar%
#' @export
select_band = function(ppp, inten.formula = NULL, data.covariate = NULL,
                       a = 0.025, band.range, correct = TRUE, A1 = NULL, A2 = A1,
                       endpt = 1.5, equal = TRUE, kern = bartlett_uni){
  # To please R CMD check
  r = j = k = NULL
  # Crate frequency index array (freq.idx.df) to sum in likelihood calculation
  if (is.null(A1) | is.null(A2)){
    A1 = diff(ppp$window$xrange); A2 = diff(ppp$window$yrange)
  }
  freq.list = generate_freq(A1, A2, ext.factor = NULL, return.comb = TRUE, endpt = endpt)
  freq.idx.df = cbind(omega1.ind = rep(seq_along(freq.list$omega1), length(freq.list$omega2)),
                      omega2.ind = rep(seq_along(freq.list$omega2), each = length(freq.list$omega1)),
                      val = NA)
  # Bandwidth to evaluate
  if (equal){
    band.mat = cbind(b1 = band.range, b2 = band.range)
  }else{
    band.mat = cbind(b1 = (A2/A1)*band.range, b2 = band.range)
  }

  # Matrix to save the LOOCV result
  cv = vector("list", 3)
  cv[[3]] = matrix(c(band.range, rep(NA, length(band.range))),
                   ncol = length(band.range), nrow = 2, byrow = TRUE)
  rownames(cv[[3]]) = c("Bandwidth", "Whittle likelihood")

  # Enumerate all possible combination of spectra to compute
  if (!spatstat.geom::is.multitype(ppp)){ # Univariate point pattern
    cate = 1
  }else{ # Multivariate point pattern
    cate = levels(spatstat.geom::marks(ppp))
  }
  cate.comb = pair_cate(cate)
  period.mat.list = period.mat.ext.list = vector("list", nrow(cate.comb))


  # Function to compute the Whittle likelihood without the summation part
  ind_Whittle = function(w, period.est, period.est.loo, cate.comb){

    I.val = sapply(period.est, function(x){x[w[2], w[1]]})
    F.val = sapply(period.est.loo, function(x){x[w[2], w[1]]})

    I.mat = as.matrix(stats::xtabs(val ~ ., cbind(data.frame(cate.comb), val = I.val)))
    F.mat = as.matrix(stats::xtabs(val ~ ., cbind(data.frame(cate.comb), val = F.val)))
    if (nrow(cate.comb) > 1){
      # For univariate case, just skip below code because it will be still the same
      I.mat = I.mat + Conj(t(I.mat)) - diag(diag(I.mat))
      F.mat = F.mat + Conj(t(F.mat)) - diag(diag(F.mat))
    }

    attr(F.mat, "class") = NULL
    attr(F.mat, "call") = NULL
    F.mat = fixpd_fun(F.mat)

    # return(sum(diag(I.mat %*% solve(F.mat))) + log(det(F.mat)))
    return(Re(sum(diag(I.mat %*% solve(F.mat)))) +
             log(prod(Re(eigen(F.mat, only.values = TRUE)$values))))
  }

  # Function to compute the Whittle likelihood for a specific bandwidth
  sum_Whittle = function(band.mat.row, cate.comb,
                         freq.list, freq.ext.list, freq.idx.df,
                         period.raw, period.ext, kern. = kern){

    period.smooth.mat.loo.list = vector("list", nrow(cate.comb))
    for (r in 1:nrow(cate.comb)){
      # Compute the leave-one-out kernel spectral estimator: F matrix
      val.smoothed = apply(freq.list$omega.comb, MARGIN = 1, smoother,
                           period.mat = period.ext[[r]],
                           w.k1 = freq.ext.list$omega1, w.k2 = freq.ext.list$omega2,
                           b1 = band.mat.row[1], b2 = band.mat.row[2],
                           loo = TRUE, kernel_uni = kern.)
      period.smooth.mat.loo.list[[r]] = matrix(val.smoothed,
                                               ncol=length(freq.list$omega1),
                                               nrow=length(freq.list$omega2), byrow = TRUE)
    }

    # Compute the Whittle likelihood (use na.rm = T to treat NaN as zero)
    freq.idx.df[, 3] = apply(freq.idx.df[, 1:2], MARGIN = 1, ind_Whittle,
                             period.est=period.raw, period.est.loo=period.smooth.mat.loo.list,
                             cate.comb=cate.comb)
    return(sum(freq.idx.df[, 3], na.rm = TRUE))
  }


  ### Whether to correct the edge effect in kernel smoothing by extending frequency grid
  if (correct){
    freq.ext.list = generate_freq(A1, A2, ext.factor = 2, return.comb = TRUE, endpt = endpt)
  }else{
    freq.ext.list = freq.list
  }


  ### Compute I matrix for each frequency
  # CRAN limits the number of cores available to packages to 2
  # https://stackoverflow.com/questions/50571325/r-cran-check-fail-when-using-parallel-functions
  chk = Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && chk == "TRUE") { # Use 2 cores in CRAN
    num_workers = 2L
  } else {
    num_workers = parallel::detectCores() - 1
  }
  cl = parallel::makeCluster(num_workers)
  doParallel::registerDoParallel(cl)

  DFT.list = foreach::foreach(k = seq_along(cate),
                     .export = c('periodogram', 'generate_freq', 'taper'),
                     .packages = 'spatstat') %dopar%{
                       periodogram(i = cate[k], j = cate[k], ppp = ppp,
                                  inten.formula = inten.formula, data.covariate = data.covariate,
                                  a = a, A1 = A1, A2 = A2, ext.factor = 2, return.DFT = TRUE, endpt = endpt,
                                  generate_freq. = generate_freq, taper. = taper)$DFT
                     }
  names(DFT.list) = cate

  for (r in 1:nrow(cate.comb)){
    # Calculate the periodogram for both extended and original frequency grids
    period.mat.ext.list[[r]] = DFT.list[[cate.comb[r,1]]]*Conj(DFT.list[[cate.comb[r,2]]])
    if (cate.comb[r, 1] == cate.comb[r, 2]){
      period.mat.ext.list[[r]] = Re(period.mat.ext.list[[r]])
    }
    period.mat.list[[r]] = period.mat.ext.list[[r]][dplyr::between(freq.ext.list$omega2, min(freq.list$omega2), max(freq.list$omega2)),
                                                    dplyr::between(freq.ext.list$omega1, min(freq.list$omega1), max(freq.list$omega1))]
  }

  ### For each bandwidth, conduct LOOCV
  cv[[3]][2,] = foreach::foreach(j = 1:nrow(band.mat),
                        .combine = 'c',
                        .export = c("smoother", 'fixpd_fun'),
                        .packages = 'dplyr') %dopar% {
    sum_Whittle(band.mat.row = band.mat[j,],
                    cate.comb = cate.comb,
                    freq.list = freq.list,
                    freq.ext.list = freq.ext.list,
                    freq.idx.df = freq.idx.df,
                    period.raw = period.mat.list,
                    period.ext = period.mat.ext.list,
                kern. = kern)
  }
  parallel::stopCluster(cl)

  ### parRapply() approach for LOOCV: this causes bug in Mac OS
  # cl = parallel::makeCluster(parallel::detectCores() - 1)
  # parallel::clusterExport(cl, c("Bartlett.fun","KK.fun","smooth.fun","between.fun"))
  # cv[[3]][2,] = parallel::parRapply(cl, band.mat, sum.Whittle.fun,
  #                         cate.comb = cate.comb,
  #                         freq.list = freq.list,
  #                         freq.ext.list = freq.ext.list,
  #                         freq.idx.df = freq.idx.df,
  #                         period.raw = period.mat.list,
  #                         period.ext = period.mat.ext.list)
  # parallel::stopCluster(cl)

  cv[[1]] = cv[[3]][1, which.min(cv[[3]][2,])]
  cv[[2]] = min(cv[[3]][2,])
  names(cv) = c("Optimal bandwidth", "Likelihood", "Result")

  if (which.min(cv[[3]][2,]) %in% c(1, length(band.range))){
    warning(paste0("The optimal bandwidth lies on the endpoint of the `band.range`. ",
                   "We suggest extending the range of `band.range` and executing the bandwidth search again."),
            call. = FALSE)
  }
  return(cv)
}
