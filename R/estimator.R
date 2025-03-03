#' Periodogram of a univatiate point pattern
#'
#' @description
#'
#' Given a multivariate point pattern `ppp`, compute the periodogram for each
#' marginal point pattern.
#'
#' @param i Mark index. An element in `levels(spatstat.geom::marks(ppp))`.
#' @param j Mark index. An element in `levels(spatstat.geom::marks(ppp))`.
#' @param ppp A point pattern of class `"ppp"`.
#' @param inten.formula A [`formula`] syntax in character format specifying the
#' log-liner model for the intensity function, which is passed to [`ppm()`]. The
#' default is constant intensity `inten.formula = "~1"`.
#' @param data.covariate Optional. The values of spatial covariates passed to
#' the `data` argument in [`ppm()`].
#' @param a Taper coefficient, a value within \eqn{(0,1/2)}. If `a = 0`, then
#' taper is not applied. Default is `a = 0.025`.
#' @param A1,A2 Optional. Side lengths of the observation window.
#' @param ext.factor Optional. If `NULL` (default), the frequency grid for the
#' periodogram is not extended. Please keep this `NULL` unless you know what you
#' are doing.
#' @param endpt A positive value indicating the scale factor of the endpoint frequency.
#' @param return.DFT If `TRUE`, also return the centered discrete Fourier
#' transform (DFT).
#' @param generate_freq.,taper.,H.h.lambda.1. Functions in their respective .R files.
#'
#' @returns
#' A matrix if `return.DFT = FALSE`; a list consisting of the periodogram and the
#' centered DFT if `return.DFT = TRUE`.
periodogram = function(i, j, ppp,
                       inten.formula = "~1", data.covariate = NULL,
                       a = 0.025, return.DFT = FALSE, A1 = NULL, A2 = A1,
                       ext.factor = NULL, endpt = 1.5,
                       generate_freq. = generate_freq, taper. = taper, H.h.lambda.1. = H.h.lambda.1){
  # Shift the observation window to the origin
  ppp = spatstat.geom::shift(ppp, origin = "centroid")

  if (!is.null(data.covariate)){
    data.covariate = lapply(data.covariate, spatstat.geom::shift, origin = "centroid")
  }

  if (is.null(A1) | is.null(A2)){
    A1 = diff(ppp$window$xrange); A2 = diff(ppp$window$yrange)
  }
  freq.list = generate_freq.(A1, A2, ext.factor = ext.factor, return.comb = TRUE,
                             endpt = endpt)

  center_DFT = function(i){
    pppi = ppp[spatstat.geom::marks(ppp) == i] # Select data from ith process

    if (is.null(inten.formula)){ # Estimate intensity function by kernel smoothing
      inten.fitted = as.function(spatstat.explore::density.ppp(spatstat.geom::unmark(pppi)))
    } else if (is.list(inten.formula)){ # Use population intensity functions
      stopifnot("The list should include the population intensity functions for all univariate processes."= length(inten.formula) == length(levels(spatstat.geom::marks(ppp))))
      names(inten.formula) = levels(spatstat.geom::marks(ppp))
      inten.fitted = inten.formula[[i]]
    } else if (is.character(inten.formula)){
      # Estimate intensity function by log-linear model
      inten.ppm = spatstat.model::ppm(Q = stats::as.formula(paste("unmark(pppi)", inten.formula)),
                                data = data.covariate)
      inten.fitted = as.function(stats::predict(inten.ppm))
    }

    # Compute DFT
    const = ((2*pi)^(-2/2))*(spatstat.geom::area(pppi)^(-1/2))/(1+a*(5/(4*pi^2) - 4/3))
    V1 = outer(freq.list$omega1, pppi$x, function(w, x, a, A) taper.(x/A, a)*exp(-1i*x*w),
               a=a, A=A1)
    V2 = outer(freq.list$omega2, pppi$y, function(w, x, a, A) taper.(x/A, a)*exp(-1i*x*w),
               a=a, A=A2)
    J_h.woconst = V2 %*% t(V1)
    mat.left = outer(freq.list$omega2, freq.list$omega1[freq.list$omega1 < 0], H.h.lambda.1.,
                     a=a, taper=taper., A1=A1, A2=A2, inten.fitted=inten.fitted)
    mat.center = outer(freq.list$omega2, 0, H.h.lambda.1.,
                       a=a, taper=taper., A1=A1, A2=A2, inten.fitted=inten.fitted)
    mat.right = Conj(matrix(rev(as.vector(mat.left)), ncol = ncol(mat.left), nrow = nrow(mat.left)))
    C_h.woconst = spatstat.geom::area(pppi)*cbind(mat.left, mat.center, mat.right)
    # Center the DFT
    DFT = const * (J_h.woconst - C_h.woconst)
    attr(DFT, "inten.fitted") = inten.fitted
    return(DFT)
  }


  if (!spatstat.geom::is.multitype(ppp)){
    spatstat.geom::marks(ppp) = 1
  }

  if (i == j){

    DFT = center_DFT(i) # Centered DFT
    colnames(DFT) = freq.list$omega1
    rownames(DFT) = freq.list$omega2
    period.mat = Mod(DFT)^2 # Periodogram for ith process
    attr(period.mat, "freq.list") = freq.list

    if (return.DFT){
      # Return both DFT and periodogram
      # DFT of both ith and jth processes can be used to compute the cross-spectrum
      period.DFT = list(Periodogram = period.mat, DFT = DFT)
      attr(period.DFT, "freq.list") = freq.list
      return(period.DFT)
    }else{
      # Only return the periodogram
      return(period.mat)
    }
  }
  else{
    # Calculate the cross-spectrum directly
    period.mat = center_DFT(i)*Conj(center_DFT(j))
    colnames(period.mat) = freq.list$omega1
    rownames(period.mat) = freq.list$omega2
    attr(period.mat, "freq.list") = freq.list
    return(period.mat)
  }
}



#' Kernel spectral estimator of a multivatiate (multitype) point pattern
#'
#' @description
#'
#' Computes the kernel spectral estimator for a multivatiate (multitype) point pattern.
#'
#' @param i Mark index. An element in `levels(spatstat.geom::marks(ppp))`.
#' @param j Mark index. An element in `levels(spatstat.geom::marks(ppp))`.
#' @param ppp A point pattern of class `"ppp"`.
#' @param inten.formula A [`formula`] syntax in character format specifying the
#' log-liner model for the intensity function, which is passed to [`ppm()`]. The
#' default is constant intensity `inten.formula = "~1"`.
#' @param data.covariate Optional. The values of spatial covariates passed to
#' the `data` argument in [`ppm`].
#' @param bandwidth A positive value indicating the bandwidth of kernel,
#' determined by [`select_band()`].
#' @param correct Logical. If `TRUE` (default), conduct edge correction when
#' computing the kernel spectral estimator.
#' @param a Taper coefficient, a value within \eqn{(0,1/2)}. If `a = 0`, then
#' taper is not applied. Default is `a = 0.025`.
#' @param A1,A2 Optional. Side lengths of the observation window.
#' @param endpt A positive value indicating the scale factor of the endpoint frequency.
#' @param equal Logical. If `TRUE`, then use the same bandwidth for both x and y direction.
#' @param kern Univariate scaled kernel function. The default is Barrlett kernel.
#'
#' @returns
#' A list of matrices, or a single matrix if `i` and `j` are specified.
#'
#' @details
#' The minimal required arguments are `ppp`, `inten.formula`, and `bandwidth`.
#' If you use any spatial covariate other than the Cartesian coordinates in `inten.formula`, then
#' `data.covariate` is also needed. All the other arguments can be left by default setting.
#' [`periodogram_smooth()`] computes all the pairwise (marginal and cross-) kernel
#' spectral estimators automatically when the mark indices `i` and `j` are
#' unspecified. If `i` and `j` are specified, then it only computes the result
#' for that mark combination.
#'
#' The bandwidth can be determined by [`select_band()`].
#'
#' @examples
#' library(spatstat)
#' lam <- function(x, y, m) {(x^2 + y) * ifelse(m == "A", 2, 1)}
#' spp <- rmpoispp(lambda = lam, win = square(5), types = c("A","B"))
#'
#' # Compute kernel spectral estimator with intensity fitted by log-linear model
#' # with Cartesian coordinates
#' ksde = periodogram_smooth(spp, inten.formula = "~ x + y", bandwidth = 1.2)
#' lapply(ksde, round, 2)
#' @importFrom foreach %dopar%
#' @export
periodogram_smooth = function(ppp, i = NULL, j = i,
                              inten.formula = "~1", data.covariate = NULL,
                              bandwidth, correct = TRUE, a = 0.025,
                              A1 = NULL, A2 = A1, endpt = 1.5, equal = TRUE,
                              kern = bartlett_uni){
  ### To please R CMD check
  r = j = k = NULL
  while(FALSE){
    spatstat::spatstat.family()
  }

  ### Construct frequency grids
  if (is.null(A1) | is.null(A2)){
    A1 = diff(ppp$window$xrange); A2 = diff(ppp$window$yrange)
  }
  freq.list = generate_freq(A1, A2, ext.factor = NULL, return.comb = TRUE, endpt = endpt) # Original grid
  omega.comb = freq.list$omega.comb
  if (correct){
    freq.ext.list = generate_freq(A1, A2, ext.factor = 2, return.comb = TRUE, endpt = endpt) # Grid for edge-correction
  }else{
    freq.ext.list = freq.list
  }

  if (equal){
    b = rep(bandwidth, 2)
  }else{
    b = c(A2/A1, 1)*bandwidth
  }
  names(b) = c("b1", "b2")

  ### List all pairwise combinations
  if (!(is.null(i) | is.null(j))){ # If both i and j are specified, compute the spectrum for [i,j]
    cate.comb = cbind(i, j)
    cate = 1
  }else{ # Otherwise, enumerate all combinations of spectra to compute
    if (!spatstat.geom::is.multitype(ppp)){ # Univariate point pattern
      cate = 1
    }else{ # Multivariate point pattern
      cate = levels(spatstat.geom::marks(ppp))
    }
    cate.comb = generate_cate(cate)
  }

  ### Spectral estimation
  if (length(cate) == 1){ # For univariate process, no need to use parallel computing
    # Calculate the periodogram
    period.mat.ext = periodogram(i = cate.comb[1,1], j = cate.comb[1,2], ppp = ppp,
                                inten.formula = inten.formula, data.covariate = data.covariate,
                                a = a, A1 = A1, A2 = A2, ext.factor = 2)
    # Calculate the kernel smoothed estimator
    val = apply(omega.comb, MARGIN = 1, smoother,
                period.mat = period.mat.ext,
                w.k1 = freq.ext.list$omega1, w.k2 = freq.ext.list$omega2,
                b1 = b[1], b2 = b[2])
    # Return the kernel smoothed estimator
    period.smooth.mat.list = list()
    period.smooth.mat.list[[1]] = matrix(val,
                                         ncol = length(freq.list$omega1),
                                         nrow = length(freq.list$omega2), byrow = TRUE)
    colnames(period.smooth.mat.list[[1]]) = freq.list$omega1
    rownames(period.smooth.mat.list[[1]]) = freq.list$omega2

  }else{ # For multivariate process, use parallel computing
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
                       .export = c('periodogram', 'generate_freq', 'taper', 'H.h.lambda.1'),
                       .packages = 'spatstat') %dopar%{
                         periodogram(i = cate[k], j = cate[k], ppp = ppp,
                                    inten.formula = inten.formula, data.covariate = data.covariate,
                                    a = a, A1 = A1, A2 = A2, ext.factor = 2, return.DFT = TRUE, endpt = endpt,
                                    generate_freq. = generate_freq, taper. = taper, H.h.lambda.1. = H.h.lambda.1)$DFT
                       }
    names(DFT.list) = cate

    period.smooth.mat.list = foreach::foreach(r = 1:nrow(cate.comb),
                                     .export = c("smoother")) %dopar% {
                                       # Calculate the periodogram
                                       period.mat.ext = DFT.list[[cate.comb[r,1]]]*Conj(DFT.list[[cate.comb[r,2]]])
                                       if (cate.comb[r,1] == cate.comb[r,2]){
                                         period.mat.ext = Re(period.mat.ext)
                                       }
                                       # Calculate the kernel smoothed estimator
                                       val = apply(omega.comb, MARGIN = 1, smoother,
                                                   period.mat = period.mat.ext,
                                                   w.k1 = freq.ext.list$omega1, w.k2 = freq.ext.list$omega2,
                                                   b1 = b[1], b2 = b[2], kernel_uni = kern)
                                       # Return the kernel smoothed estimator
                                       period.smooth.mat = matrix(val,
                                                                  ncol = length(freq.list$omega1),
                                                                  nrow = length(freq.list$omega2), byrow = TRUE)
                                       colnames(period.smooth.mat) = freq.list$omega1
                                       rownames(period.smooth.mat) = freq.list$omega2
                                       period.smooth.mat
                                     }
    parallel::stopCluster(cl)

    inten.fitted.list = DFT.list
  }



  names(period.smooth.mat.list) = paste0(cate.comb[,1], ", ", cate.comb[,2])
  if (nrow(cate.comb) == 1){
    out = period.smooth.mat.list[[1]]
  }else{
    out = period.smooth.mat.list
  }

  # Store some metadata
  inten.fitted.list = lapply(inten.fitted.list, attr, which = "inten.fitted")
  names(inten.fitted.list) = cate
  inten.fitted.scaled = lapply(inten.fitted.list, scale2unitarea, length.x = A1, length.y = A2)
  attr(out, "inten.fitted.scaled") = inten.fitted.scaled
  attr(out, "a") = a
  attr(out, "A1") = A1
  attr(out, "A2") = A2
  attr(out, "freq.list") = freq.list
  attr(out, "cate") = cate
  attr(out, "cate.comb") = cate.comb
  attr(out, "bandwidth") = b

  return(out)
}
