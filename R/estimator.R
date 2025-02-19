# Periodogram of a univatiate point pattern
# Return a matrix if `return.DFT = FALSE`; a list consisting of the periodogram and DFT if `return.DFT = TRUE`.
periodogram = function(i, j, ppp,
                       inten.formula = "~1", data.covariate = NULL,
                       a = 0.025, return.DFT = FALSE, A1 = NULL, A2 = A1,
                       ext.factor = NULL, endpt = 1.5,
                       generate_freq. = generate_freq, taper. = taper){
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

  exp_term = function(w, x, a, A){
    return(taper.(x/A, a)*exp(-1i*x*w))
  }

  H.h.lambda.1 = Vectorize(function(w2, w1, a, A1, A2, inten.fitted){
    surface = spatstat.geom::as.im(X = function(x,y,a,A1,A2,w1,w2) taper.(x,a)*taper.(y,a)*inten.fitted(x*A1,y*A2)*exp(-1i*(A1*x*w1+A2*y*w2)),
                    W = spatstat.geom::owin(xrange = c(-.5,.5), yrange = c(-.5,.5)),
                    a = a, A1 = A1, A2 = A2, w1 = w1, w2 = w2)
    return(spatstat.geom::integral.im(surface)) # This part is the computation bottleneck, which takes most of the time
  }, vectorize.args = c("w2","w1") # Vectorize the argument w1 & w2 for `outer()` to use
  )

  Jh_fun = function(i){
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
    V1 = outer(freq.list$omega1, pppi$x, exp_term, a=a, A=A1)
    V2 = outer(freq.list$omega2, pppi$y, exp_term, a=a, A=A2)
    J_h.woconst = V2 %*% t(V1)
    mat.left = outer(freq.list$omega2, freq.list$omega1[freq.list$omega1 < 0], H.h.lambda.1,
                     a=a, A1=A1, A2=A2, inten.fitted=inten.fitted)
    mat.center = outer(freq.list$omega2, 0, H.h.lambda.1,
                       a=a, A1=A1, A2=A2, inten.fitted=inten.fitted)
    mat.right = Conj(matrix(rev(as.vector(mat.left)), ncol = ncol(mat.left), nrow = nrow(mat.left)))
    C_h.woconst = spatstat.geom::area(pppi)*cbind(mat.left, mat.center, mat.right)

    DFT = const * (J_h.woconst - C_h.woconst)
    attr(DFT, "inten.fitted") = inten.fitted
    return(DFT)
  }


  if (!spatstat.geom::is.multitype(ppp)){
    spatstat.geom::marks(ppp) = 1
  }

  if (i == j){

    DFT = Jh_fun(i) # Centered DFT
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
    period.mat = Jh_fun(i)*Conj(Jh_fun(j))
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
#' @param i Mark index. An element in `marks(ppp)`.
#' @param j Mark index. An element in `marks(ppp)`.
#' @param ppp A point pattern of class `"ppp"`.
#' @param inten.formula A [`formula`] syntax in character format specifying the log-liner model for the
#'  intensity function, which is passed to [`ppm`]. The default is constant intensity `inten.formula = "~1"`.
#' @param data.covariate Optional. The values of spatial covariates passed to the `data` argument in [`ppm`].
#' @param bandwidth A positive value indicating the bandwidth of kernel, determined by [select_band()].
#' @param correct Logical. If `TRUE` (default), conduct edge correction when computing the kernel spectral estimator.
#' @param a Taper coefficient, a value within unit interval. If `a = 1`, then no data taper is used.
#'  Default is `a = 0.025`.
#' @param A1,A2 Optional. Side lengths of the observation window.
#' @param endpt A positive value indicating the scale factor of the endpoint frequency.
#' @param equal Logical. If `TRUE`, then use the same bandwidth for both x and y direction.
#' @param kern Univariate scaled kernel function, e.g., Barrlett kernel (default).
#'
#' @return
#' A list of matrices, or a single matrix if `i` and `j` are specified.
#'
#' @details
#' The minimal required arguments are `ppp`, `inten.formula`, and `bandwidth`.
#' If you use any spatial covariate (other than the Cartesian coordinates) in `inten.formula`, then
#' `data.covariate` is also needed. All the other arguments can be left by default setting.
#' `periodogram_smooth()` computes all the pairwise (marginal and cross-) kernel spectral estimators
#' automatically when the mark indices `i` and `j` are unspecified. If `i` and `j` are specified,
#' then it only computes the result for that mark combination.
#'
#' The bandwidth can be determined by [select_band()].
#'
#' @examples
#' library(spatstat)
#' lam <- function(x, y, m) {(x^2 + y) * ifelse(m == "A", 2, 1)}
#' spp <- rmpoispp(lambda = lam, win = square(5), types = c("A","B"))
#'
#' # Compute kernel spectral estimator with intensity fitted by log-linear model
#' # with Cartesian coordinates
#' periodogram_smooth(spp, inten.formula = "~ x + y", bandwidth = 1.2)
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

  ### List all combinations
  if (!(is.null(i) | is.null(j))){ # If both i and j are specified, compute the spectrum for [i,j]
    cate.comb = cbind(i, j)
    cate = 1
  }else{ # Otherwise, enumerate all combinations of spectra to compute
    if (!spatstat.geom::is.multitype(ppp)){ # Univariate point pattern
      cate = 1
    }else{ # Multivariate point pattern
      cate = levels(spatstat.geom::marks(ppp))
    }
    cate.comb = pair_cate(cate)
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
                       .export = c('periodogram', 'generate_freq', 'taper'),
                       .packages = 'spatstat') %dopar%{
                         periodogram(i = cate[k], j = cate[k], ppp = ppp,
                                    inten.formula = inten.formula, data.covariate = data.covariate,
                                    a = a, A1 = A1, A2 = A2, ext.factor = 2, return.DFT = TRUE, endpt = endpt,
                                    generate_freq. = generate_freq, taper. = taper)$DFT
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


#' Maximum coherence and partial coherence matrix
#'
#' @description
#'
#' Computes the coherence or partial coherence matrix for a multivatiate (multitype) point pattern.
#'
#' @param sp.est List. The kernel spectral estimate from [periodogram_smooth()].
#' @param ppp A point pattern of class `"ppp"`.
#' @param type If `"partial"` (default), calculate the partial coherence. If `"normal"`, calculate the coherence.
#' @param all Logical. If `TRUE`, extract the maximum (partial) coherence across all frequencies, which is not recommended.
#'
#' @return
#' A matrix with each entry storing the maximum (partial) coherence of the two individual point processes.
#'
#' @examples
#' library(spatstat)
#' lam <- function(x, y, m) {(x^2 + y) * ifelse(m == "A", 2, 1)}
#' spp <- rmpoispp(lambda = lam, win = square(5), types = c("A","B"))
#'
#' # Compute kernel spectral estimator with fitted intensity by log-linear model:
#' # with Cartesian coordinates
#' spectra <- periodogram_smooth(spp, inten.formula = "~ x + y", bandwidth = 1.2)
#' coh.partial <- coherence(spectra, spp) # Compute the maximum partial coherence
#' attr(coh.partial, "CohTable") # Print the partial coherence values for all frequencies
#'
#' @importFrom foreach %dopar%
#' @export
coherence = function(sp.est, ppp, type = "partial", all = FALSE){

  # To please R CMD check
  r = NULL

  H.list = Hmatrix(sp.est = sp.est, ppp = ppp)
  cate = attr(sp.est, "cate")
  omega.comb = attr(sp.est, "freq.list")$omega.comb[, 1:2]

  if (all){ # Compute (partial) coherence for all frequencies.
    cat(paste0("Use all", nrow(omega.comb), "frequencies to pick the maximum",
               ifelse(type == "partial", " partial", ""),
               " coherence."))
  }else{ # Compute (partial) coherence only for frequencies which share no overlapping smoothing neighbors
    omega1 = attr(sp.est, "freq.list")$omega1
    omega2 = attr(sp.est, "freq.list")$omega2
    b = attr(sp.est, "bandwidth")
    origin.idx = c((length(omega1)+1)/2, (length(omega2)+1)/2)
    omega.diff = c(omega1[2]-omega1[1], omega2[2]-omega2[1])
    idx.diff = ceiling(2*b/omega.diff)
    omega1 = omega1[c(rev(seq(from = origin.idx[1], to = 1, by = -idx.diff[1])[-1]),
                      seq(from = origin.idx[1], to = length(omega1), by = idx.diff[1]))]
    omega2 = omega2[c(rev(seq(from = origin.idx[2], to = 1, by = -idx.diff[2])[-1]),
                      seq(from = origin.idx[2], to = length(omega2), by = idx.diff[2]))]
    omega.comb = cbind(omega1 = rep(omega1, length(omega2)),
                       omega2 = rep(omega2, each = length(omega1)))

    cat(paste0("Number of frequencies to pick the maximum",
               ifelse(type == "partial", " partial", ""),
               " coherence:"),
        nrow(omega.comb),
        "(this value should not be too small).")
  }

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
  out = foreach::foreach(r = 1:nrow(omega.comb), .combine = "rbind",
                         .export = c("IRspec", "CohByFreq", "fixpd_fun"),
                         .inorder = T) %dopar% {
                  sp.IR = IRspec(w1 = omega.comb[r,1], w2 = omega.comb[r,2],
                                  sp.est = sp.est, H.list = H.list)
                  coh.mat = CohByFreq(sp.IR = sp.IR, H.list = H.list, type = type)
                  return(coh.mat[upper.tri(coh.mat, diag = F)])
                }
  parallel::stopCluster(cl)


  coh.w = cbind(omega.comb, out)
  max.coh = apply(coh.w[, -c(1:2), drop = FALSE], MARGIN = 2, max)
  max.coh.mat = diag(1, length(cate))
  max.coh.mat[upper.tri(max.coh.mat, diag = F)] = max.coh
  max.coh.mat = max.coh.mat + t(max.coh.mat) - diag(diag(max.coh.mat))
  row.names(max.coh.mat) = colnames(max.coh.mat) = cate

  attr(max.coh.mat, 'CohTable') = coh.w
  class(max.coh.mat) = "SuppAttr" # Define a new class for customized print method which does not print attribute

  return(max.coh.mat)
}
