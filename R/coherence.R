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

#' @export
print.SuppAttr = function(x, ...){
  print(matrix(as.numeric(x),
               nrow = attributes(x)$dim[1],
               ncol = attributes(x)$dim[2],
               dimnames = list(attributes(x)$dimnames[[1]],
                               attributes(x)$dimnames[[2]])),
        ...)
}
