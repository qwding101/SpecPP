#' Compute the spectrum of the intensity reweighted process
#'
#' @param w1,w2 Frequency vector (only allow frequency values evaluated in `sp.est`)
#' @param sp.est A list of kernel spectral estimate matrices of the point pattern.
#' @param i,j Optional. Index of the multivariate point process.
#' @param H.list Optional. The list from `Hmatrix()`.
#' @param ppp Point pattern.

IRspec = function(w1, w2, sp.est, i = NULL, j = NULL, H.list = NULL, ppp){

  sp.est.full = sp.est
  A1 = attr(sp.est, "A1")
  A2 = attr(sp.est, "A2")
  cate = attr(sp.est, "cate")
  cate.comb = attr(sp.est, "cate.comb")
  omega1 = attr(sp.est, "freq.list")$omega1
  omega2 = attr(sp.est, "freq.list")$omega2
  inten.fitted.scaled = attr(sp.est, "inten.fitted.scaled")

  if (is.null(i) + is.null(j) == 2){ # Return a m by m matrix, computing all pairwise local spectra
    i = j = seq_along(cate)
  }else if (is.null(i) + is.null(j) == 1){
    stop("Please specify the index for the other process, or leave both i and j unspecified.", call. = FALSE)
  }else{ # Only compute the local spectrum of the (i, j)-th pair
    cate = unique(c(i, j))
    subset.comb.idx = apply(cate.comb, 1, function(x) all(x %in% c(i, j)))
    cate.comb = cate.comb[subset.comb.idx, , drop = FALSE]
    sp.est = sp.est[subset.comb.idx]
  }

  ### Estimate the inverse Fourier transform of L2
  ## F.hat matrix given a frequency pair (w1, w2)
  w.idx = c(which(omega1 == w1), which(omega2 == w2))
  F.val = sapply(sp.est, function(x){x[w.idx[2], w.idx[1]]})
  F.mat = as.matrix(stats::xtabs(val ~ ., cbind(data.frame(cate.comb), val = F.val)))
  F.mat = F.mat + Conj(t(F.mat)) - diag(diag(F.mat))
  attr(F.mat, "class") = NULL; attr(F.mat, "call") = NULL
  ## Obtain the subcomponents in L2 if they were not supplied
  if (is.null(H.list)){
    H.list = Hmatrix(sp.est = sp.est.full, ppp = ppp)
  }
  ## Inverse Fourier transform of L2
  invFT.L2.est = (H.list$H_h2*F.mat[i, j] - (2*pi)^(-2)*H.list$H_h2l.val[i, j])/H.list$H_h2ll.val[i, j]


  ### Return the result, which is a function of w1, w2
  unit.mat = diag(1, length(inten.fitted.scaled))
  rownames(unit.mat) = colnames(unit.mat) = names(inten.fitted.scaled)
  out = (2*pi)^(-2)*unit.mat[i,j] + invFT.L2.est

  attr(out, "cate") = names(inten.fitted.scaled)
  attr(out, "eval.freq") = c(w1, w2)

  return(out)
}
