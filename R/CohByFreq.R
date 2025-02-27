#' (Partial) coherence for a given frequency
#'
#' @param w1,w2 Frequency vector (only allow frequency values evaluated in `sp.est`)
#' @param sp.est A list of kernel spectral estimate matrices of the point pattern.
#' @param type If `type = 'partial'`, compute partial coherence. If `type = 'normal'`,
#'  compute coherence.
#' @param i,j Optional. Index of the multivariate point process.
#' @param sp.IR Optional. Spectrum estimate of the intensity reweighted process, which is
#' calculated by [`IRspec()`]. If this argument is specified, then `w1`, `w2`, `H.list`,
#' `ppp` are not required.
#' @param H.list Optional. A list from [`Hmatrix()`]. This argument is only required
#' if `sp.IR` is unspecified.
#' @param ppp Optional. Point pattern. This argument is only required if `H.list`
#'  is unspecified.

CohByFreq = function(w1, w2, sp.est, type = "partial", i = NULL, j = NULL,
                     sp.IR = NULL, H.list = NULL, ppp = NULL){

  # To please R CMD check
  r = NULL

  if (is.null(sp.IR)){
    sp.IR = IRspec(w1, w2, sp.est, H.list = H.list, ppp = ppp)
    sp.IR = fixpd_fun(sp.IR)
    nom = switch(type,
                 partial = solve(sp.IR),
                 normal = sp.IR,
                 stop("`type` should be 'partial' or 'normal'.", call. = FALSE))
    cate = attr(sp.est, "cate")
  }else if (is.matrix(sp.IR)){
    sp.IR = fixpd_fun(sp.IR)
    nom = switch(type,
                 partial = solve(sp.IR),
                 normal = sp.IR,
                 stop("`type` should be 'partial' or 'normal'.", call. = FALSE))
    cate = attr(sp.IR, "cate")
  }else{
    stop("`sp.IR` should be a matrix.", call. = FALSE)
  }

  # D = diag(1/diag(nom))
  # coh.mat = D %*% (nom)^2 %*% D
  D = Mod(diag(1/diag(nom)))
  coh.mat = D %*% Mod(nom)^2 %*% D

  colnames(coh.mat) = cate
  row.names(coh.mat) = cate

  if (is.null(i) + is.null(j) == 2){ # Return a m by m matrix
    return(coh.mat)
  }else if (is.null(i) + is.null(j) == 1){
    stop("Please specify the index for the other process, or leave both i and j unspecified.",
         call. = FALSE)
  }else{ # Only compute the coherence of the (i, j)-th pair
    return(coh.mat[i, j])
  }
}
