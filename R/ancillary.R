fixpd_fun = function(mat){
  min.eigen = min(eigen(mat, only.values = TRUE)$values)
  if (min.eigen <= 0){
    mat = mat + (1e-3 - min.eigen)*diag(1, nrow(mat))
  }
  return(mat)
}

pair_cate = function(cate){
  pos1 = pos2 = c()
  for (i in seq_along(cate)){
    pos1 = append(pos1, rep(cate[i], length(cate)-i+1))
    pos2 = append(pos2, cate[i:length(cate)])
  }
  return(cbind(pos1, pos2))
}

#' Frequency grid generator
#'
#' @param A1,A2 Side lengths of the observational window.
#' @param ext.factor A positive value indicating the extension factor of frequency.
#' If `NULL`, no extension is conducted.
#' @param return.comb Logical. If `TURE`, also return the `data.frame`.
#' @param endpt A positive value indicating the scale factor of the endpoint frequency.
#'
#' @return A list of two frequency vectors. If `return.comb = TRUE`, also return a
#' `data.frame` containing all combinations of frequencies.
#' @details
#' This function generates the frequency grid which is used in `periodogram()` and [periodogram_smooth()].
#' In [periodogram_smooth()], to correct the edge effect, the `ext.factor` argument in [generate_freq()]
#' is set to `2` by default. Suppose the original frequency ranges from `-1.5*π` to `1.5*π`. Then the
#' extended frequency ranges from `-1.5*π*ext.factor` to `1.5*π*ext.factor`.
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


#' Univariate data taper
#'
#' @param x Cartesian coordinate of events. This argument is vectorized.
#' @param a A value in unit interval.
taper = function(x, a){
  # x is vectorized, while a is a value
  if (a == 0){
    result = rep(1, length(x))
  }else{
    result = rep(NA, length(x))
    idx = x >= -0.5 & x <= a - 0.5
    result[idx] = ((x[idx]+0.5)/a) - ((1/(2*pi))*sin(2*pi*(x[idx]+0.5)/a))
    result[x > a-0.5 & x < 0.5-a] = 1
    idx = x >= 0.5-a & x <= 0.5
    result[idx] = ((0.5-x[idx])/a) - ((1/(2*pi))*sin(2*pi*(0.5-x[idx])/a))
  }
  return(result)
}


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

scale2unitarea = function(im, length.x, length.y, return.fun = TRUE){
  if (!spatstat.geom::is.im(im)){
    im = spatstat.geom::as.im.funxy(im)
  }
  im.new = im
  im.new$xrange = im.new$yrange = c(-.5, .5)
  im.new$xstep = im$xstep/length.x
  im.new$ystep = im$ystep/length.y
  im.new$xcol = im$xcol/length.x
  im.new$ycol = im$ycol/length.y

  if (return.fun){
    return(spatstat.geom::as.function.im(im.new))
  }else{
    return(im.new)
  }
}

Hmatrix = function(sp.est, ppp){
  # sp.est: a list of kernel spectral estimate matrices calculated from `ppp`
  # ppp: point pattern

  inten.fitted.scaled = attr(sp.est, "inten.fitted.scaled")
  a = attr(sp.est, "a")
  cate = attr(sp.est, "cate")
  cate.comb = attr(sp.est, "cate.comb")

  ppp = spatstat.geom::shift(ppp, origin = "centroid")
  A1 = diff(ppp$window$xrange); A2 = diff(ppp$window$yrange)
  ppp.li = spatstat.geom::split.ppp(ppp)

  H_h2 = (a*(5/(4*(pi^2)) - 4/3) + 1)^2
  H_h2l.val = matrix(0, nrow = length(cate), ncol = length(cate), dimnames = list(cate, cate))
  H_h2ll.val = matrix(NA, nrow = length(cate), ncol = length(cate), dimnames = list(cate, cate))

  for (k in seq_along(cate)){
    H_h2l.val[k,k] = sum((taper(spatstat.geom::coords.ppp(ppp.li[[cate[k]]])$x/A1, a)^2)*
                           (taper(spatstat.geom::coords.ppp(ppp.li[[cate[k]]])$y/A2, a)^2))/
      spatstat.geom::area(ppp)
  }

  for (r in 1:nrow(cate.comb)){
    Xi = cate.comb[r, 1]
    Xj = cate.comb[r, 2]
    Xi.x = spatstat.geom::coords.ppp(ppp.li[[Xi]])$x/A1
    Xi.y = spatstat.geom::coords.ppp(ppp.li[[Xi]])$y/A2
    Xj.x = spatstat.geom::coords.ppp(ppp.li[[Xj]])$x/A1
    Xj.y = spatstat.geom::coords.ppp(ppp.li[[Xj]])$y/A2

    H_h2ll.val[Xi, Xj] = 0.5*(sum((taper(Xi.x, a)^2)*(taper(Xi.y, a)^2)*inten.fitted.scaled[[Xj]](x = Xi.x, y = Xi.y)) +
                                sum((taper(Xj.x, a)^2)*(taper(Xj.y, a)^2)*inten.fitted.scaled[[Xi]](x = Xj.x, y = Xj.y)))/spatstat.geom::area(ppp)
    H_h2ll.val[Xj, Xi] = H_h2ll.val[Xi, Xj]
  }

  return(list(H_h2 = H_h2, H_h2l.val = H_h2l.val, H_h2ll.val = H_h2ll.val))
}

IRspec = function(w1, w2, sp.est, i = NULL, j = NULL, H.list = NULL, ppp){
  # w1, w2: frequency vector (only allow frequency values evaluated in `sp.est`)
  # sp.est: a list of kernel spectral estimate matrices of the original point process
  # i, j: pair index (optional)
  # H.list: the result from `Hmatrix()` (optional)
  # ppp: point pattern (this argument is only required if you don't specify `H.list`)

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

CohByFreq = function(w1, w2, sp.est, type = "partial", i = NULL, j = NULL,
                     sp.IR = NULL, H.list = NULL, ppp = NULL){
  # w1, w2: frequency vector (only allow frequency values evaluated in `sp.est`)
  # sp.est: a list of kernel spectral estimate matrices of the point process
  # partial: if TRUE, compute partial coherence; otherwise, compute coherence
  # i, j: pair index (optional)
  # sp.IR: spectrum estimate of the IR process (if this is provided, then w1, w2, H.list, ppp are not required)
  # H.list: the result from `Hmatrix()` (this argument is only required if you don't specify `sp.IR`)
  # ppp: point pattern (this argument is only required if you don't specify `H.list`)

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

#' @export
print.SuppAttr = function(x, ...){
  print(matrix(as.numeric(x),
               nrow = attributes(x)$dim[1],
               ncol = attributes(x)$dim[2],
               dimnames = list(attributes(x)$dimnames[[1]],
                               attributes(x)$dimnames[[2]])),
        ...)
}
