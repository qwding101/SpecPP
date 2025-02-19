#' Compute the the H matrix used in coherence.R
#'
#' @param sp.est A list of kernel spectral estimate matrices of `ppp`.
#' @param ppp Point pattern.

Hmatrix = function(sp.est, ppp){

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
