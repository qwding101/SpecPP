#' Generate all pairwise combinations of
#'
#' @param cate Vector containing the categorical mark of a point process, which
#' can be obtained from `levels(spatstat.geom::marks(ppp))`, where `ppp` is the
#' point pattern.

generate_cate = function(cate){
  pos1 = pos2 = c()
  for (i in seq_along(cate)){
    pos1 = append(pos1, rep(cate[i], length(cate) - i + 1))
    pos2 = append(pos2, cate[i:length(cate)])
  }
  return(cbind(pos1, pos2))
}
