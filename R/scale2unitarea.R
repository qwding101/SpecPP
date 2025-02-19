#' Recale the rectangular image to unit square
#' @param im Image object.
#' @param length.x,length.y Side lengths.
#' @param return.fun Logical. If `TRUE`, convert the scaled image to function.
#' Otherwise, return the image itself.

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
