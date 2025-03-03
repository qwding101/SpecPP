#' Plot individual spectral density estimate
#'
#' @description
#' Use `ggplot2` package to plot the heatmap of a single periodogram or kernel
#' spectral density estimate.
#'
#' @param period Matrix. Output from [`periodogram()`] or [`periodogram_smooth()`].
#' @param freq.list List. Output from [`generate_freq()`] with `return.comb = TRUE`.
#' @param type If `type = "Re"` (default), plot the real part of the estimate.
#' If `type = "Im"`, plot the imaginary part. This argument is useful for cross-
#' spectrum estimate, which is complex-valued.
#' @param title Character. The title of the plot.
#' @param palette Color schemes for visualization. Type `?ggplot2::scale_fill_distiller()`
#'  to check available options.
#' @param legend.range Specify `c(lower, upper)` to manually set the lower and upper
#' bound of the spectral density to visualize.
#' @seealso [`plot_pairs()`]
#' @examples
#' library(spatstat)
#' lam <- function(x, y, m) {(x^2 + y) * ifelse(m == "A", 2, 1)}
#' set.seed(227823)
#' spp <- rmpoispp(lambda = lam, win = square(5), types = c("A","B"))
#' KSDE.list <- periodogram_smooth(spp, inten.formula = "~ x + y", bandwidth = 1.15)
#' names(KSDE.list)
#'
#' freq.list <- attr(KSDE.list, "freq.list")
#' plot_period(KSDE.list[[1]], freq.list = freq.list,
#'             title = "Kernel spectral estimate: type A")
#' plot_period(KSDE.list[[2]], freq.list = freq.list,
#'             title = "Kernel spectral estimate: cross-spectrum (Real part)")
#' plot_period(KSDE.list[[2]], freq.list = freq.list, type = "Im",
#'             title = "Kernel spectral estimate: cross-spectrum (Imaginary part)")
#' plot_period(KSDE.list[[3]], freq.list = freq.list,
#'             title = "Kernel spectral estimate: type B")
#' @importFrom ggplot2 .data
#' @export
plot_period = function(period, freq.list = NULL, type = "Re", # remove=NULL,
                       title = NULL, palette = "Spectral", legend.range = NULL){

  if (is.matrix(period)){
    if (is.null(freq.list)){
      omega.comb = attr(period, "freq.list")$omega.comb
    }else{
      omega.comb = freq.list$omega.comb
    }

    period = switch(type,
                    Re = Re(period),
                    Im = Im(period),
                    stop("`type` should be 'Re' or 'Im'.", call. = FALSE))
    period.long = as.data.frame(omega.comb)
    # period.long$dist = sqrt(period.long$omega1^2 + period.long$omega2^2)
    period.long$Density = as.vector(t(period))
  }
  else if (is.data.frame(period)){ # Used for plotting coherence
    period.long = period
  }
  else{
    stop("`period` should be a matrix or long-format data.frame.", call. = FALSE)
  }

  # if (!is.null(remove)){
  #   period.long$Density[dplyr::between(period.long$omega1, -remove[1], remove[1]) &
  #                         dplyr::between(period.long$omega2, -remove[2], remove[2])] = NA
  # }

  ggplot2::ggplot(period.long, ggplot2::aes(.data$omega1, .data$omega2)) +
    ggplot2::geom_raster(ggplot2::aes(fill = .data$Density)) +
    ggplot2::scale_fill_distiller(palette = palette, name = "", direction = 1,
                                  limit = legend.range, na.value = "transparent") +
    ggplot2::labs(title = title, x = expression(omega[1]), y = expression(omega[2])) +
    ggplot2::theme_classic()
}
