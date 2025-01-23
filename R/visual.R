#' Plot individual spectral estimate
#'
#' @description
#' Use `ggplot2` package to plot the heatmap of a single periodogram or kernel spectral estimate.
#'
#' @param period Matrix. Output from `periodogram()` or [periodogram_smooth()].
#' @param freq.list List. Output from `generate_freq()` with `return.comb = TRUE`.
#' @param remove An experimental feature. Please ignore.
#' @param title Character. The title of the plot.
#' @param palette Color schemes for visualization. Type `?scale_fill_distiller()`
#'  to check available options.
#' @param legend.range Input `c(lower, upper)` to manually set the lower and upper limit of the spectral density to visualize.
#' @seealso [plot_pairs()]
#' @examples
#' library(spatstat)
#' lam <- function(x, y, m) {(x^2 + y) * ifelse(m == "A", 2, 1)}
#' spp <- rmpoispp(lambda = lam, win = square(5), types = c("A","B"))
#' period.list <- periodogram_smooth(spp, inten.formula = "~ x + y", bandwidth = 1.2)
#' names(period.list)
#'
#' # Only plot the real part of the periodogram
#' period.list <- lapply(period.list, Re)
#' freq.list <- generate_freq(A1 = diff(spp$window$xrange), A2 = diff(spp$window$yrange),
#'                            return.comb = TRUE)
#' plot_period(period.list[[1]], freq.list = freq.list,
#'             title = "Kernel spectral estimate: type A")
#' plot_period(period.list[[2]], freq.list = freq.list,
#'             title = "Kernel spectral estimate: cross-spectrum")
#' plot_period(period.list[[3]], freq.list = freq.list,
#'             title = "Kernel spectral estimate: type B")
#' @importFrom ggplot2 .data
#' @export
plot_period = function(period, freq.list = NULL, remove=NULL,
                       title = NULL, palette = "Spectral", legend.range = NULL){

  if (is.matrix(period)){
    if (is.null(freq.list)){
      omega.comb = attr(period, "freq.list")$omega.comb
    }else{
      omega.comb = freq.list$omega.comb
    }
    period.long = as.data.frame(omega.comb)
    period.long$dist = sqrt(period.long$omega1^2 + period.long$omega2^2)
    period.long$Density = as.vector(t(period))
  }
  else if (is.data.frame(period)){
    period.long = period
  }
  else{
    stop("The periodogram should be a matrix or long-format data.frame.", call. = FALSE)
  }

  if (!is.null(remove)){
    period.long$Density[dplyr::between(period.long$omega1, -remove[1], remove[1]) &
                        dplyr::between(period.long$omega2, -remove[2], remove[2])] = NA
  }

  ggplot2::ggplot(period.long, ggplot2::aes(.data$omega1, .data$omega2)) +
    ggplot2::geom_raster(ggplot2::aes(fill = .data$Density)) +
    ggplot2::scale_fill_distiller(palette = palette, name = "", direction = 1,
                         limit = legend.range, na.value = "transparent") +
    ggplot2::labs(title = title, x = expression(omega[1]), y = expression(omega[2])) +
    ggplot2::theme_classic()
}



#' Plot all possible spectral estimates
#'
#' @description
#' Given a multivariate (multitype) point pattern, use `ggplot2` package to
#' visualize all (marginal and joint) kernel spectral estimates.
#'
#' @param est.list List. Output from [periodogram_smooth()].
#' @param ppp A point pattern of class `"ppp"`.
#' @param xnorm Logical. If `TRUE` (default), plot the radially-averaged spectral
#' estimates. Otherwise, plot the raw values by heatmap.
#' @param freq.list List. Output from `generate_freq()` with `return.comb = TRUE`.
#' @param shared.legend Logical. Whether to share the legend across all plots.
#' @param remove An experimental feature. Please ignore.
#'
#' @examples
#' library(spatstat)
#' lam <- function(x, y, m) {(x^2 + y) * ifelse(m == "A", 2, 1)}
#' spp <- rmpoispp(lambda = lam, win = square(5), types = c("A","B"))
#' period.list <- periodogram_smooth(spp, inten.formula = "~ x + y", bandwidth = 1.2)
#'
#' # Only plot the real part of the periodogram
#' period.list <- lapply(period.list, Re)
#' freq.list <- generate_freq(A1 = diff(spp$window$xrange), A2 = diff(spp$window$yrange),
#'                            return.comb = TRUE)
#' plot_pairs(est.list = period.list, ppp = spp, freq.list = freq.list)
#' @importFrom ggplot2 .data
#' @export
plot_pairs = function(est.list, ppp, xnorm = TRUE, freq.list = NULL,
                      shared.legend = TRUE, remove = NULL){
  cate = levels(spatstat.geom::marks(ppp))
  plot.layout = matrix(NA, ncol = length(cate), nrow = length(cate))
  plot.layout[lower.tri(plot.layout, diag = TRUE)] = seq_along(est.list)

  titles = rep("", length(plot.layout))
  titles[diag(plot.layout)] = cate

  legend.range = switch(shared.legend + 1, # Basically `ifelse()` but ifelse() cannot return NULL
                        NULL,
                        c(min(sapply(est.list, min)),
                          max(sapply(est.list, max)))
  )

  period.plot.li = vector("list", length(est.list))


  if (xnorm){
    if (is.null(freq.list)){
      omega.comb = attr(est.list, "freq.list")$omega.comb
    }else{
      omega.comb = freq.list$omega.comb
    }
    period.long = as.data.frame(omega.comb)
    period.long$dist = sqrt(period.long$omega1^2 + period.long$omega2^2)

    for (r in seq_along(est.list)){

      period.long$Density = as.vector(t(est.list[[r]]))
      period.plot.li[[r]] = dplyr::group_by(period.long, .data$dist) %>%
        dplyr::summarise(avg.density = mean(.data$Density)) %>%
        ggplot2::ggplot(ggplot2::aes(.data$dist, .data$avg.density)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
        ggplot2::theme_classic() +
        ggplot2::labs(title = titles[r], x = NULL, y = NULL) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                       panel.background = ggplot2::element_rect(fill = "transparent"),
                       plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank()) +
        ggplot2::coord_cartesian(ylim = legend.range)
    }

    result = cowplot::plot_grid(plotlist = period.plot.li[as.vector(plot.layout)],
                       nrow = length(cate), ncol = length(cate), byrow = FALSE,
                       scale = 0.85)
    result = result +
      cowplot::draw_label(expression(paste("||", bold("\u03c9"), "||")),
                          x=0.5, y=0, vjust=-0.5, size=12, angle=0) +
      cowplot::draw_label("Radially averaged density", x=0, y=0.5, vjust= 1.5,
                          size=12, angle=90)

  }else{
    for (r in seq_along(est.list)){
      period.plot.li[[r]] = plot_period(est.list[[r]],
                                        freq.list = attr(est.list, "freq.list"),
                                        remove = remove,
                                        title = names(est.list)[r],
                                        legend.range = legend.range) +
        ggplot2::labs(title = titles[r], x = NULL, y = NULL) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    }
    if (shared.legend){
      period.plot.li = lapply(period.plot.li, function(x) x + ggplot2::guides(fill = "none"))
      common.legend = cowplot::get_legend(plot_period(est.list[[1]],
                                                      legend.range = legend.range))
      result = cowplot::plot_grid(plotlist = period.plot.li[as.vector(plot.layout)],
                                  nrow = length(cate), ncol = length(cate),
                                  byrow = F, scale = 0.85)
      result = cowplot::plot_grid(result, common.legend, ncol = 2, rel_widths = c(1, .1))

    }else{
      result = cowplot::plot_grid(plotlist = period.plot.li[as.vector(plot.layout)],
                                  nrow = length(cate), ncol = length(cate), byrow = FALSE,
                                  scale = 0.85)
    }
    # https://stackoverflow.com/questions/49577461/adding-x-and-y-laxis-label-to-ggplot-grid-build-with-cowplot
    result = result +
      cowplot::draw_label(expression(omega[1]), x=0.5, y=0, vjust=-0.5, angle=0) +
      cowplot::draw_label(expression(omega[2]), x=0, y=0.5, vjust= 1.5, angle=90)
  }

  return(result)
}



#' Plot coherence and partial coherence
#'
#' @description
#' Use `ggplot2` package to visualize the coherence and partial coherence.
#'
#' @param sp.est List. The kernel spectral estimate from [periodogram_smooth()].
#' @param coh.mat Coherence matrix from [coherence()] with `type = "normal"`.
#' @param partial.coh.mat Partial coherence matrix from [coherence()] with `type = "partial"`.
#' @param xnorm Logical. If `TRUE` (default), plot the radial averaged values. If `FALSE`, plot
#' the raw coherence and partial coherence values via heatmap.
#' @param ylim A numeric vector `c(lower, upper)` to specify the range to draw for the
#' radially-averaged plot. Not required if `xnorm = FALSE`.
#'
#' @seealso [coherence()]
#' @examples
#' library(spatstat)
#' lam <- function(x, y, m) {(x^2 + y) * ifelse(m == "A", 2, 1)}
#' spp <- rmpoispp(lambda = lam, win = square(5), types = c("A","B"))
#'
#' # Compute kernel spectral estimator with fitted intensity by log-linear model: HERE
#' # with Cartesian coordinates
#' spectra <- periodogram_smooth(spp, inten.formula = "~ x + y", bandwidth = 1.2)
#' coh.partial <- coherence(spectra, spp)
#' coh <- coherence(spectra, spp, type = "normal")
#' plot_coher(spectra, coh, coh.partial)
#' @importFrom ggplot2 .data
#' @export
plot_coher = function(sp.est, coh.mat, partial.coh.mat, xnorm = TRUE, ylim = NULL){

  cate = attr(sp.est, "cate")
  cate.comb = data.frame(attr(sp.est, "cate.comb"))
  coh.allw = data.frame(cbind(#attr(sp.est, "freq.list")$omega.comb[, 1:2],
    attr(coh.mat, "CohTable"),
    attr(partial.coh.mat, "CohTable")[, -(1:2)])
  )
  coh.allw$dist = sqrt(coh.allw$omega1^2 + coh.allw$omega2^2)
  colnames(coh.allw) = c("omega1", "omega2",
                         paste0("C", 1:(nrow(cate.comb)-length(cate))),
                         paste0("PC", 1:(nrow(cate.comb)-length(cate))), "dist")

  plot.order = diag(rep(999, length(cate)))
  plot.order[upper.tri(plot.order)] = paste0("C", 1:(nrow(cate.comb)-length(cate)))
  plot.order[lower.tri(plot.order)] = t(plot.order)[lower.tri(plot.order)]
  plot.order[lower.tri(plot.order)] = paste0("P", plot.order[lower.tri(plot.order)])
  diag(plot.order) = cate
  plot.order = as.vector(plot.order)


  plot.list = vector("list", length = length(plot.order))

  if (xnorm){

    for (k in seq_along(plot.order)){
      if (plot.order[k] %in% cate){
        plot.list[[k]] = ggplot2::ggplot() +
          ggplot2::annotate("text", x = 1, y = 1, size = 5, label = plot.order[k]) +
          ggplot2::theme_void()
      }else{
        df = coh.allw[, c("omega1","omega2", "dist", plot.order[k])]
        colnames(df)[4] = "Value"
        plot.list[[k]] = dplyr::group_by(df, .data$dist) %>%
          dplyr::summarise(avg.value = mean(.data$Value)) %>%
          ggplot2::ggplot(ggplot2::aes(.data$dist, .data$avg.value)) +
          ggplot2::geom_line() +
          ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
          ggplot2::theme_classic() +
          ggplot2::labs(x = NULL, y = NULL) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                panel.background = ggplot2::element_rect(fill='transparent'),
                plot.background = ggplot2::element_rect(fill='transparent', color=NA),
                panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank()) +
          ggplot2::coord_cartesian(ylim = ylim)
      }
    }
    common.title = cowplot::ggdraw() +
      cowplot::draw_label("Radially averaged (partial) coherence",
                           fontface = 'bold', x = 0, hjust = 0)
    result = cowplot::plot_grid(plotlist = plot.list,
                                nrow = length(cate), ncol = length(cate),
                                byrow = FALSE, scale = 0.8)
    result = result +
      cowplot::draw_label(expression(paste("||", bold("\u03c9"), "||")),
                          x=0.5, y=0, vjust=-0.5, size=12, angle=0) +
      cowplot::draw_label("Coherence", x=1, y=0.5, vjust=1.5, angle=-90) +
      cowplot::draw_label("Partial coherence", x=0, y=0.5, vjust= 1.5, angle=90)
    result = cowplot::plot_grid(common.title, result, ncol = 1, rel_heights = c(0.05, 1))

  }else{

    for (k in seq_along(plot.order)){
      if (plot.order[k] %in% cate){
        plot.list[[k]] = ggplot2::ggplot() +
          ggplot2::annotate("text", x = 1, y = 1, size = 5, label = plot.order[k]) +
          ggplot2::theme_void()
      }else{
        df = coh.allw[, c("omega1","omega2", plot.order[k])]
        colnames(df)[3] = "Density"
        plot.list[[k]] = plot_period(df, palette = "Spectral") +
          ggplot2::labs(x = NULL, y = NULL)
      }
    }
    result = cowplot::plot_grid(plotlist = plot.list,
                                nrow = length(cate), ncol = length(cate),
                                byrow = F, scale = 0.8)
    # https://stackoverflow.com/questions/49577461/adding-x-and-y-laxis-label-to-ggplot-grid-build-with-cowplot
    result = result +
      cowplot::draw_label("Coherence", x=1, y=0.5, vjust=1.5, angle=-90) +
      cowplot::draw_label("Partial coherence", x=0, y=0.5, vjust= 1.5, angle=90)
  }

  return(result)
}
