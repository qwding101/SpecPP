#' Plot all possible spectral density estimates
#'
#' @description
#' Given a multivariate (multitype) point pattern, use `ggplot2` package to
#' visualize all (marginal and joint) spectral density estimates.
#'
#' @param est.list List. The kernel spectral density estimate from [`periodogram_smooth()`].
#' @param ppp A point pattern of class `"ppp"`.
#' @param xnorm Logical. If `TRUE` (default), plot the radially-averaged spectral
#' estimates. Otherwise, plot the raw values by heatmap.
#' @param type If `type = "Re"` (default), plot the real part of the estimates.
#' If `type = "Im"`, plot the imaginary part.
#' @param shared.legend Logical. Whether to share the legend across all plots.
#' @param remove An experimental feature. Please ignore.
#'
#' @examples
#' library(spatstat)
#' lam <- function(x, y, m) {(x^2 + y) * ifelse(m == "A", 2, 1)}
#' set.seed(227823)
#' spp <- rmpoispp(lambda = lam, win = square(5), types = c("A","B"))
#' KSDE.list <- periodogram_smooth(spp, inten.formula = "~ x + y", bandwidth = 1.15)
#'
#' plot_pairs(est.list = KSDE.list, ppp = spp)
#' plot_pairs(est.list = KSDE.list, ppp = spp, type = "Im")
#' plot_pairs(est.list = KSDE.list, ppp = spp, xnorm = FALSE)
#' plot_pairs(est.list = KSDE.list, ppp = spp, xnorm = FALSE, type = "Im")
#' @importFrom ggplot2 .data
#' @export
plot_pairs = function(est.list, ppp, xnorm = TRUE, type = "Re",
                      shared.legend = TRUE, remove = NULL){

  cate = levels(spatstat.geom::marks(ppp))
  plot.layout = matrix(NA, ncol = length(cate), nrow = length(cate))
  plot.layout[lower.tri(plot.layout, diag = TRUE)] = seq_along(est.list)
  titles = rep("", length(plot.layout))
  titles[diag(plot.layout)] = cate
  freq.list = attr(est.list, "freq.list")

  est.list = switch(type,
                    Re = lapply(est.list, Re),
                    Im = lapply(est.list, Im),
                    stop("`type` should be 'Re' or 'Im'.", call. = FALSE))
  legend.range = switch(shared.legend + 1, # Basically `ifelse()` but ifelse() cannot return NULL
                        NULL,
                        c(min(sapply(est.list, min)),
                          max(sapply(est.list, max)))
  )

  period.plot.li = vector("list", length(est.list))


  if (xnorm){

    period.long = as.data.frame(freq.list$omega.comb)
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
      cowplot::draw_label(ifelse(type == "Re",
                                 paste("Radially averaged density", "(Real part)"),
                                 paste("Radially averaged density", "(Imaginary part)")),
                          x=0, y=0.5, vjust= 1.5,
                          size=12, angle=90)

  }else{
    for (r in seq_along(est.list)){
      period.plot.li[[r]] = plot_period(est.list[[r]],
                                        #freq.list = attr(est.list, "freq.list"),
                                        freq.list = freq.list,
                                        remove = remove,
                                        title = names(est.list)[r],
                                        legend.range = legend.range) +
        ggplot2::labs(title = titles[r], x = NULL, y = NULL) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    }
    if (shared.legend){
      # common.legend = cowplot::get_legend(period.plot.li[[1]])
      # Above code can work but pops out unnecessary warning message
      common.legend = cowplot::get_plot_component(period.plot.li[[1]], "guide-box-right")
      period.plot.li = lapply(period.plot.li, function(x) x + ggplot2::guides(fill = "none"))
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
