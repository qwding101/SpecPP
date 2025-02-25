#' Plot coherence and partial coherence
#'
#' @description
#' Use `ggplot2` package to visualize the coherence and partial coherence.
#'
#' @param sp.est List. The kernel spectral density estimate from [`periodogram_smooth()`].
#' @param coh.mat Coherence matrix from [`coherence()`] with `type = "normal"`.
#' @param partial.coh.mat Partial coherence matrix from [`coherence()`] with `type = "partial"`.
#' @param xnorm Logical. If `TRUE` (default), plot the radial averaged values.
#' If `FALSE`, plot the raw coherence and partial coherence values via heatmap.
#' @param ylim A numeric vector `c(lower, upper)` to specify the range to draw
#' for the radially-averaged plot. Not required if `xnorm = FALSE`.
#'
#' @seealso [`coherence()`]
#' @examples
#' library(spatstat)
#' lam <- function(x, y, m) {(x^2 + y) * ifelse(m == "A", 2, 1)}
#' set.seed(227823)
#' spp <- rmpoispp(lambda = lam, win = square(5), types = c("A","B"))
#'
#' KSDE.list <- periodogram_smooth(spp, inten.formula = "~ x + y", bandwidth = 1.15)
#' coh.partial <- coherence(KSDE.list, spp)
#' coh <- coherence(KSDE.list, spp, type = "normal")
#' plot_coher(KSDE.list, coh, coh.partial)
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
