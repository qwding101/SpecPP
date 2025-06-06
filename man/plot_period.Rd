% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_period.R
\name{plot_period}
\alias{plot_period}
\title{Plot individual spectral density estimate}
\usage{
plot_period(
  period,
  freq.list = NULL,
  type = "Re",
  title = NULL,
  palette = "Spectral",
  legend.range = NULL
)
}
\arguments{
\item{period}{Matrix. Output from \code{\link[=periodogram]{periodogram()}} or \code{\link[=periodogram_smooth]{periodogram_smooth()}}.}

\item{freq.list}{List. Output from \code{\link[=generate_freq]{generate_freq()}} with \code{return.comb = TRUE}.}

\item{type}{If \code{type = "Re"} (default), plot the real part of the estimate.
If \code{type = "Im"}, plot the imaginary part. This argument is useful for cross-
spectrum estimate, which is complex-valued.}

\item{title}{Character. The title of the plot.}

\item{palette}{Color schemes for visualization. Type \code{?ggplot2::scale_fill_distiller()}
to check available options.}

\item{legend.range}{Specify \code{c(lower, upper)} to manually set the lower and upper
bound of the spectral density to visualize.}
}
\description{
Use \code{ggplot2} package to plot the heatmap of a single periodogram or kernel
spectral density estimate.
}
\examples{
library(spatstat)
lam <- function(x, y, m) {(x^2 + y) * ifelse(m == "A", 2, 1)}
set.seed(227823)
spp <- rmpoispp(lambda = lam, win = square(5), types = c("A","B"))
KSDE.list <- periodogram_smooth(spp, inten.formula = "~ x + y", bandwidth = 1.15)
names(KSDE.list)

freq.list <- attr(KSDE.list, "freq.list")
plot_period(KSDE.list[[1]], freq.list = freq.list,
            title = "Kernel spectral estimate: type A")
plot_period(KSDE.list[[2]], freq.list = freq.list,
            title = "Kernel spectral estimate: cross-spectrum (Real part)")
plot_period(KSDE.list[[2]], freq.list = freq.list, type = "Im",
            title = "Kernel spectral estimate: cross-spectrum (Imaginary part)")
plot_period(KSDE.list[[3]], freq.list = freq.list,
            title = "Kernel spectral estimate: type B")
}
\seealso{
\code{\link[=plot_pairs]{plot_pairs()}}
}
