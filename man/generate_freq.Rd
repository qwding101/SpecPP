% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_freq.R
\name{generate_freq}
\alias{generate_freq}
\title{Frequency grid generator}
\usage{
generate_freq(A1, A2, ext.factor = NULL, return.comb = FALSE, endpt = 1.5)
}
\arguments{
\item{A1, A2}{Side lengths of the observational window.}

\item{ext.factor}{A positive value indicating the extension factor of frequency.
If \code{NULL}, no extension is conducted.}

\item{return.comb}{Logical. If \code{TURE}, also return the \code{data.frame}.}

\item{endpt}{A positive value indicating the scale factor of the endpoint frequency.}
}
\value{
A list of two frequency vectors. If \code{return.comb = TRUE}, also return a
\code{data.frame} containing all combinations of frequencies.
}
\description{
Frequency grid generator
}
\details{
This function generates the frequency grid which is used in \code{periodogram()} and \code{\link[=periodogram_smooth]{periodogram_smooth()}}.
In \code{\link[=periodogram_smooth]{periodogram_smooth()}}, to correct the edge effect, the \code{ext.factor} argument in \code{\link[=generate_freq]{generate_freq()}}
is set to \code{1.5} by default. Suppose the original frequency ranges from \eqn{-1.5\pi} to \eqn{1.5\pi}.
Then the extended frequency ranges from \eqn{-1.5\pi\times\text{ext.factor}} to
\eqn{-1.5\pi\times\text{ext.factor}}.
}
\examples{
generate_freq(A1 = 8, A2 = 5, return.comb = TRUE)
}
