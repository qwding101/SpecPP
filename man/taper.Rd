% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taper.R
\name{taper}
\alias{taper}
\title{Univariate data taper}
\usage{
taper(x, a)
}
\arguments{
\item{x}{Cartesian coordinate of events. This argument is vectorized.}

\item{a}{A value within \eqn{(0,1/2)}. If \code{a = 0}, then taper is not applied.}
}
\description{
This function implements
\deqn{
h(x) =
\begin{cases}
(x+\frac{1}{2})/a - \frac{1}{2\pi}\sin\left(2\pi(x+\frac{1}{2})/a\right), &
-\frac{1}{2} \leq x \leq -\frac{1}{2}+a \\
1, & -\frac{1}{2}+a \leq x \leq \frac{1}{2}-a \\
h(-x), & \frac{1}{2}-a \leq x \leq \frac{1}{2}
\end{cases}
,}
where \eqn{a\in (0,1/2)}. If \eqn{a = 0}, set \eqn{h(x) = 1}.
}
