
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SpecPP

<!-- badges: start -->
<!-- badges: end -->

An R package for spectral analysis of multivariate spatial point
patterns in a rectangular region ([Ding et al.,
2025](https://arxiv.org/abs/2502.09948)). Given a multivariate point
pattern and a parametric model of the intensity function, this package
facilitates the characterization of frequency-domain features from your
data through kernel smoothing. It includes functions for spectral
density estimation, bandwidth selection, visualization, and coherence
analysis.

## Installation

You can install the development version of SpecPP from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("qwding101/SpecPP")
```

## Tutorial

[This
tutorial](https://qwding101.github.io/SpecPP/articles/lansing.html)
demonstrates how to use SpecPP to analyze multivariate point pattern
data.

## References

- Ding, Q. W., Yang, J., & Shin, J. (2025). Pseudo-spectra of
  multivariate inhomogeneous spatial point processes. arXiv preprint
  [arXiv:2502.09948](https://arxiv.org/abs/2502.09948).
- Yang, J., & Guan, Y. (2024). Fourier analysis of spatial point
  processes. arXiv preprint
  [arXiv:2401.06403](https://arxiv.org/abs/2401.06403).

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
