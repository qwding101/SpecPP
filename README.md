
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SpecPP

<!-- badges: start -->
<!-- badges: end -->

This R package performs spectral analysis of multivariate spatial point
patterns observed in a rectangular region ([Ding et al.,
2025](https://arxiv.org/abs/2502.09948)). By specifying a parametric
model for the intensity, it extracts the frequency‐domain
characteristics of your data by smoothing the periodogram. Key features
of this tool include spectral density estimation, cross-validated
bandwidth selection, spectrum visualization, and coherence analysis.

## Installation

You can install the development version of SpecPP from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("qwding101/SpecPP")
```

## FAQ

**Q:** What is a *(spatial) point pattern*?

**A:** Please refer to [link
1](https://geographicdata.science/book/notebooks/08_point_pattern_analysis.html)
or [link
2](https://documentation.sas.com/doc/en/pgmsascdc/v_063/statug/statug_spp_overview02.htm)
for an introduction to point patterns. Point patterns are very common in
daily life, e.g.,
[earthquakes](https://earthquake.usgs.gov/earthquakes/map/?extent=-89.1006,-517.5&extent=89.1006,66.09375&range=month&magnitude=4.5&baseLayer=terrain&settings=true),
[crimes](https://www.crimemapping.com/map/agency/98), [disease
outbreaks](https://www.healthmap.org/en/), [geotagged
photos](https://www.kaggle.com/datasets/ifeanyichukwunwobodo/tokyo-geotagged-flickr-images).

**Q:** How to use this R package?

**A:** [This
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
