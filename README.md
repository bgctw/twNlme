
<!-- 
README.md is generated from README.Rmd. Please edit that file
rmarkdown::render(
  "README.Rmd", output_format = rmarkdown::github_document(html_preview = FALSE)) 
-->
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/twNlme)](http://cran.r-project.org/package=twNlme) [![Travis-CI Build Status](https://travis-ci.org/bgctw/twNlme.svg?branch=master)](https://travis-ci.org/bgctw/twNlme)

Overview
--------

The `twNlme` package helps with estimating uncertainty of predictions of nonlinear mixed effects models.

While the residual variance decreases with averaging over many predictions, the uncertainty due to random effects does not decrease in general.

Installation
------------

``` r
# From CRAN (in future)
#install.packages("twNlme")

# First install dependencies
install.packages("nlme")
# Install from github
# install.packages("devtools")
devtools::install_github("bgctw/twNlme")
```

Usage
-----

See the [Tree biomass vignette](https://github.com/bgctw/twNlme/blob/master/vignettes/TreeBiomassExample.md) and other [package vignettes](https://github.com/bgctw/twNlme/blob/master/vignettes/) (\*.md) for examples.
