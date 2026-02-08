
# BFun

<!-- badges: start -->

[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![License:
GNU-3](https://img.shields.io/badge/license-GNU--3-blue.svg)](https://cran.r-project.org/web/licenses/GNU-3)
[![](https://img.shields.io/github/last-commit/JSartini/BFun.svg)](https://github.com/JSartini/BFun/commits/main)
<!-- badges: end -->

The goal of BFun is to provide tools necessary to fit Bayesian
Functional (Mixed) Models where random effects are represented using
stochastic eigenfunctions. This type of joint modeling accounts for all
potential sources of uncertainty, providing valid inferences even when
signal to noise ratio is relatively small. The package is oriented
around a Markov Chain Monte Carlo approach, but provides an experimental
variational implementation as well (still undergoing testing).

For the statistical insights behind these models, refer to (Sartini et
al. 2025) and (Sartini, Zeger, and Crainiceanu 2025).

## Installation

You can install BFun from GitHub using the `remotes` package as follows:

``` r
install.packages("remotes")
library(remotes)
remotes::install_github("JSartini/BFun")
```

## Example

This is a basic example showing how to fit a Bayesian functional model
using BFun on the DTI data from the `refund` package. The model is a
simple multilevel FPCA here, as no covariates are included in the
formula.

``` r
library(BFun)
library(refund)
library(tidyr)

# Testing data
data(DTI)
fit_data = drop_na(DTI)

# Fit the model
Mod = bfmm(cca ~ 1, data = fit_data, 
           id = "ID", visit = "visit", 
           n_chains = 4, n_cores = 4, 
           n_iter = 1000, n_burnin = 1000, 
           alpha = 0.05)

# Raw samples for inference on any desired functional
Mod$Samples

# Estimates
Mod$Estimates

# Uncertainty bounds at given alpha level
Mod$Inferences

# Diagnostics (variance explained, MCMC convergence)
Mod$Diagnostics
```

# References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-sartini_bayesian_2025" class="csl-entry">

Sartini, Joseph, Scott Zeger, and Ciprian Crainiceanu. 2025. “Bayesian
Multivariate Sparse Functional PCA.” arXiv.
<https://doi.org/10.48550/arXiv.2509.03512>.

</div>

<div id="ref-Sartini08122025" class="csl-entry">

Sartini, Joseph, Xinkai Zhou, Liz Selvin, Scott Zeger, and Ciprian M.
Crainiceanu. 2025. “Fast Bayesian Functional Principal Components
Analysis.” *Journal of Computational and Graphical Statistics* 0 (ja):
1–20. <https://doi.org/10.1080/10618600.2025.2592768>.

</div>

</div>
