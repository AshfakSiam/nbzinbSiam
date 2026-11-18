# nbzinbSiam

<!-- badges: start -->
[![R-CMD-check](https://github.com/AshfakSiam/nbzinbSiam/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AshfakSiam/nbzinbSiam/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**nbzinbSiam** fits two regression models for overdispersed count data,
with the estimation coded by hand:

| Function | Model | Algorithm | Reference implementation |
|---|---|---|---|
| `nbreg()` | Negative binomial (NB2), log link | IRLS for the coefficients, alternated with Newton-Raphson for log(theta) | `MASS::glm.nb()` |
| `zinbreg()` | Zero-inflated negative binomial | EM algorithm whose M-steps are IRLS fits | `pscl::zeroinfl(..., dist = "negbin")` |

No general-purpose optimizer (`optim()`, `nlm()`, `glm()`, ...) is used
anywhere in the fitting code. Standard errors for the ZINB model come from
the observed information matrix, computed with an analytic Hessian (or,
optionally, numerically with `numDeriv::hessian()`).

This package was written for the STA 651 final project (Fall 2026,
Miami University).

**Website:** <https://ashfaksiam.github.io/nbzinbSiam/>, with the function
reference and the full project report
([PDF version](https://github.com/AshfakSiam/nbzinbSiam/blob/main/report/report.pdf)).

## Installation

```r
# install.packages("devtools")
devtools::install_github("AshfakSiam/nbzinbSiam")
```

## Usage

```r
library(nbzinbSiam)
data("bioChemists", package = "pscl")

# Negative binomial regression
fit_nb <- nbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
summary(fit_nb)

# Zero-inflated negative binomial regression:
# count-model covariates before |, zero-model covariates after it
fit_zinb <- zinbreg(art ~ fem + mar + kid5 + phd + ment | ment,
                    data = bioChemists)
summary(fit_zinb)

# Model comparison
AIC(fit_nb, fit_zinb)

# Predictions: overall mean, count-model mean, structural-zero probability
nd <- data.frame(fem = "Women", mar = "Married", kid5 = 0, phd = 3,
                 ment = c(0, 10, 20))
predict(fit_zinb, newdata = nd, type = "response")
predict(fit_zinb, newdata = nd, type = "count")
predict(fit_zinb, newdata = nd, type = "zero")
```

Both fitting functions accept an offset, either as `offset(log(exposure))`
in the formula or through the `offset` argument, and a `control` list
(`nbreg_control()`, `zinbreg_control()`) with the tolerance, the maximum
number of iterations and a `trace` option.

## Available methods

Objects of class `"nbreg"` and `"zinbreg"` support `print()`, `summary()`,
`coef()`, `vcov()`, `logLik()` (so `AIC()` and `BIC()` work), `fitted()`,
`residuals()` (response and Pearson), `nobs()` and `predict()` (with
`newdata`). For `zinbreg`, `predict()` has `type = "response"`
(the overall mean, (1 - pi) * mu), `"count"` (the count-model mean mu) and
`"zero"` (the structural-zero probability pi).

## Agreement with the reference implementations

On `pscl::bioChemists`, the package reproduces `MASS::glm.nb()` and
`pscl::zeroinfl()`; the test suite in `tests/testthat/` checks this
automatically:

| Quantity | `nbreg` vs `glm.nb` | `zinbreg` vs `zeroinfl` |
|---|---|---|
| Coefficients | < 1e-5 | < 1e-3 |
| Standard errors | < 1e-5 | < 1e-3 |
| theta (NB) / log(theta) (ZINB) | < 1e-4 | < 1e-3 |
| Log-likelihood | < 1e-4 | < 1e-3 |

The benchmark and simulation scripts behind the report are in
`inst/benchmark/`.

## Notes on `R CMD check`

`devtools::check()` finishes with 0 errors, 0 warnings and 0 notes.
