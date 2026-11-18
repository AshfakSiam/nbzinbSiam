# nbzinbSiam

**nbzinbSiam** fits two regression models for overdispersed count data,
with the estimation coded by hand:

| Function | Model | Algorithm | Reference implementation |
|----|----|----|----|
| [`nbreg()`](https://ashfaksiam.github.io/nbzinbSiam/reference/nbreg.md) | Negative binomial (NB2), log link | IRLS for the coefficients, alternated with Newton-Raphson for log(theta) | [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) |
| [`zinbreg()`](https://ashfaksiam.github.io/nbzinbSiam/reference/zinbreg.md) | Zero-inflated negative binomial | EM algorithm whose M-steps are IRLS fits | `pscl::zeroinfl(..., dist = "negbin")` |

No general-purpose optimizer
([`optim()`](https://rdrr.io/r/stats/optim.html),
[`nlm()`](https://rdrr.io/r/stats/nlm.html),
[`glm()`](https://rdrr.io/r/stats/glm.html), …) is used anywhere in the
fitting code. Standard errors for the ZINB model come from the observed
information matrix, computed with an analytic Hessian (or, optionally,
numerically with
[`numDeriv::hessian()`](https://rdrr.io/pkg/numDeriv/man/hessian.html)).

This package was written for the STA 651 final project (Fall 2026, Miami
University).

**Website:** <https://ashfaksiam.github.io/nbzinbSiam/>, with the
function reference and the full project report ([PDF
version](https://github.com/AshfakSiam/nbzinbSiam/blob/main/report/report.pdf)).

## Installation

``` r

# install.packages("devtools")
devtools::install_github("AshfakSiam/nbzinbSiam")
```

## Usage

``` r

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

Both fitting functions accept an offset, either as
`offset(log(exposure))` in the formula or through the `offset` argument,
and a `control` list
([`nbreg_control()`](https://ashfaksiam.github.io/nbzinbSiam/reference/nbreg_control.md),
[`zinbreg_control()`](https://ashfaksiam.github.io/nbzinbSiam/reference/zinbreg_control.md))
with the tolerance, the maximum number of iterations and a `trace`
option.

## Available methods

Objects of class `"nbreg"` and `"zinbreg"` support
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html) (so
[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) work),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html) (response and
Pearson), [`nobs()`](https://rdrr.io/r/stats/nobs.html) and
[`predict()`](https://rdrr.io/r/stats/predict.html) (with `newdata`).
For `zinbreg`, [`predict()`](https://rdrr.io/r/stats/predict.html) has
`type = "response"` (the overall mean, (1 - pi) \* mu), `"count"` (the
count-model mean mu) and `"zero"` (the structural-zero probability pi).

## Agreement with the reference implementations

On [`pscl::bioChemists`](https://rdrr.io/pkg/pscl/man/bioChemists.html),
the package reproduces
[`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) and
[`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html); the
test suite in `tests/testthat/` checks this automatically:

| Quantity                       | `nbreg` vs `glm.nb` | `zinbreg` vs `zeroinfl` |
|--------------------------------|---------------------|-------------------------|
| Coefficients                   | \< 1e-5             | \< 1e-3                 |
| Standard errors                | \< 1e-5             | \< 1e-3                 |
| theta (NB) / log(theta) (ZINB) | \< 1e-4             | \< 1e-3                 |
| Log-likelihood                 | \< 1e-4             | \< 1e-3                 |

The benchmark and simulation scripts behind the report are in
`inst/benchmark/`.

## Notes on `R CMD check`

[`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
finishes with 0 errors, 0 warnings and 0 notes.
