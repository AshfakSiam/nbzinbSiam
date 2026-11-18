# Methods for zinbreg objects

Standard S3 methods for models fitted by
[`zinbreg()`](https://ashfaksiam.github.io/nbzinbSiam/reference/zinbreg.md).

## Usage

``` r
# S3 method for class 'zinbreg'
print(x, digits = max(3, getOption("digits") - 3), ...)

# S3 method for class 'zinbreg'
summary(object, ...)

# S3 method for class 'summary.zinbreg'
print(x, digits = max(3, getOption("digits") - 3), ...)

# S3 method for class 'zinbreg'
coef(object, model = c("full", "count", "zero"), ...)

# S3 method for class 'zinbreg'
vcov(object, model = c("full", "count", "zero"), ...)

# S3 method for class 'zinbreg'
logLik(object, ...)

# S3 method for class 'zinbreg'
nobs(object, ...)

# S3 method for class 'zinbreg'
fitted(object, ...)

# S3 method for class 'zinbreg'
residuals(object, type = c("response", "pearson"), ...)

# S3 method for class 'zinbreg'
predict(object, newdata, type = c("response", "count", "zero"), ...)
```

## Arguments

- x, object:

  An object of class `"zinbreg"` (or `"summary.zinbreg"` for the `print`
  method of a summary).

- digits:

  Number of significant digits to print.

- ...:

  Further arguments (currently ignored).

- model:

  For `coef` and `vcov`: `"full"` (both parts, the default), `"count"`
  or `"zero"`.

- type:

  For `predict`: `"response"` (the mean \\(1 - \pi)\mu\\, the default),
  `"count"` (the count-model mean \\\mu\\) or `"zero"` (the
  zero-inflation probability \\\pi\\). For `residuals`: `"response"`
  (`y - (1 - pi) mu`, the default) or `"pearson"` (the response residual
  divided by the ZINB standard deviation).

- newdata:

  An optional data frame with the covariates at which to predict. If
  omitted, predictions are for the data used in the fit.

## Value

- `print` returns its argument invisibly.

- `summary` returns an object of class `"summary.zinbreg"` with separate
  coefficient tables (estimates, standard errors, z values, p-values)
  for the count and zero parts; the count table includes `Log(theta)`.

- `coef` returns the coefficients, named `count_*` and `zero_*` when
  `model = "full"`, and `vcov` their variance-covariance matrix.

- `logLik` returns an object of class `"logLik"` whose `df` counts both
  sets of coefficients plus `theta`, so that
  [`AIC()`](https://rdrr.io/r/stats/AIC.html) and
  [`BIC()`](https://rdrr.io/r/stats/AIC.html) work.

- `fitted`, `residuals` and `predict` return numeric vectors.

- `nobs` returns the number of observations.

## Examples

``` r
set.seed(1)
n <- 500
d <- data.frame(x = rnorm(n), z = rnorm(n))
mu <- exp(1 + 0.5 * d$x)
p <- plogis(-1 + 0.8 * d$z)
d$y <- ifelse(runif(n) < p, 0, rnbinom(n, size = 2, mu = mu))
fit <- zinbreg(y ~ x | z, data = d)

fit
#> 
#> Call:
#> zinbreg(formula = y ~ x | z, data = d)
#> 
#> Count model coefficients (negative binomial with log link):
#> (Intercept)           x 
#>      1.1342      0.4486 
#> 
#> Zero-inflation model coefficients (binomial with logit link):
#> (Intercept)           z 
#>     -0.8172      0.7150 
#> 
#> Theta: 2.057
#> Log-likelihood: -944.6338 on 5 Df
#> 
summary(fit)
#> 
#> Call:
#> zinbreg(formula = y ~ x | z, data = d)
#> 
#> Count model coefficients (negative binomial with log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  1.13416    0.06076  18.666  < 2e-16 ***
#> x            0.44855    0.05387   8.326  < 2e-16 ***
#> Log(theta)   0.72138    0.19560   3.688 0.000226 ***
#> 
#> Zero-inflation model coefficients (binomial with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -0.8172     0.1672  -4.887 1.02e-06 ***
#> z             0.7150     0.1393   5.133 2.86e-07 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Theta: 2.057  (SE 0.4024)
#> Log-likelihood: -944.6338 on 5 Df
#> AIC: 1899.268   Number of observations: 500
#> Number of EM iterations: 36
#> 
coef(fit)
#> count_(Intercept)           count_x  zero_(Intercept)            zero_z 
#>         1.1341602         0.4485549        -0.8171633         0.7149913 
coef(fit, model = "zero")
#> (Intercept)           z 
#>  -0.8171633   0.7149913 
AIC(fit)
#> [1] 1899.268
BIC(fit)
#> [1] 1920.341
nd <- data.frame(x = c(-1, 0, 1), z = c(1, 0, -1))
predict(fit, newdata = nd)
#>        1        2        3 
#> 1.043145 2.156204 4.003188 
predict(fit, newdata = nd, type = "count")
#>        1        2        3 
#> 1.984973 3.108562 4.868155 
predict(fit, newdata = nd, type = "zero")
#>         1         2         3 
#> 0.4744792 0.3063661 0.1776787 
```
