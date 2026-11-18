# Methods for nbreg objects

Standard S3 methods for models fitted by
[`nbreg()`](https://ashfaksiam.github.io/nbzinbSiam/reference/nbreg.md).

## Usage

``` r
# S3 method for class 'nbreg'
print(x, digits = max(3, getOption("digits") - 3), ...)

# S3 method for class 'nbreg'
summary(object, ...)

# S3 method for class 'summary.nbreg'
print(x, digits = max(3, getOption("digits") - 3), ...)

# S3 method for class 'nbreg'
coef(object, ...)

# S3 method for class 'nbreg'
vcov(object, ...)

# S3 method for class 'nbreg'
logLik(object, ...)

# S3 method for class 'nbreg'
nobs(object, ...)

# S3 method for class 'nbreg'
fitted(object, ...)

# S3 method for class 'nbreg'
residuals(object, type = c("response", "pearson"), ...)

# S3 method for class 'nbreg'
predict(object, newdata, type = c("response", "link"), ...)
```

## Arguments

- x, object:

  An object of class `"nbreg"` (or `"summary.nbreg"` for the `print`
  method of a summary).

- digits:

  Number of significant digits to print.

- ...:

  Further arguments (currently ignored).

- type:

  For `predict`: `"response"` (the mean, the default) or `"link"` (the
  linear predictor). For `residuals`: `"response"` (`y - mu`, the
  default) or `"pearson"` (`(y - mu) / sqrt(mu + mu^2 / theta)`).

- newdata:

  An optional data frame with the covariates at which to predict. If
  omitted, predictions are for the data used in the fit.

## Value

- `print` returns its argument invisibly.

- `summary` returns an object of class `"summary.nbreg"` whose
  `coefficients` component is a table of estimates, standard errors, z
  values and p-values.

- `coef` returns the coefficient vector and `vcov` its
  variance-covariance matrix.

- `logLik` returns an object of class `"logLik"` whose `df` counts the
  regression coefficients plus `theta`, so that
  [`AIC()`](https://rdrr.io/r/stats/AIC.html) and
  [`BIC()`](https://rdrr.io/r/stats/AIC.html) work.

- `fitted`, `residuals` and `predict` return numeric vectors.

- `nobs` returns the number of observations.

## Examples

``` r
set.seed(1)
d <- data.frame(x = rnorm(200))
d$y <- rnbinom(200, size = 2, mu = exp(0.5 + 0.8 * d$x))
fit <- nbreg(y ~ x, data = d)

fit
#> 
#> Call:
#> nbreg(formula = y ~ x, data = d)
#> 
#> Coefficients:
#> (Intercept)           x 
#>      0.5265      0.7584 
#> 
#> Theta: 2.061  (SE 0.4497)
#> Log-likelihood: -360.7766 on 3 Df
#> 
summary(fit)
#> 
#> Call:
#> nbreg(formula = y ~ x, data = d)
#> 
#> Coefficients (negative binomial with log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.52647    0.07790   6.758  1.4e-11 ***
#> x            0.75838    0.08003   9.477  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Theta: 2.061  (SE 0.4497)
#> Log-likelihood: -360.7766 on 3 Df
#> AIC: 727.5532   Number of observations: 200
#> Number of outer iterations: 3
#> 
coef(fit)
#> (Intercept)           x 
#>   0.5264662   0.7583828 
vcov(fit)
#>              (Intercept)            x
#> (Intercept)  0.006068088 -0.002269413
#> x           -0.002269413  0.006404308
logLik(fit)
#> 'log Lik.' -360.7766 (df=3)
AIC(fit)
#> [1] 727.5532
BIC(fit)
#> [1] 737.4481
nobs(fit)
#> [1] 200
head(fitted(fit))
#>         1         2         3         4         5         6 
#> 1.0527170 1.9459263 0.8982919 5.6762999 2.1735443 0.9086794 
head(residuals(fit, type = "pearson"))
#>           1           2           3           4           5           6 
#>  1.54407703  1.05604101  0.97006203  3.10283516 -0.08212136 -0.79412006 
predict(fit, newdata = data.frame(x = c(-1, 0, 1)))
#>         1         2         3 
#> 0.7930123 1.6929391 3.6141219 
```
