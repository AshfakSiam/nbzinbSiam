# Negative binomial regression

Fits a negative binomial (NB2) regression model with log link, similar
to [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html). The
coefficients are estimated by iteratively reweighted least squares
(IRLS) and the dispersion parameter `theta` by Newton-Raphson on the log
scale; the two steps are alternated until the log-likelihood converges.
No general-purpose optimizer is used.

## Usage

``` r
nbreg(formula, data, offset, control = nbreg_control())
```

## Arguments

- formula:

  A model formula such as `y ~ x1 + x2`. An offset can be included as
  `offset(log(exposure))`.

- data:

  A data frame containing the variables in `formula`.

- offset:

  An optional numeric vector with an offset on the log scale (for
  example `log(exposure)`).

- control:

  A list of control parameters from
  [`nbreg_control()`](https://ashfaksiam.github.io/nbzinbSiam/reference/nbreg_control.md).

## Value

An object of class `"nbreg"`: a list containing

- coefficients:

  Estimated regression coefficients.

- vcov:

  Variance-covariance matrix of the coefficients.

- theta:

  Estimated dispersion parameter.

- SE.theta:

  Standard error of `theta`.

- loglik:

  Log-likelihood at the estimates.

- iterations:

  Number of outer iterations used.

- converged:

  `TRUE` if the algorithm converged.

- fitted.values:

  Fitted means.

- linear.predictors:

  Linear predictor, including the offset.

- call:

  The matched call.

together with the response `y`, the `offset`, the model frame `model`,
the `terms`, and the numbers `n` and `df.residual`.

## Details

The model is \\E(Y_i) = \mu_i\\, \\\log \mu_i = x_i^\top \beta + o_i\\
and \\Var(Y_i) = \mu_i + \mu_i^2 / \theta\\, where \\o_i\\ is an
optional offset. Standard errors of the coefficients come from \\(X^\top
\hat{W} X)^{-1}\\ with `theta` treated as known, as in
[`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html).

## Examples

``` r
set.seed(1)
x <- rnorm(200)
y <- rnbinom(200, size = 2, mu = exp(0.5 + 0.8 * x))
d <- data.frame(y = y, x = x)

fit <- nbreg(y ~ x, data = d)
fit$coefficients
#> (Intercept)           x 
#>   0.5264662   0.7583828 
fit$theta
#> [1] 2.060895
```
