# Zero-inflated negative binomial regression

Fits a zero-inflated negative binomial (ZINB) regression model, similar
to `pscl::zeroinfl(..., dist = "negbin")`. The model is a mixture: with
probability \\\pi_i\\ observation \\i\\ is a structural zero, and with
probability \\1 - \pi_i\\ it comes from a negative binomial distribution
with mean \\\mu_i\\ and dispersion \\\theta\\, where \\\log \mu_i =
x_i^\top \beta + o_i\\ and \\\mathrm{logit}\\ \pi_i = z_i^\top \gamma\\.

## Usage

``` r
zinbreg(formula, data, offset, control = zinbreg_control())
```

## Arguments

- formula:

  A two-part formula `y ~ x1 + x2 | z1 + z2`: count-model covariates
  before the `|` and zero-model covariates after it. Without `|`, the
  same covariates are used in both parts. An offset for the count model
  can be included as `offset(log(exposure))` before the `|`.

- data:

  A data frame containing the variables in `formula`.

- offset:

  An optional numeric vector with an offset for the count model, on the
  log scale.

- control:

  A list of control parameters from
  [`zinbreg_control()`](https://ashfaksiam.github.io/nbzinbSiam/reference/zinbreg_control.md).

## Value

An object of class `"zinbreg"`: a list containing

- coefficients:

  A list with components `count` and `zero`.

- vcov:

  Variance-covariance matrix of all estimates \\(\beta, \gamma,
  \log\theta)\\, from the inverse of the observed information matrix.

- theta:

  Estimated dispersion parameter of the count model.

- SE.logtheta:

  Standard error of \\\log\theta\\.

- SE.theta:

  Standard error of \\\theta\\, by the delta method: \\\theta \times
  SE(\log\theta)\\.

- loglik:

  Observed-data log-likelihood at the estimates.

- loglik_path:

  Log-likelihood at the start and after every EM iteration (it never
  decreases).

- iterations:

  Number of EM iterations used.

- converged:

  `TRUE` if the algorithm converged.

- fitted.values:

  Fitted means \\(1 - \pi_i) \mu_i\\.

- mu:

  Count-model means \\\mu_i\\.

- pi:

  Zero-inflation probabilities \\\pi_i\\.

- call:

  The matched call.

together with the response `y`, the count-model `offset`, the model
frame `model`, the `terms`, and the numbers `n` and `df.residual`.

## Details

The model is fitted by an EM algorithm. The E-step computes the
posterior probability that each zero is structural. The M-step fits a
logistic regression for \\\gamma\\ and a weighted negative binomial
regression for \\(\beta, \theta)\\, both by IRLS (with Newton-Raphson
for \\\theta\\). No general-purpose optimizer is used.

Standard errors are obtained by inverting the observed information
matrix of the observed-data log-likelihood with respect to \\(\beta,
\gamma, \log\theta)\\ at the estimates. By default its second
derivatives are computed analytically; see
[`zinbreg_control()`](https://ashfaksiam.github.io/nbzinbSiam/reference/zinbreg_control.md).

## Examples

``` r
set.seed(1)
n <- 500
d <- data.frame(x = rnorm(n), z = rnorm(n))
mu <- exp(1 + 0.5 * d$x)
p <- plogis(-1 + 0.8 * d$z)
d$y <- ifelse(runif(n) < p, 0, rnbinom(n, size = 2, mu = mu))

fit <- zinbreg(y ~ x | z, data = d)
fit$coefficients
#> $count
#> (Intercept)           x 
#>   1.1341602   0.4485549 
#> 
#> $zero
#> (Intercept)           z 
#>  -0.8171633   0.7149913 
#> 
fit$theta
#> [1] 2.05727
fit$iterations
#> [1] 36
```
