# Control parameters for zinbreg()

Control parameters for zinbreg()

## Usage

``` r
zinbreg_control(
  tol = 1e-08,
  maxit = 1000,
  trace = FALSE,
  hessian = c("analytic", "numerical")
)
```

## Arguments

- tol:

  Convergence tolerance. The EM algorithm stops when the observed-data
  log-likelihood changes by less than `tol` between two iterations. The
  IRLS and Newton-Raphson loops inside each M-step use the same
  tolerance.

- maxit:

  Maximum number of EM iterations.

- trace:

  If `TRUE`, print the log-likelihood and theta at every EM iteration.

- hessian:

  How the Hessian for the standard errors is computed: `"analytic"` (the
  default) uses the second derivatives derived by hand; `"numerical"`
  uses
  [`numDeriv::hessian()`](https://rdrr.io/pkg/numDeriv/man/hessian.html)
  on the log-likelihood.

## Value

A list with components `tol`, `maxit`, `trace` and `hessian`.

## Examples

``` r
zinbreg_control(tol = 1e-10, maxit = 5000)
#> $tol
#> [1] 1e-10
#> 
#> $maxit
#> [1] 5000
#> 
#> $trace
#> [1] FALSE
#> 
#> $hessian
#> [1] "analytic"
#> 
zinbreg_control(hessian = "numerical")
#> $tol
#> [1] 1e-08
#> 
#> $maxit
#> [1] 1000
#> 
#> $trace
#> [1] FALSE
#> 
#> $hessian
#> [1] "numerical"
#> 
```
