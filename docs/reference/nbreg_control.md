# Control parameters for nbreg()

Control parameters for nbreg()

## Usage

``` r
nbreg_control(tol = 1e-08, maxit = 100, trace = FALSE)
```

## Arguments

- tol:

  Convergence tolerance. The outer loop stops when the log-likelihood
  changes by less than `tol`. The inner IRLS and Newton-Raphson loops
  use the same tolerance for the change in the parameters.

- maxit:

  Maximum number of iterations, for the outer loop and for each inner
  loop.

- trace:

  If `TRUE`, print the log-likelihood and theta at every outer
  iteration.

## Value

A list with components `tol`, `maxit` and `trace`.

## Examples

``` r
nbreg_control(tol = 1e-10, maxit = 200)
#> $tol
#> [1] 1e-10
#> 
#> $maxit
#> [1] 200
#> 
#> $trace
#> [1] FALSE
#> 
```
