# Internal helpers shared by nbreg() and zinbreg().
# None of these are exported; users never call them directly.

#' Negative binomial (NB2) log-likelihood
#'
#' Implements the log-likelihood in Section 2.1 of the project description:
#' sum_i a_i [ lgamma(y_i + theta) - lgamma(theta) - log(y_i!)
#'             + theta log(theta) + y_i log(mu_i)
#'             - (theta + y_i) log(theta + mu_i) ].
#'
#' @param y Vector of observed counts.
#' @param mu Vector of fitted means (same length as y, all > 0).
#' @param theta Dispersion parameter (a single number > 0).
#' @param weights Prior weights a_i (default: all 1).
#' @return A single number: the weighted log-likelihood.
#' @noRd
nb_loglik <- function(y, mu, theta, weights = rep(1, length(y))) {
  sum(weights * (lgamma(y + theta) - lgamma(theta) - lgamma(y + 1) +
                   theta * log(theta) + y * log(mu) -
                   (theta + y) * log(theta + mu)))
}

#' Weighted least squares
#'
#' Solves beta = (X' W X)^{-1} X' W z, the update used in every IRLS step,
#' where W = diag(w).
#'
#' @param X Design matrix (n rows, p columns).
#' @param z Working response (length n).
#' @param w Working weights (length n, all >= 0).
#' @return The coefficient vector (length p).
#' @noRd
wls_solve <- function(X, z, w) {
  XtWX <- crossprod(X, w * X)   # X' W X  (p x p)
  XtWz <- crossprod(X, w * z)   # X' W z  (p x 1)
  drop(solve(XtWX, XtWz))
}

