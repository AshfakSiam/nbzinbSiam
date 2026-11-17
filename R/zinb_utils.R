# Building blocks for zinbreg(): the ZINB log-likelihood, the E-step, and
# the weighted logistic IRLS used in the zero-part M-step.

#' ZINB observed-data log-likelihood
#'
#' Implements the log-likelihood in Section 2.2 of the project description:
#'   sum over y_i = 0:  log[ pi_i + (1 - pi_i) f_NB(0; mu_i, theta) ]
#'   sum over y_i > 0:  log(1 - pi_i) + log f_NB(y_i; mu_i, theta)
#' with f_NB(0; mu, theta) = (theta / (theta + mu))^theta.
#'
#' Probabilities are handled on the log scale so that the result stays
#' accurate when pi_i is very close to 0 or 1.
#'
#' @param y Vector of counts.
#' @param mu Count-model means (length n).
#' @param eta_zero Zero-model linear predictor z_i' gamma (length n), so
#'   that pi_i = plogis(eta_zero).
#' @param theta Dispersion parameter of the count model.
#' @return A single number.
#' @noRd
zinb_loglik <- function(y, mu, eta_zero, theta) {
  log_pi <- stats::plogis(eta_zero, log.p = TRUE)            # log(pi)
  log_1mpi <- stats::plogis(eta_zero, lower.tail = FALSE,
                            log.p = TRUE)                    # log(1 - pi)
  log_f <- stats::dnbinom(y, size = theta, mu = mu, log = TRUE)

  zero <- y == 0
  # log(pi + (1 - pi) f0) computed as log(exp(a) + exp(b)) without overflow.
  a <- log_pi[zero]
  b <- log_1mpi[zero] + log_f[zero]
  m <- pmax(a, b)
  ll_zero <- m + log(exp(a - m) + exp(b - m))

  sum(ll_zero) + sum(log_1mpi[!zero] + log_f[!zero])
}

#' E-step: posterior probability of a structural zero
#'
#' Implements the E-step in Section 3.2 of the project description:
#'   tau_i = pi_i / (pi_i + (1 - pi_i) f_NB(0; mu_i, theta))  if y_i = 0,
#'   tau_i = 0                                                 if y_i > 0.
#'
#' @inheritParams zinb_loglik
#' @return A vector of tau_i values in [0, 1].
#' @noRd
zinb_estep <- function(y, mu, eta_zero, theta) {
  pi0 <- stats::plogis(eta_zero)
  f0 <- (theta / (theta + mu))^theta
  tau <- pi0 / (pi0 + (1 - pi0) * f0)
  tau[y > 0] <- 0
  tau
}

#' IRLS for a logistic regression with a fractional response
#'
#' Implements the zero-part M-step in Section 3.2 of the project
#' description. It maximizes
#'   sum_i a_i [ tau_i log(pi_i) + (1 - tau_i) log(1 - pi_i) ]
#' over gamma, with logit(pi_i) = z_i' gamma, by repeated weighted least
#' squares with
#'   working weights   w_i = a_i * pi_i * (1 - pi_i)
#'   working response  zz_i = z_i' gamma + (tau_i - pi_i) / (pi_i (1 - pi_i)).
#' The response tau_i may be any number in [0, 1], not only 0 or 1.
#'
#' @param Z Design matrix of the zero model.
#' @param tau Response values in [0, 1].
#' @param weights Prior weights a_i (default: all 1).
#' @param start Optional starting value for gamma.
#' @param tol Convergence tolerance for the change in gamma.
#' @param maxit Maximum number of iterations.
#' @return A list with coefficients, eta, the number of iterations and a
#'   convergence flag.
#' @noRd
logit_irls <- function(Z, tau, weights = rep(1, length(tau)), start = NULL,
                       tol = 1e-8, maxit = 50) {
  # Objective: the expected complete-data log-likelihood of the zero part.
  objective <- function(eta) {
    sum(weights * (tau * stats::plogis(eta, log.p = TRUE) +
                     (1 - tau) * stats::plogis(eta, lower.tail = FALSE,
                                               log.p = TRUE)))
  }

  if (is.null(start)) start <- rep(0, ncol(Z))
  gamma <- start
  eta <- drop(Z %*% gamma)
  obj <- objective(eta)

  converged <- FALSE
  iter <- 0
  while (iter < maxit) {
    iter <- iter + 1

    # Working weights and working response. pi is kept away from exactly
    # 0 and 1 so that the division below is always defined.
    p <- stats::plogis(eta)
    p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
    w <- weights * p * (1 - p)
    zz <- eta + (tau - p) / (p * (1 - p))

    gamma_new <- wls_solve(Z, zz, w)
    eta_new <- drop(Z %*% gamma_new)
    obj_new <- objective(eta_new)

    # Step-halving, as in nb_irls().
    half <- 0
    while ((!is.finite(obj_new) || obj_new < obj - 1e-8 * abs(obj)) &&
           half < 10) {
      half <- half + 1
      gamma_new <- (gamma_new + gamma) / 2
      eta_new <- drop(Z %*% gamma_new)
      obj_new <- objective(eta_new)
    }

    change <- max(abs(gamma_new - gamma))
    gamma <- gamma_new
    eta <- eta_new
    obj <- obj_new

    if (change <= tol * (1 + max(abs(gamma)))) {
      converged <- TRUE
      break
    }
  }

  names(gamma) <- colnames(Z)
  list(coefficients = gamma, eta = eta, iter = iter, converged = converged)
}
