#' Control parameters for nbreg()
#'
#' @param tol Convergence tolerance. The outer loop stops when the
#'   log-likelihood changes by less than `tol`. The inner IRLS and
#'   Newton-Raphson loops use the same tolerance for the change in the
#'   parameters.
#' @param maxit Maximum number of iterations, for the outer loop and for
#'   each inner loop.
#' @param trace If `TRUE`, print the log-likelihood and theta at every outer
#'   iteration.
#'
#' @return A list with components `tol`, `maxit` and `trace`.
#' @examples
#' nbreg_control(tol = 1e-10, maxit = 200)
#' @export
nbreg_control <- function(tol = 1e-8, maxit = 100, trace = FALSE) {
  if (!is.numeric(tol) || length(tol) != 1 || tol <= 0) {
    stop("'tol' must be a single positive number.")
  }
  if (!is.numeric(maxit) || length(maxit) != 1 || maxit < 1) {
    stop("'maxit' must be a single positive integer.")
  }
  list(tol = tol, maxit = as.integer(maxit), trace = isTRUE(trace))
}

#' NB fitting engine: alternate IRLS for beta and Newton-Raphson for theta
#'
#' Implements the outer loop in Section 3.1(c) of the project description
#' and the standard errors in Section 3.1(d). It works on a design matrix,
#' so it can be reused by the ZINB M-step with prior weights 1 - tau.
#'
#' @param X Design matrix.
#' @param y Vector of counts.
#' @param offset Offset vector.
#' @param weights Prior weights a_i.
#' @param control A list from nbreg_control().
#' @param start Optional starting value for beta.
#' @param theta_start Optional starting value for theta.
#' @return A list with the estimates, their variances, the log-likelihood
#'   and convergence information.
#' @noRd
nb_fit <- function(X, y, offset = rep(0, length(y)),
                   weights = rep(1, length(y)), control = nbreg_control(),
                   start = NULL, theta_start = NULL) {
  tol <- control$tol
  maxit <- control$maxit
  theta_max <- 1e7    # numerical cap: beyond this the model is Poisson
  theta_warn <- 1e4   # above this the fit is flagged as "Poisson limit"

  # Starting values (Section 3.1(c)): a Poisson fit for beta, obtained as
  # the NB IRLS with a huge theta, then a moment estimate for theta.
  if (is.null(start)) {
    start <- nb_irls(X, y, theta = theta_max, offset = offset,
                     weights = weights, tol = tol, maxit = maxit)$coefficients
  }
  beta <- start
  mu <- exp(drop(X %*% beta) + offset)
  theta <- theta_start
  if (is.null(theta)) {
    theta <- sum(weights) / sum(weights * (y / mu - 1)^2)
  }
  theta <- min(theta, theta_max)
  ll <- nb_loglik(y, mu, theta, weights)

  converged <- FALSE
  iter <- 0
  while (iter < maxit) {
    iter <- iter + 1

    # (a) IRLS for beta with theta held fixed.
    fit_beta <- nb_irls(X, y, theta, offset = offset, weights = weights,
                        start = beta, tol = tol, maxit = maxit)
    beta <- fit_beta$coefficients
    mu <- fit_beta$mu

    # (b) Newton-Raphson for theta with beta held fixed.
    theta <- nb_theta_newton(y, mu, weights = weights, start = theta,
                             tol = tol, maxit = maxit)$theta
    theta <- min(theta, theta_max)

    ll_new <- nb_loglik(y, mu, theta, weights)
    if (control$trace) {
      cat(sprintf("iter %3d  logLik = %.8f  theta = %.6f\n",
                  iter, ll_new, theta))
    }

    done <- abs(ll_new - ll) < tol
    ll <- ll_new
    if (done) {
      converged <- TRUE
      break
    }
  }

  # Standard errors (Section 3.1(d)):
  #   Var(beta) = (X' W X)^{-1} with theta treated as known,
  #   SE(theta) = sqrt(-1 / H(theta)).
  eta <- drop(X %*% beta) + offset
  w <- weights * theta * mu / (theta + mu)
  vcov <- solve(crossprod(X, w * X))
  dimnames(vcov) <- list(colnames(X), colnames(X))
  H <- nb_theta_hess(theta, y, mu, weights)
  se_theta <- if (is.finite(H) && H < 0) sqrt(-1 / H) else NA_real_

  list(coefficients = beta, theta = theta, se_theta = se_theta, vcov = vcov,
       loglik = ll, iter = iter, converged = converged,
       eta = eta, mu = mu, theta_large = theta >= theta_warn)
}

#' Negative binomial regression
#'
#' Fits a negative binomial (NB2) regression model with log link, similar to
#' `MASS::glm.nb()`. The coefficients are estimated by iteratively reweighted
#' least squares (IRLS) and the dispersion parameter `theta` by
#' Newton-Raphson on the log scale; the two steps are alternated until the
#' log-likelihood converges. No general-purpose optimizer is used.
#'
#' The model is \eqn{E(Y_i) = \mu_i}, \eqn{\log \mu_i = x_i^\top \beta + o_i}
#' and \eqn{Var(Y_i) = \mu_i + \mu_i^2 / \theta}, where \eqn{o_i} is an
#' optional offset. Standard errors of the coefficients come from
#' \eqn{(X^\top \hat{W} X)^{-1}} with `theta` treated as known, as in
#' `MASS::glm.nb()`.
#'
#' @param formula A model formula such as `y ~ x1 + x2`. An offset can be
#'   included as `offset(log(exposure))`.
#' @param data A data frame containing the variables in `formula`.
#' @param offset An optional numeric vector with an offset on the log scale
#'   (for example `log(exposure)`).
#' @param control A list of control parameters from [nbreg_control()].
#'
#' @return An object of class `"nbreg"`: a list containing
#' \describe{
#'   \item{coefficients}{Estimated regression coefficients.}
#'   \item{vcov}{Variance-covariance matrix of the coefficients.}
#'   \item{theta}{Estimated dispersion parameter.}
#'   \item{SE.theta}{Standard error of `theta`.}
#'   \item{loglik}{Log-likelihood at the estimates.}
#'   \item{iterations}{Number of outer iterations used.}
#'   \item{converged}{`TRUE` if the algorithm converged.}
#'   \item{fitted.values}{Fitted means.}
#'   \item{linear.predictors}{Linear predictor, including the offset.}
#'   \item{call}{The matched call.}
#' }
#' together with the response `y`, the `offset`, the model frame `model`,
#' the `terms`, and the numbers `n` and `df.residual`.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(200)
#' y <- rnbinom(200, size = 2, mu = exp(0.5 + 0.8 * x))
#' d <- data.frame(y = y, x = x)
#'
#' fit <- nbreg(y ~ x, data = d)
#' fit$coefficients
#' fit$theta
#' @export
nbreg <- function(formula, data, offset, control = nbreg_control()) {
  cl <- match.call()

  # Build the model frame the same way glm() does.
  mf <- match.call(expand.dots = FALSE)
  keep <- match(c("formula", "data", "offset"), names(mf), 0L)
  mf <- mf[c(1L, keep)]
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  y <- stats::model.response(mf)
  X <- stats::model.matrix(mt, mf)
  off <- stats::model.offset(mf)
  if (is.null(off)) off <- rep(0, NROW(X))

  # Informative errors for bad input.
  if (is.null(y) || !is.numeric(y) || is.matrix(y)) {
    stop("The response must be a numeric vector of counts.")
  }
  if (any(y < 0) || any(abs(y - round(y)) > 1e-8)) {
    stop("The response must contain non-negative integers (counts).")
  }
  if (length(y) <= ncol(X)) {
    stop("There are not enough observations to fit this model.")
  }
  if (!is.list(control) || is.null(control$tol) || is.null(control$maxit)) {
    stop("'control' must be created by nbreg_control().")
  }

  fit <- nb_fit(X, y, offset = off, control = control)

  if (!fit$converged) {
    warning("nbreg() did not converge in ", control$maxit,
            " iterations; consider increasing 'maxit'.", call. = FALSE)
  }
  if (fit$theta_large) {
    warning("theta is very large: the data show no overdispersion and the ",
            "fit is essentially a Poisson regression.", call. = FALSE)
  }

  out <- list(coefficients = fit$coefficients,
              vcov = fit$vcov,
              theta = fit$theta,
              SE.theta = fit$se_theta,
              loglik = fit$loglik,
              iterations = fit$iter,
              converged = fit$converged,
              fitted.values = fit$mu,
              linear.predictors = fit$eta,
              y = y,
              offset = off,
              n = length(y),
              df.residual = length(y) - ncol(X),
              call = cl,
              formula = formula,
              terms = mt,
              model = mf,
              control = control)
  class(out) <- "nbreg"
  out
}

