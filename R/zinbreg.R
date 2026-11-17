#' Control parameters for zinbreg()
#'
#' @param tol Convergence tolerance. The EM algorithm stops when the
#'   observed-data log-likelihood changes by less than `tol` between two
#'   iterations. The IRLS and Newton-Raphson loops inside each M-step use the
#'   same tolerance.
#' @param maxit Maximum number of EM iterations.
#' @param trace If `TRUE`, print the log-likelihood and theta at every EM
#'   iteration.
#'
#' @return A list with components `tol`, `maxit` and `trace`.
#' @examples
#' zinbreg_control(tol = 1e-10, maxit = 5000)
#' @export
zinbreg_control <- function(tol = 1e-8, maxit = 1000, trace = FALSE) {
  if (!is.numeric(tol) || length(tol) != 1 || tol <= 0) {
    stop("'tol' must be a single positive number.")
  }
  if (!is.numeric(maxit) || length(maxit) != 1 || maxit < 1) {
    stop("'maxit' must be a single positive integer.")
  }
  list(tol = tol, maxit = as.integer(maxit), trace = isTRUE(trace))
}

#' ZINB fitting engine: EM algorithm with IRLS M-steps
#'
#' Implements Section 3.2 of the project description.
#'
#' Starting values: an NB fit for (beta, theta), and a logistic regression
#' of I(y = 0) on Z for gamma. Then repeat
#'   E-step:            tau_i = P(structural zero | y_i)
#'   M-step, zero part: logistic IRLS of tau on Z
#'   M-step, count part: NB fit with prior weights 1 - tau
#' until the observed-data log-likelihood changes by less than tol.
#'
#' @param X Count-model design matrix.
#' @param Z Zero-model design matrix.
#' @param y Vector of counts.
#' @param offset Count-model offset.
#' @param control A list from zinbreg_control().
#' @return A list with the estimates, the log-likelihood at every iteration
#'   and convergence information.
#' @noRd
zinb_fit <- function(X, Z, y, offset = rep(0, length(y)),
                     control = zinbreg_control()) {
  tol <- control$tol
  inner <- nbreg_control(tol = tol, maxit = 100)

  # Starting values.
  nb0 <- nb_fit(X, y, offset = offset, control = inner)
  beta <- nb0$coefficients
  theta <- nb0$theta
  gamma <- logit_irls(Z, as.numeric(y == 0), tol = tol,
                      maxit = 100)$coefficients

  mu <- exp(drop(X %*% beta) + offset)
  eta_zero <- drop(Z %*% gamma)
  ll <- zinb_loglik(y, mu, eta_zero, theta)
  ll_path <- ll

  converged <- FALSE
  iter <- 0
  while (iter < control$maxit) {
    iter <- iter + 1

    # E-step.
    tau <- zinb_estep(y, mu, eta_zero, theta)

    # M-step, zero part: logistic regression with fractional response tau.
    gamma <- logit_irls(Z, tau, start = gamma, tol = tol,
                        maxit = 100)$coefficients
    eta_zero <- drop(Z %*% gamma)

    # M-step, count part: the NB fit of Section 3.1 with prior weights
    # 1 - tau (structural zeros carry no information about the count model).
    nb <- nb_fit(X, y, offset = offset, weights = 1 - tau, control = inner,
                 start = beta, theta_start = theta)
    beta <- nb$coefficients
    theta <- nb$theta
    mu <- nb$mu

    ll_new <- zinb_loglik(y, mu, eta_zero, theta)
    ll_path <- c(ll_path, ll_new)
    if (control$trace) {
      cat(sprintf("EM iter %4d  logLik = %.10f  theta = %.6f\n",
                  iter, ll_new, theta))
    }

    done <- abs(ll_new - ll) < tol
    ll <- ll_new
    if (done) {
      converged <- TRUE
      break
    }
  }

  names(beta) <- colnames(X)
  names(gamma) <- colnames(Z)

  # Standard errors from the observed information matrix.
  vcov <- zinb_vcov(beta, gamma, theta, X, Z, y, offset)

  list(count = beta, zero = gamma, theta = theta, vcov = vcov, loglik = ll,
       loglik_path = ll_path, iter = iter, converged = converged,
       mu = mu, pi = stats::plogis(eta_zero), eta_zero = eta_zero,
       theta_large = theta >= 1e4)
}

#' Zero-inflated negative binomial regression
#'
#' Fits a zero-inflated negative binomial (ZINB) regression model, similar
#' to `pscl::zeroinfl(..., dist = "negbin")`. The model is a mixture: with
#' probability \eqn{\pi_i} observation \eqn{i} is a structural zero, and with
#' probability \eqn{1 - \pi_i} it comes from a negative binomial distribution
#' with mean \eqn{\mu_i} and dispersion \eqn{\theta}, where
#' \eqn{\log \mu_i = x_i^\top \beta + o_i} and
#' \eqn{\mathrm{logit}\, \pi_i = z_i^\top \gamma}.
#'
#' The model is fitted by an EM algorithm. The E-step computes the posterior
#' probability that each zero is structural. The M-step fits a logistic
#' regression for \eqn{\gamma} and a weighted negative binomial regression
#' for \eqn{(\beta, \theta)}, both by IRLS (with Newton-Raphson for
#' \eqn{\theta}). No general-purpose optimizer is used.
#'
#' Standard errors are obtained by inverting the observed information
#' matrix of the observed-data log-likelihood with respect to
#' \eqn{(\beta, \gamma, \log\theta)} at the estimates.
#'
#' @param formula A two-part formula `y ~ x1 + x2 | z1 + z2`: count-model
#'   covariates before the `|` and zero-model covariates after it. Without
#'   `|`, the same covariates are used in both parts. An offset for the
#'   count model can be included as `offset(log(exposure))` before the `|`.
#' @param data A data frame containing the variables in `formula`.
#' @param offset An optional numeric vector with an offset for the count
#'   model, on the log scale.
#' @param control A list of control parameters from [zinbreg_control()].
#'
#' @return An object of class `"zinbreg"`: a list containing
#' \describe{
#'   \item{coefficients}{A list with components `count` and `zero`.}
#'   \item{vcov}{Variance-covariance matrix of all estimates
#'     \eqn{(\beta, \gamma, \log\theta)}, from the inverse of the observed
#'     information matrix.}
#'   \item{theta}{Estimated dispersion parameter of the count model.}
#'   \item{SE.logtheta}{Standard error of \eqn{\log\theta}.}
#'   \item{SE.theta}{Standard error of \eqn{\theta}, by the delta method:
#'     \eqn{\theta \times SE(\log\theta)}.}
#'   \item{loglik}{Observed-data log-likelihood at the estimates.}
#'   \item{loglik_path}{Log-likelihood at the start and after every EM
#'     iteration (it never decreases).}
#'   \item{iterations}{Number of EM iterations used.}
#'   \item{converged}{`TRUE` if the algorithm converged.}
#'   \item{fitted.values}{Fitted means \eqn{(1 - \pi_i) \mu_i}.}
#'   \item{mu}{Count-model means \eqn{\mu_i}.}
#'   \item{pi}{Zero-inflation probabilities \eqn{\pi_i}.}
#'   \item{call}{The matched call.}
#' }
#' together with the response `y`, the count-model `offset`, the model
#' frame `model`, the `terms`, and the numbers `n` and `df.residual`.
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' d <- data.frame(x = rnorm(n), z = rnorm(n))
#' mu <- exp(1 + 0.5 * d$x)
#' p <- plogis(-1 + 0.8 * d$z)
#' d$y <- ifelse(runif(n) < p, 0, rnbinom(n, size = 2, mu = mu))
#'
#' fit <- zinbreg(y ~ x | z, data = d)
#' fit$coefficients
#' fit$theta
#' fit$iterations
#' @export
zinbreg <- function(formula, data, offset, control = zinbreg_control()) {
  cl <- match.call()
  if (missing(data)) data <- environment(formula)

  # Two-part formula. Without "|", use the same covariates in both parts.
  f <- formula
  rhs <- f[[3]]
  if (!(is.call(rhs) && identical(rhs[[1]], as.name("|")))) {
    f[[3]] <- call("|", rhs, rhs)
  }
  F <- Formula::as.Formula(f)

  # Model frame for all variables in both parts.
  mf <- match.call(expand.dots = FALSE)
  keep <- match(c("formula", "data", "offset"), names(mf), 0L)
  mf <- mf[c(1L, keep)]
  mf$formula <- F
  mf$data <- data
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  y <- stats::model.response(mf)
  mtX <- stats::terms(F, data = data, rhs = 1L)
  mtZ <- stats::delete.response(stats::terms(F, data = data, rhs = 2L))
  X <- stats::model.matrix(mtX, mf)
  Z <- stats::model.matrix(mtZ, mf)
  off <- count_offset(mtX, mf)

  # Informative errors for bad input.
  if (is.null(y) || !is.numeric(y) || is.matrix(y)) {
    stop("The response must be a numeric vector of counts.")
  }
  if (any(y < 0) || any(abs(y - round(y)) > 1e-8)) {
    stop("The response must contain non-negative integers (counts).")
  }
  if (!any(y == 0)) {
    stop("The response has no zeros, so a zero-inflated model cannot be fitted.")
  }
  if (length(y) <= ncol(X) + ncol(Z) + 1) {
    stop("There are not enough observations to fit this model.")
  }
  if (!is.list(control) || is.null(control$tol) || is.null(control$maxit)) {
    stop("'control' must be created by zinbreg_control().")
  }

  fit <- zinb_fit(X, Z, y, offset = off, control = control)

  if (!fit$converged) {
    warning("zinbreg() did not converge in ", control$maxit,
            " EM iterations; consider increasing 'maxit'.", call. = FALSE)
  }
  if (fit$theta_large) {
    warning("theta is very large: the count model shows no overdispersion ",
            "and is essentially a Poisson model.", call. = FALSE)
  }

  se_logtheta <- sqrt(fit$vcov["Log(theta)", "Log(theta)"])
  out <- list(coefficients = list(count = fit$count, zero = fit$zero),
              vcov = fit$vcov,
              theta = fit$theta,
              SE.logtheta = se_logtheta,
              SE.theta = fit$theta * se_logtheta,
              loglik = fit$loglik,
              loglik_path = fit$loglik_path,
              iterations = fit$iter,
              converged = fit$converged,
              fitted.values = (1 - fit$pi) * fit$mu,
              mu = fit$mu,
              pi = fit$pi,
              y = y,
              offset = off,
              n = length(y),
              df.residual = length(y) - ncol(X) - ncol(Z) - 1,
              call = cl,
              formula = F,
              terms = list(count = mtX, zero = mtZ),
              model = mf,
              control = control)
  class(out) <- "zinbreg"
  out
}

#' Count-model offset from a model frame
#'
#' Adds up the offset() terms of the count-model formula and the offset
#' given through the `offset` argument (stored by model.frame() as the
#' column "(offset)").
#'
#' @param mt Terms of the count model.
#' @param mf Model frame.
#' @return A numeric vector of length nrow(mf).
#' @noRd
count_offset <- function(mt, mf) {
  off <- rep(0, nrow(mf))
  idx <- attr(mt, "offset")
  if (!is.null(idx)) {
    vars <- vapply(as.list(attr(mt, "variables"))[-1L],
                   function(v) paste(deparse(v), collapse = ""), "")
    for (v in vars[idx]) off <- off + mf[[v]]
  }
  if (!is.null(mf[["(offset)"]])) off <- off + mf[["(offset)"]]
  off
}
