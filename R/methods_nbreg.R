#' Methods for nbreg objects
#'
#' Standard S3 methods for models fitted by [nbreg()].
#'
#' @param x,object An object of class `"nbreg"` (or `"summary.nbreg"` for
#'   the `print` method of a summary).
#' @param newdata An optional data frame with the covariates at which to
#'   predict. If omitted, predictions are for the data used in the fit.
#' @param type For `predict`: `"response"` (the mean, the default) or
#'   `"link"` (the linear predictor). For `residuals`: `"response"`
#'   (`y - mu`, the default) or `"pearson"`
#'   (`(y - mu) / sqrt(mu + mu^2 / theta)`).
#' @param digits Number of significant digits to print.
#' @param ... Further arguments (currently ignored).
#'
#' @return
#' * `print` returns its argument invisibly.
#' * `summary` returns an object of class `"summary.nbreg"` whose
#'   `coefficients` component is a table of estimates, standard errors,
#'   z values and p-values.
#' * `coef` returns the coefficient vector and `vcov` its
#'   variance-covariance matrix.
#' * `logLik` returns an object of class `"logLik"` whose `df` counts the
#'   regression coefficients plus `theta`, so that `AIC()` and `BIC()` work.
#' * `fitted`, `residuals` and `predict` return numeric vectors.
#' * `nobs` returns the number of observations.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(200))
#' d$y <- rnbinom(200, size = 2, mu = exp(0.5 + 0.8 * d$x))
#' fit <- nbreg(y ~ x, data = d)
#'
#' fit
#' summary(fit)
#' coef(fit)
#' vcov(fit)
#' logLik(fit)
#' AIC(fit)
#' BIC(fit)
#' nobs(fit)
#' head(fitted(fit))
#' head(residuals(fit, type = "pearson"))
#' predict(fit, newdata = data.frame(x = c(-1, 0, 1)))
#'
#' @name nbreg-methods
NULL

#' @rdname nbreg-methods
#' @export
print.nbreg <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(format(x$coefficients, digits = digits), quote = FALSE)
  cat("\nTheta: ", format(x$theta, digits = digits),
      "  (SE ", format(x$SE.theta, digits = digits), ")\n", sep = "")
  cat("Log-likelihood: ", format(x$loglik, digits = digits + 3),
      " on ", length(x$coefficients) + 1, " Df\n\n", sep = "")
  invisible(x)
}

#' @rdname nbreg-methods
#' @export
summary.nbreg <- function(object, ...) {
  est <- object$coefficients
  se <- sqrt(diag(object$vcov))
  z <- est / se
  p <- 2 * stats::pnorm(-abs(z))
  tab <- cbind(Estimate = est, `Std. Error` = se,
               `z value` = z, `Pr(>|z|)` = p)

  df <- length(est) + 1
  out <- list(call = object$call,
              coefficients = tab,
              theta = object$theta,
              SE.theta = object$SE.theta,
              loglik = object$loglik,
              df = df,
              aic = -2 * object$loglik + 2 * df,
              n = object$n,
              iterations = object$iterations,
              converged = object$converged)
  class(out) <- "summary.nbreg"
  out
}

#' @rdname nbreg-methods
#' @export
print.summary.nbreg <- function(x, digits = max(3, getOption("digits") - 3),
                                ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nCoefficients (negative binomial with log link):\n")
  stats::printCoefmat(x$coefficients, digits = digits,
                      signif.legend = TRUE)
  cat("\nTheta: ", format(x$theta, digits = digits),
      "  (SE ", format(x$SE.theta, digits = digits), ")\n", sep = "")
  cat("Log-likelihood: ", format(x$loglik, digits = digits + 3),
      " on ", x$df, " Df\n", sep = "")
  cat("AIC: ", format(x$aic, digits = digits + 3),
      "   Number of observations: ", x$n, "\n", sep = "")
  cat("Number of outer iterations: ", x$iterations,
      if (!x$converged) "  (did NOT converge)", "\n\n", sep = "")
  invisible(x)
}

#' @rdname nbreg-methods
#' @importFrom stats coef
#' @export
coef.nbreg <- function(object, ...) {
  object$coefficients
}

#' @rdname nbreg-methods
#' @importFrom stats vcov
#' @export
vcov.nbreg <- function(object, ...) {
  object$vcov
}

#' @rdname nbreg-methods
#' @importFrom stats logLik
#' @export
logLik.nbreg <- function(object, ...) {
  structure(object$loglik,
            df = length(object$coefficients) + 1,
            nobs = object$n,
            class = "logLik")
}

#' @rdname nbreg-methods
#' @importFrom stats nobs
#' @export
nobs.nbreg <- function(object, ...) {
  object$n
}

#' @rdname nbreg-methods
#' @importFrom stats fitted
#' @export
fitted.nbreg <- function(object, ...) {
  object$fitted.values
}

#' @rdname nbreg-methods
#' @importFrom stats residuals
#' @export
residuals.nbreg <- function(object, type = c("response", "pearson"), ...) {
  type <- match.arg(type)
  y <- object$y
  mu <- object$fitted.values
  res <- y - mu
  if (type == "pearson") {
    res <- res / sqrt(mu + mu^2 / object$theta)
  }
  res
}

#' @rdname nbreg-methods
#' @importFrom stats predict
#' @export
predict.nbreg <- function(object, newdata, type = c("response", "link"),
                          ...) {
  type <- match.arg(type)

  if (missing(newdata) || is.null(newdata)) {
    eta <- object$linear.predictors
  } else {
    # Rebuild the design matrix for the new data, using the factor levels
    # and contrasts of the original fit.
    tt <- stats::delete.response(object$terms)
    xlev <- stats::.getXlevels(object$terms, object$model)
    X_old <- stats::model.matrix(object$terms, object$model)
    mf <- stats::model.frame(tt, newdata, xlev = xlev)
    X <- stats::model.matrix(tt, mf, contrasts.arg = attr(X_old, "contrasts"))

    # Offset: from offset() terms in the formula and/or the offset argument.
    off <- stats::model.offset(mf)
    if (is.null(off)) off <- rep(0, nrow(X))
    if (!is.null(object$call$offset)) {
      off <- off + eval(object$call$offset, newdata)
    }
    eta <- drop(X %*% object$coefficients) + off
  }

  if (type == "response") exp(eta) else eta
}

