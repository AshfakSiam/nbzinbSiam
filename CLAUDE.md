# nbzinbSiam

STA 651 final project (Fall 2026, Dr. Donghyung Lee). An R package that fits
negative binomial (NB) and zero-inflated negative binomial (ZINB) regression
with hand-written estimation code.

## Hard rules for code in R/

The estimation must be implemented by hand:

- NB: IRLS for beta, alternated with Newton-Raphson on log(theta).
- ZINB: EM algorithm whose M-steps are IRLS fits.

PROHIBITED inside R/ (fitting code):

- optim, optimize, optimise, nlm, nlminb, uniroot
- glm, glm.fit, lm, lm.fit, lm.wfit
- MASS::glm.nb, MASS::theta.ml, any MASS fitting code
- pscl::zeroinfl, countreg, glmmTMB, gamlss, VGAM
- optimx, maxLik, stats4::mle, bbmle, or any other optimizer

ALLOWED inside R/:

- solve, qr, chol, crossprod, %*%
- lgamma, digamma, trigamma, dnbinom, plogis, qlogis
- model.frame, model.matrix, model.response, model.offset
- the Formula package for two-part formulas
- numDeriv::hessian, for standard errors only

The reference functions (glm, glm.nb, zeroinfl) may and should be used in
tests/, inst/benchmark/ and the report, never in R/.

## Conventions

- Exported functions: nbreg(), zinbreg(), nbreg_control(), zinbreg_control().
- S3 classes: "nbreg" and "zinbreg".
- Document with roxygen2; never edit man/ or NAMESPACE by hand.
- Tests live in tests/testthat/ and compare against MASS and pscl, guarded
  with skip_if_not_installed().
- Small commits with clear messages.

