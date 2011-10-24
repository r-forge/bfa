\name{predict.bfa}
\alias{predict.bfa}
\title{Posterior predictive and univariate conditional posterior predictive distributions
Posterior predictive and univariate conditional posterior predictive distributions, currently
implemented only for factor copula models}
\usage{
  \method{predict}{bfa} (object, resp.var = NA, cond.vars =
  NA, numeric.as.factor = TRUE, ...)
}
\arguments{
  \item{object}{\code{bfa} model object}

  \item{resp.var}{Response variable (see details). If NA,
  return draws from the joint posterior predictive.
  Otherwise}

  \item{cond.vars}{Conditioning variables; either a list
  like list(X1=val1, X2=val2) with X1, X2 variables in the
  original data frame, or a P length vector with either the
  conditioning value or NA (for marginalized variables).
  Ignored if resp.var is NA}

  \item{numeric.as.factor}{Treat numeric variables as
  ordinal when conditioning}

  \item{...}{Ignored}
}
\description{
  Posterior predictive and univariate conditional posterior
  predictive distributions Posterior predictive and
  univariate conditional posterior predictive
  distributions, currently implemented only for factor
  copula models
}

