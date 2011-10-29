\name{bfa_copula}
\alias{bfa_copula}
\title{Initialize and fit a copula factor model}
\usage{
  bfa_copula(x, data = NULL, num.factor = 1, restrict = NA,
  normal.dist = NA, center.data = TRUE, scale.data = TRUE,
  nsim = 0, nburn = 0, thin = 1, print.status = 500,
  keep.scores = FALSE, keep.loadings = TRUE, loading.prior
  = c("gdp", "pointmass", "normal"), factor.scales = FALSE,
  px = TRUE, coda = "loadings", coda.scale = TRUE, imh =
  FALSE, imh.iter = 500, imh.burn = 500, ...)
}
\arguments{
  \item{x}{A formula, matrix or bfa object.}

  \item{data}{The data if x is a formula}

  \item{num.factor}{Number of factors}

  \item{restrict}{A matrix or list giving restrictions on
  factor loadings. A matrix should be the same size as the
  loadings matrix. Acceptable values are 0 (identically 0),
  1 (unrestricted), or 2 (strictly positive). List elements
  should be character vectors of the form c('variable',1,
  ">0") where 'variable' is the manifest variable, 1 is the
  factor, and ">0" is the restriction. Acceptable
  restrictions are ">0" or "0".}

  \item{normal.dist}{A character vector specifying which
  variables should be treated as observed Gaussian.
  Defaults to none (a completely semiparametric copula
  model)}

  \item{center.data}{Center continuous variables}

  \item{scale.data}{Scale continuous variables}

  \item{nsim}{Number of iterations past burn-in}

  \item{nburn}{Number of initial (burn-in) iterations to
  discard}

  \item{thin}{Keep every thin'th MCMC sample (i.e. save
  nsim/thin samples)}

  \item{print.status}{How often to print status messages to
  console}

  \item{keep.scores}{Save samples of factor scores}

  \item{keep.loadings}{Save samples of factor loadings}

  \item{loading.prior}{Specify point mass ("pointmass",
  default) or normal priors ("normal")}

  \item{factor.scales}{Include a separate scale parameter
  for each factor}

  \item{px}{Use parameter expansion for
  ordinal/copula/mixed factor models (strongly
  recommended)}

  \item{coda}{Create \code{mcmc} objects to allow use of
  functions from the \code{coda} package: "all" for
  loadings and scores, "loadings" or "scores" for one or
  the other, or "none" for neither}

  \item{coda.scale}{Put the loadings on the correlation
  scale when creating \code{mcmc} objects}

  \item{imh}{Use Independence Metropolis-Hastings step for
  discrete margins. If FALSE, use the semiparametric
  extended rank likelihood}

  \item{imh.iter}{Iterations used to build IMH proposal}

  \item{imh.burn}{Burn-in before collecting samples used to
  build IMH proposal (total burn-in is
  nburn+imh.iter+imh.burn)}

  \item{...}{Prior parameters and other (experimental)
  arguments (see details)}
}
\value{
  An S3 \code{bfa} object \code{model}, with posterior
  samples/summaries.
}
\description{
  This function performs a specified number of MCMC
  iterations and returns an object containing summary
  statistics from the MCMC samples as well as the actual
  samples if keep.scores or keep.loadings are \code{TRUE}.
  Default behavior is to save only the loadings.
}
\details{
  Additional parameters: \itemize{ \item loadings.var:
  Factor loading prior variance \item tau.a, tau.b: Gamma
  hyperparameters (scale=1/b) for factor precisions (if
  factor.scales=T) \item rho.a, rho.b: Beta hyperparameters
  for point mass prior \item sigma2.a, sigma2.b: Gamma
  hyperparameters for error precisions \item gdp.alpha,
  gdp.beta: GDP prior parameters }
}
\examples{
\dontrun{
require(MASS)
data(UScereal)
UScereal$shelf = factor(UScereal$shelf, ordered=TRUE)
UScereal$vitamins = factor(UScereal$vitamins, ordered=TRUE,
levels=c("none", "enriched", "100\%"))
fit_cop = bfa_copula(~., data=UScereal[,-1], num.factor=2, nsim=5000, nburn=500, thin=2,
normal.dist=rep(0,10), rest=list(c("sugars", 2, "0")),
loading.prior="gdp", keep.scores=T, init.fa=FALSE)
plot_loadings(fit_cop)
biplot(fit_cop, cex=c(0.8, 0.8))
plot(get_coda(fit_cop))
plot(get_coda(fit_cop, loadings=F, scores=T))
HPDinterval(fit_cop)
}
}

