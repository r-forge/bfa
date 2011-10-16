\name{fit_bfa}
\alias{fit_bfa}
\title{Perform MCMC model fitting for a bfa model...}
\usage{fit_bfa(model, nsim, nburn, thin=1, print.status=500, keep.scores=FALSE, keep.loadings=TRUE,
    loading.prior="pointmass", factor.scales=FALSE, px=TRUE, coda="loadings",
    coda.scale=TRUE, ...)
}
\description{Perform MCMC model fitting for a bfa model}
\details{This function performs a specified number of MCMC iterations and
returns an object containing summary statistics from the MCMC samples
as well as the actual samples if keep.scores or keep.loadings are \code{TRUE}.
Default behavior is to save only the loadings. 

Prior parameters:
loadings.var: Factor loading prior variance
tau.a, tau.b: Gamma hyperparameters (scale=1/b) for factor precisions (if factor.scales=T)
rho.a, rho.b: Beta hyperparameters for point mass prior
sigma2.a, sigma2.b: Gamma hyperparameters for error precisions
gdp.alpha, gdp.beta: GDP prior parameters}
\value{The S3 \code{bfa} object \code{model}, now with posterior samples/summaries.}
\arguments{\item{model}{an object of type bfa, as returned by bfa(data)}
\item{nsim}{number of iterations past burn-in}
\item{nburn}{number of initial (burn-in) iterations to discard}
\item{thin}{keep every thin'th MCMC sample (i.e. save nsim/thin samples)}
\item{print.status}{how often to print status messages to console}
\item{keep.scores}{save samples of factor scores}
\item{keep.loadings}{save samples of factor loadings}
\item{loading.prior}{Specify point mass ("pointmass", default) or normal priors ("normal")}
\item{factor.scales}{Include a separate scale parameter for each factor}
\item{px}{Use parameter expansion for ordinal/copula/mixed factor models (recommended)}
\item{coda}{create \code{mcmc} objects to allow use of functions from the 
\code{coda} package: "all" for loadings and scores, "loadings" or "scores" for one or the
other, or "none" for neither}
\item{coda.scale}{put the loadings on the correlation scale when creating \code{mcmc} objects}
\item{...}{Prior parameters and other (experimental) arguments (see details)}
}

