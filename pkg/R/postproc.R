<<<<<<< HEAD
=======
#' @include rng.R
NULL
>>>>>>> stable
#' Invert the matrix lam*t(lam)+diag(u) via Woodbury identity
#' @param lam p by k matrix
#' @param u p dimensional vector
#' @export
woodbury = function(lam, u) {
  k = ncol(lam)
<<<<<<< HEAD
  if(is.null(k)) { k=1; lam=matrix(lam, nrow=length(lam)) }
  Ui = diag(1/u)
  if(k==1) {
    out = Ui - Ui %*% lam  %*% t(lam) %*% Ui/c(1+ t(lam) %*% Ui %*% lam)
  } else {
    out = Ui - Ui %*% lam %*% solve(diag(1, k)+t(lam) %*% Ui %*% lam) %*% t(lam) %*% Ui
  }
  return(out)
}

#' Compute ``regression coefficients" for variable (index) against the rest
#' @param model \code{bfa} model object
#' @param index index of the response variable
=======
  if(is.null(k)) { 
    k=1 
    lam=matrix(lam, nrow=length(u)) 
  }
  Ui = diag(1/u)
  
  M = diag(1, k)+t(lam) %*% Ui %*% lam
  out = Ui - Ui %*% lam %*% solve(M, t(lam)) %*% Ui
  
  return(out)
}

#' Compute regression coefficients for variables in (index) against the rest
#' @param model \code{bfa} model object
#' @param index index(es) of the response 
>>>>>>> stable
#' @export
reg.samp = function(model, index = c(1)) {
  pl = model$post.loadings
  ns = dim(pl)[3]
  out = array(NA, dim=c(length(index),model$P-length(index),ns))
  dimnames(out) = list(model$varlabel[index], model$varlabel[-index], NULL)
  for(i in 1:ns) {
<<<<<<< HEAD
    u = 1/(1 + rowSums(pl[, , i]^2))
    pl[, , i]  = pl[, , i]*sqrt(u)
    mat = t(pl[-index, , i])%*%woodbury(pl[-index, , i], u[-index])
    out[,,i] = pl[index,,i]%*%mat
=======
    if(attr(m, "type") == "gauss") {
      u = model$post.sigma2[i,]
    } else {
      u = 1/(1 + rowSums(pl[, , i]^2))
      pl[, , i]  = pl[, , i]*sqrt(u)
    }
    mat = t(pl[-index, , i])%*%woodbury(pl[-index, , i], u[-index])
    out[, , i] = pl[index, , i]%*%mat
>>>>>>> stable
  }
  return(out)
}

#' Compute correlation matrix at each MCMC iteration
#' @param model \code{bfa} model object
#' @export
corr.samp = function(model) {
  set.seed(1)
  pl=model$post.loadings
  ns = dim(pl)[3]
  out = array(NA, dim=c(model$P, model$P, ns))
  for (i in 1:ns) {
    u = 1/(1 + rowSums(pl[, , i]^2))
    pl[, , i]  = pl[, , i]*sqrt(u) #/sqrt(1 + rowSums(pl[, , i]^2))
    out[, , i] = pl[, , i]%*%t(pl[, , i])+diag(u)
  }
  return(out)
<<<<<<< HEAD
=======
}

#' Posterior predictive and univariate conditional posterior predictive distributions
#' Posterior predictive and univariate conditional posterior predictive distributions, currently
#' implemented only for factor copula models
#' @param object \code{bfa} model object
#' @param resp.var Response variable (see details). If NA, return draws from the joint
#' posterior predictive. Otherwise 
#' @param cond.vars Conditioning variables; either a list like list(X1=val1, X2=val2) with
#' X1, X2 variables in the original data frame, or a P length vector with either the conditioning
#' value or NA (for marginalized variables). Ignored if resp.var is NA
#' @param numeric.as.factor Treat numeric variables as ordinal when conditioning
#' @param ... Ignored
#' @method predict bfa
#' @export

predict.bfa = function(object, resp.var=NA, cond.vars=NA,
                       numeric.as.factor=TRUE, ...) {
  n.mcmc=dim(object$post.loadings)[3]
  
  if(!is.na(resp.var)) {
    y.idx = which(resp.var==colnames(object$original.data))
    
    if(is.list(cond.vars)) {
      x.var = names(cond.vars)
      x.idx = which(colnames(object$original.data) %in% x.var)
    } else {
      x.idx = which(!is.na(cond.vars))
    }
    
    lo = hi = rep(NA, length(x.var))
    ji=1
    for(j in x.idx) {
      if(is.factor(object$original.data[,j]) || numeric.as.factor) {
        cv = as.numeric(factor(cond.vars[[ji]], levels=levels(object$original.data[,j])))
        f = ecdf(object$original.data[,j])
        lo[ji] = qnorm(f(cv-1)); hi[ji] = qnorm(f(cv))
      } else {
        f = ecdf(object$original.data[,j])
        u = f(cond.vars[[ji]])
        if(u==1) u = (object$N-1)/object$N
        lo[ji] = hi[ji] = qnorm(u)
      }
      ji=ji+1
    }
    y.uq = sort(unique(object$original.data[,y.idx]))
    y.samp   = matrix(NA, nrow=n.mcmc, ncol=length(y.uq))
    y.cutpts = qnorm(ecdf(object$original.data[,y.idx])(y.uq))

    for(i in 1:n.mcmc) {
      A = matrix(object$post.loadings[,,i], nrow=object$P)
      u = 1/(1+rowSums(A^2))
      A = A*sqrt(u)
      k = ncol(A)
      eta = rnorm(k)
      xs = matrix(rnorm(length(x.idx)), ncol=1)
      for(h in 1:length(x.idx)) {
        j = x.idx[h]
        xs[h] = rtnorm(lo[h], hi[h], crossprod(A[j,], eta), sqrt(u[j]))
      }
      miv = woodbury(A[x.idx,], u[x.idx])
      rc = tcrossprod(A[y.idx[1],], A[x.idx,])%*%miv
      mu.y = rc%*%xs
      var.y = 1-A[y.idx[1],] %*% t(A[x.idx,]) %*% miv %*% A[x.idx,]%*% matrix(A[y.idx[1],],nrow=k)
      y.samp[i,] = pnorm(y.cutpts, mu.y, sqrt(var.y))
    }
  } else {
    y.samp = matrix(NA, nrow=n.mcmc, ncol=object$P)
    for(i in 1:n.mcmc) {
      A = matrix(object$post.loadings[,,i], nrow=object$P)
      ru = sqrt(1+rowSums(A^2))
      A = A/ru
      k = ncol(A)
      eta = rnorm(k)
      z = A%*%eta+rnorm(object$P, 0, 1/ru)
      y.samp[i,] = z
    }
    y.samp = data.frame(pnorm(y.samp))
    colnames(y.samp) = object$varlabel
    for(j in 1:object$P) {
      y.samp[,j] = quantile(object$ranks[j,]+1, y.samp[,j], type=1)
      y.samp[,j] = object$original.data[match(y.samp[,j], object$ranks[j,]+1), j]
    }
  }
  return(y.samp)
>>>>>>> stable
}