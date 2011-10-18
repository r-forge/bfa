#' Invert the matrix lam*t(lam)+diag(u) via Woodbury identity
#' @param lam p by k matrix
#' @param u p dimensional vector
#' @export
woodbury = function(lam, u) {
  k = ncol(lam)
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
#' @export
reg.samp = function(model, index = c(1)) {
  pl = model$post.loadings
  ns = dim(pl)[3]
  out = array(NA, dim=c(length(index),model$P-length(index),ns))
  dimnames(out) = list(model$varlabel[index], model$varlabel[-index], NULL)
  for(i in 1:ns) {
    u = 1/(1 + rowSums(pl[, , i]^2))
    pl[, , i]  = pl[, , i]*sqrt(u)
    mat = t(pl[-index, , i])%*%woodbury(pl[-index, , i], u[-index])
    out[,,i] = pl[index,,i]%*%mat
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
}