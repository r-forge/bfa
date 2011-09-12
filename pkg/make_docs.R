library(roxygen)
setwd('~/Dropbox/')
f = paste('~/Dropbox/bfa/man/',list.files('~/Dropbox/bfa/man/'), sep='')
lapply(f, unlink)

roxygenize('bfa', 'bfa', copy.package=FALSE, overwrite=TRUE, 
           unlink.target=FALSE, use.Rd2=TRUE)
