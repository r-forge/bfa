library(roxygen2)
dir='~/Documents/svn/bfa'
setwd(dir)

f = paste('~/Documents/svn/bfa/pkg/man/',list.files('~/Documents/svn/bfa/pkg/man/'), sep='')
lapply(f, unlink)

roxygenize('pkg', 'pkg', copy.package=FALSE, overwrite=TRUE, 
           unlink.target=TRUE)
