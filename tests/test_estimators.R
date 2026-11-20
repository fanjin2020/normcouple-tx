args <- commandArgs(TRUE)
.libPaths(c(args,.libPaths()))
script <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
source(file.path(dirname(dirname(normalizePath(script))),'pipeline','estimators.R'))
set.seed(183)
a <- matrix(rexp(280*7,.3),280,7)
for(k in seq_len(7)) a[((k-1)*40+1):(k*40),k] <- a[((k-1)*40+1):(k*40),k]+50
a <- sweep(a,2,colSums(a),'/')
truth <- matrix(rexp(7*8),7,8)
truth <- sweep(truth,2,colSums(truth),'/')
y <- a%*%truth
fit <- linear_svr_fractions(a,y)
stopifnot(all(is.finite(fit)),all(fit>=0),max(abs(rowSums(fit)-1))<1e-12,
          sqrt(mean((fit-t(truth))^2))<.08)
again <- linear_svr_fractions(a,sweep(y,2,seq_len(ncol(y))*100,'*'))
stopifnot(max(abs(fit-again))<1e-8)
bad <- try(linear_svr_fractions(a,y*0),silent=TRUE)
stopifnot(inherits(bad,'try-error'))
cat('LINEAR_SVR_TESTS_PASS rmse=',sqrt(mean((fit-t(truth))^2)), '\n')
