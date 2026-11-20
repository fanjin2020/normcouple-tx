# Linear nu-support-vector regression uses LIBSVM, not a least-squares solver.
# Its nonnegative closure is postprocessing, not an NNLS optimization.
linear_svr_fractions <- function(reference, mixtures, nus=c(.25,.5,.75)) {
  stopifnot(as.character(packageVersion("e1071"))=="1.7.17")
  reference <- as.matrix(reference)
  mixtures <- as.matrix(mixtures)
  stopifnot(nrow(reference)==nrow(mixtures), ncol(reference)==7L,
            all(is.finite(reference)), all(is.finite(mixtures)),
            all(reference>=0), all(mixtures>=0),
            all(colSums(reference)>0), all(colSums(mixtures)>0))
  reference <- sweep(reference,2,colSums(reference),"/")
  mixtures <- sweep(mixtures,2,colSums(mixtures),"/")
  scale_x <- sd(as.vector(reference))
  stopifnot(is.finite(scale_x), scale_x>0)
  x <- (reference-mean(reference))/scale_x
  answer <- matrix(NA_real_,ncol(mixtures),ncol(reference))
  selected <- errors <- numeric(ncol(mixtures))
  for(i in seq_len(ncol(mixtures))) {
    y <- mixtures[,i]
    scale_y <- sd(y)
    stopifnot(is.finite(scale_y),scale_y>0)
    y <- (y-mean(y))/scale_y
    candidates <- lapply(nus,function(nu) {
      fit <- e1071::svm(x=x,y=y,type="nu-regression",kernel="linear",
                        cost=1,nu=nu,scale=FALSE,tolerance=1e-3,cachesize=100)
      weight <- as.numeric(t(fit$coefs)%*%fit$SV)
      weight <- pmax(weight,0)
      if(!all(is.finite(weight)) || sum(weight)<=0)
        return(list(weight=rep(NA_real_,ncol(x)),error=Inf))
      weight <- weight/sum(weight)
      list(weight=weight,error=mean((as.numeric(x%*%weight)-y)^2))
    })
    loss <- vapply(candidates,function(z)z$error,numeric(1))
    best <- which.min(loss)
    stopifnot(is.finite(loss[best]))
    answer[i,] <- candidates[[best]]$weight
    selected[i] <- nus[best]
    errors[i] <- loss[best]
  }
  colnames(answer) <- colnames(reference)
  attr(answer,"nu") <- selected
  attr(answer,"reconstruction_mse") <- errors
  answer
}

cached_music_engine <- function(folder, classes, commit) {
  cached <- readRDS(file.path(folder,"music_basis.rds"))
  stopifnot(identical(cached$commit,commit),identical(cached$classes,classes),
            identical(cached$reference_md5,readLines(file.path(folder,"reference.rds.md5"))))
  basis <- function(x,non.zero=TRUE,markers=NULL,clusters,samples,select.ct=NULL,
                    cell_size=NULL,ct.cov=FALSE,verbose=TRUE) {
    stopifnot(isTRUE(non.zero),identical(clusters,"cell_type"),
              identical(samples,"donor_id"),identical(select.ct,classes),
              is.null(cell_size),!ct.cov)
    value <- cached$basis
    chosen <- if(is.null(markers)) rownames(value$Disgn.mtx) else
      intersect(unlist(markers),rownames(value$Disgn.mtx))
    for(name in c("Disgn.mtx","M.theta","Sigma"))
      value[[name]] <- value[[name]][chosen,,drop=FALSE]
    value
  }
  engine <- MuSiC::music_prop
  environment(engine) <- list2env(list(music_basis=basis),parent=environment(MuSiC::music_prop))
  stopifnot(identical(body(engine),body(MuSiC::music_prop)),
            identical(formals(engine),formals(MuSiC::music_prop)))
  engine
}
