args <- commandArgs(TRUE)
.libPaths(c(args[2], .libPaths()))
needed <- c('jsonlite', 'dtangle', 'edgeR')
for (pkg in needed) if (!requireNamespace(pkg, quietly=TRUE)) stop(paste('Missing R package:', pkg))
stopifnot(as.character(packageVersion('dtangle')) == '2.0.9')
versions <- setNames(lapply(needed, function(p) as.character(packageVersion(p))), needed)
versions$R <- R.version.string
if (args[1] == '--check') {
  jsonlite::write_json(versions, args[3], auto_unbox=TRUE, pretty=TRUE)
  quit(status=0)
}
root <- args[1]
plan <- jsonlite::fromJSON(file.path(root, 'jobs.json'), simplifyVector=FALSE)
ax <- jsonlite::fromJSON(file.path(plan$reference_folder, 'axes.json'))
stopifnot(length(intersect(ax$test_rows, ax$marker_rows)) == 0)
marks <- lapply(seq_len(7), function(k) which(ax$marker_owner == k))
readmat <- function(name, nr, nc) {
  x <- readBin(file.path(root,name), 'double', n=nr*nc, size=8, endian='little')
  if (length(x) != nr*nc) stop(paste('Truncated binary:', name))
  matrix(x, nrow=nr, ncol=nc)
}
fit <- function(y, a, scale, n) {
  warnings <- character()
  tryCatch({
    if (any(!is.finite(y)) || any(!is.finite(a))) stop('nonfinite_input')
    p <- withCallingHandlers(dtangle::dtangle(Y=log2(y*scale+1), references=log2(t(a)*scale+1),
      markers=marks, n_markers=lengths(marks), data_type='rna-seq')$estimates,
      warning=function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning')})
    if (!identical(dim(p), c(as.integer(n),7L)) || any(!is.finite(p)) || any(p<0) || max(abs(rowSums(p)-1))>1e-5) stop('invalid_prediction')
    list(p=p, reason='', warnings=paste(unique(warnings), collapse='; '))
  }, error=function(e) list(p=matrix(NA_real_,n,7), reason=conditionMessage(e), warnings=paste(unique(warnings),collapse='; ')))
}
predict <- function(marker, spec, n) {
  if (any(!is.finite(spec$d)) || any(spec$d<=0))
    return(list(p=matrix(NA_real_,n,7), reason='invalid_denominator', warnings=''))
  fit(marker/spec$d, spec$A, spec$c, n)
}
rows <- list()
checks <- list()
for (ei in seq_along(plan$entries)) {
  entry <- plan$entries[[ei]]
  n <- as.integer(entry$n)
  raw <- readmat(entry$raw, ax$normalization_genes, n)
  A <- readmat(entry$A, length(ax$marker_rows), 7)
  total0 <- as.numeric(readmat(entry$total, n, 1))
  anchor <- as.numeric(entry$marker_anchor)
  stopifnot(all(raw>=0),all(is.finite(raw)),all(total0>0),all(colSums(A)>0),anchor>0,anchor<=1)
  AM <- sweep(A, 2, colSums(A), '/')
  additions <- lapply(entry$outcomes, function(case) readmat(case$additions, length(ax$test_rows), n))
  for (si in seq_along(entry$subsets)) {
    subset <- entry$subsets[[si]]
    stress <- readmat(subset$stress, length(ax$marker_rows), n)
    stopifnot(all(is.finite(stress)),all(stress>=0))
    stressed <- raw
    stressed[ax$marker_rows,] <- stressed[ax$marker_rows,] + stress
    marker <- t(stressed[ax$marker_rows,,drop=FALSE])
    dmarker <- rowSums(marker)
    specs <- list(marker_matched=list(d=dmarker,A=AM,c=1e6*anchor),
                  marker_denominator_only=list(d=dmarker/anchor,A=A,c=1e6))
    cached <- lapply(specs, function(spec) predict(marker,spec,n))
    files <- setNames(paste0(entry$key,'_',subset$name,'__',names(specs),'.csv'),names(specs))
    for (method in names(cached)) write.table(cached[[method]]$p,file.path(root,files[[method]]),
                        sep=',',row.names=FALSE,col.names=FALSE,na='nan')
    for (ci in seq_along(entry$outcomes)) {
      case <- entry$outcomes[[ci]]
      added <- additions[[ci]]
      stopifnot(all(is.finite(added)),all(added>=0))
      observed <- stressed
      observed[ax$test_rows,] <- observed[ax$test_rows,] + added
      # Cache only after asserting identical marker inputs for each direct effect.
      current_marker <- t(observed[ax$marker_rows,,drop=FALSE])
      stopifnot(identical(current_marker,marker))
      total <- total0+colSums(stress)+colSums(added)
      key <- paste(entry$key,subset$name,case$name,sep='_')
      tmm_error <- ''
      tmm_warnings <- character()
      d <- tryCatch(withCallingHandlers({
        nf <- edgeR::calcNormFactors(edgeR::DGEList(counts=observed,lib.size=total),method='TMM')
        total*nf$samples$norm.factors
      },warning=function(w) {tmm_warnings <<- c(tmm_warnings,conditionMessage(w)); invokeRestart('muffleWarning')}),
        error=function(e) {tmm_error <<- conditionMessage(e); rep(NA_real_,n)})
      value <- predict(marker,list(d=d,A=A,c=1e6),n)
      value$reason <- paste(c(value$reason,tmm_error)[nzchar(c(value$reason,tmm_error))],collapse='; ')
      value$warnings <- paste(c(value$warnings,tmm_warnings)[nzchar(c(value$warnings,tmm_warnings))],collapse='; ')
      file <- paste0(key,'__tmm.csv')
      write.table(value$p,file.path(root,file),sep=',',row.names=FALSE,col.names=FALSE,na='nan')
      rows[[length(rows)+1L]] <- data.frame(key=key,method='tmm',file=file,reason=value$reason,
        warnings=value$warnings,mean_denominator=mean(d),prediction_cached=FALSE)
      for (method in names(cached)) {
        value <- cached[[method]]
        rows[[length(rows)+1L]] <- data.frame(key=key,method=method,file=files[[method]],reason=value$reason,
          warnings=value$warnings,mean_denominator=mean(specs[[method]]$d),prediction_cached=TRUE)
      }
      checked <- ei==1L && si==1L && ci==length(entry$outcomes)
      for (method in names(cached)) {
        finite <- all(is.finite(cached[[method]]$p))
        error <- NA_real_
        if (checked) {
          fresh <- predict(current_marker,specs[[method]],n)
          stopifnot(identical(is.finite(fresh$p),is.finite(cached[[method]]$p)))
          if (finite) {
            error <- max(abs(fresh$p-cached[[method]]$p))
            stopifnot(error < 1e-12)
          }
        }
        checks[[length(checks)+1L]] <- data.frame(key=key,method=method,
          marker_input_identical=TRUE,prediction_finite=finite,fresh_repeat_checked=checked,
          fresh_repeat_max_difference=error)
      }
    }
  }
}
write.csv(do.call(rbind,rows),file.path(root,'outputs.csv'),row.names=FALSE,na='')
write.csv(do.call(rbind,checks),file.path(root,'cache_checks.csv'),row.names=FALSE,na='')
