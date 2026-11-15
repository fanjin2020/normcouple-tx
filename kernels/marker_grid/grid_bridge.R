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
marks <- lapply(seq_len(7), function(k) which(ax$marker_owner == k))
readmat <- function(name, nr, nc) {
  data <- readBin(file.path(root, name), 'double', n=nr*nc, size=8, endian='little')
  if (length(data) != nr*nc) stop(paste('Truncated binary', name))
  matrix(data, nrow=nr, ncol=nc)
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
  }, error=function(e) list(p=matrix(NA_real_, n, 7), reason=conditionMessage(e), warnings=paste(unique(warnings), collapse='; ')))
}
rows <- list()
for (entry in plan$entries) {
  n <- as.integer(entry$n)
  raw <- readmat(entry$raw, ax$normalization_genes, n)
  A <- readmat(entry$A, length(ax$marker_rows), 7)
  total0 <- as.numeric(readmat(entry$total, n, 1))
  anchor <- as.numeric(entry$marker_anchor)
  stopifnot(all(raw>=0), all(is.finite(raw)), all(total0>0), all(colSums(A)>0), anchor>0, anchor<=1)
  AM <- sweep(A, 2, colSums(A), '/')
  for (condition in entry$conditions) {
    added <- readmat(condition$additions, length(ax$test_rows), n)
    stress <- readmat(condition$stress, length(ax$marker_rows), n)
    observed <- raw
    observed[ax$test_rows,] <- observed[ax$test_rows,] + added
    observed[ax$marker_rows,] <- observed[ax$marker_rows,] + stress
    total <- total0 + colSums(added) + colSums(stress)
    marker <- t(observed[ax$marker_rows,,drop=FALSE])
    dmarker <- rowSums(marker)
    tmm_error <- ''
    tmm <- tryCatch({
      nf <- edgeR::calcNormFactors(edgeR::DGEList(counts=observed, lib.size=total), method='TMM')
      total * nf$samples$norm.factors
    }, error=function(e) {tmm_error <<- conditionMessage(e); rep(NA_real_, n)})
    specs <- list(observed=list(d=total, A=A, c=1e6), tmm=list(d=tmm, A=A, c=1e6),
      marker_matched=list(d=dmarker, A=AM, c=1e6*anchor),
      marker_denominator_only=list(d=dmarker/anchor, A=A, c=1e6))
    for (method in names(specs)) {
      spec <- specs[[method]]
      value <- if (any(!is.finite(spec$d)) || any(spec$d<=0))
        list(p=matrix(NA_real_, n, 7), reason=paste('invalid_denominator',tmm_error), warnings='') else
        fit(marker/spec$d, spec$A, spec$c, n)
      key <- paste(entry$key, condition$name, sep='_')
      file <- paste0(key, '__', method, '.csv')
      write.table(value$p, file.path(root,file), sep=',', row.names=FALSE, col.names=FALSE, na='nan')
      rows[[length(rows)+1L]] <- data.frame(key=key, method=method, file=file, reason=value$reason,
        warnings=value$warnings, mean_denominator=mean(spec$d), marker_anchor=anchor)
    }
  }
}
write.csv(do.call(rbind, rows), file.path(root,'outputs.csv'), row.names=FALSE)
