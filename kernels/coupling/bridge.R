args <- commandArgs(TRUE)
.libPaths(c(args[2], .libPaths()))
needed <- c('jsonlite','dtangle','edgeR')
for (pkg in needed) if (!requireNamespace(pkg,quietly=TRUE)) stop(paste('Missing R package:',pkg))
stopifnot(as.character(packageVersion('dtangle')) == '2.0.9')
if (args[1] == '--check') {
  v <- setNames(lapply(needed,function(p) as.character(packageVersion(p))),needed)
  v$R <- R.version.string
  jsonlite::write_json(v,args[3],auto_unbox=TRUE,pretty=TRUE)
  quit(status=0)
}
root <- args[1]
plan <- jsonlite::fromJSON(file.path(root,'jobs.json'),simplifyVector=FALSE)
ax <- plan$axes
marks <- lapply(seq_len(7),function(k) which(unlist(ax$marker_owner)==k))
mr <- as.integer(unlist(ax$marker_rows))
tr <- as.integer(unlist(ax$test_rows))
ng <- as.integer(ax$normalization_genes)
stopifnot(!anyDuplicated(mr),!anyDuplicated(tr),length(intersect(mr,tr))==0,
          all(mr>=1),all(mr<=ng),all(tr>=1),all(tr<=ng),all(lengths(marks)>0))
non_test <- setdiff(seq_len(ng),tr)
readmat <- function(name,nr,nc) {
  v <- readBin(file.path(root,name),'double',n=nr*nc,size=8,endian='little')
  stopifnot(length(v)==nr*nc)
  matrix(v,nrow=nr,ncol=nc)
}
fit <- function(marker,d,A,scale,n) {
  warnings <- character()
  ans <- tryCatch(withCallingHandlers({
    if(any(!is.finite(d))||any(d<=0)) stop('invalid_denominator')
    q <- dtangle::dtangle(Y=log2(marker/d*scale+1), references=log2(t(A)*scale+1),
      markers=marks,n_markers=lengths(marks),data_type='rna-seq')$estimates
    if(!identical(dim(q),c(as.integer(n),7L))||any(!is.finite(q))||any(q<0)||
       max(abs(rowSums(q)-1))>1e-5) stop('invalid_prediction')
    list(p=q,reason='')
  },warning=function(w) {warnings <<- c(warnings,conditionMessage(w)); invokeRestart('muffleWarning')}),
  error=function(e) list(p=matrix(NA_real_,n,7),reason=conditionMessage(e)))
  ans$warnings <- paste(unique(warnings),collapse='; ')
  ans
}
effective <- function(counts,lib) {
  if(any(!is.finite(lib))||any(lib<=0)) stop('nonpositive_TMM_library')
  z <- edgeR::calcNormFactors(edgeR::DGEList(counts=counts,lib.size=lib),method='TMM')
  lib*z$samples$norm.factors
}
rows <- list()
for (entry in plan$entries) {
  n <- as.integer(entry$n)
  raw <- readmat(entry$raw,ng,n)
  A <- readmat(entry$A,length(mr),7)
  total <- as.numeric(readmat(entry$total,n,1))
  group <- as.numeric(unlist(entry$group))
  anchor <- as.numeric(entry$marker_anchor)
  non_anchor <- as.numeric(entry$non_test_anchor)
  stopifnot(all(raw>=0),all(is.finite(raw)),length(group)==n,
            anchor>0,anchor<=1,non_anchor>0,non_anchor<=1,all(colSums(A)>0))
  non_total <- total-colSums(raw[tr,,drop=FALSE])
  stopifnot(all(non_total>0))
  marker <- t(raw[mr,,drop=FALSE])
  AM <- sweep(A,2,colSums(A),'/')
  # Recompute the complete-exclusion factor from an unchanged count/library universe.
  excluded_error <- ''
  excluded <- tryCatch(effective(raw[non_test,,drop=FALSE],non_total)/non_anchor,
                       error=function(e) {excluded_error <<- conditionMessage(e); rep(NA_real_,n)})
  for (step in unlist(entry$steps)) {
    counts <- raw
    counts[tr,] <- sweep(raw[tr,,drop=FALSE],2,exp(step*group),'*')
    library <- non_total+colSums(counts[tr,,drop=FALSE])
    stopifnot(identical(t(counts[mr,,drop=FALSE]),marker))
    specifications <- list(
      tmm=list(A=A,c=1e6),
      tmm_exclude_factor=list(A=A,c=1e6),
      tmm_exclude_both=list(A=A,c=1e6),
      marker_denominator_only=list(A=A,c=1e6),
      marker_matched=list(A=AM,c=1e6*anchor))
    for (method in names(specifications)) {
      warning_text <- character()
      error_text <- ''
      d <- tryCatch(withCallingHandlers(switch(method,
        tmm=effective(counts,library),
        tmm_exclude_factor=effective(counts[non_test,,drop=FALSE],library),
        tmm_exclude_both=excluded,
        marker_denominator_only=rowSums(marker)/anchor,
        marker_matched=rowSums(marker)),
        warning=function(w) {warning_text <<- c(warning_text,conditionMessage(w)); invokeRestart('muffleWarning')}),
        error=function(e) {error_text <<- conditionMessage(e); rep(NA_real_,n)})
      spec <- specifications[[method]]
      if(method=='tmm_exclude_both') error_text <- excluded_error
      result <- fit(marker,d,spec$A,spec$c,n)
      file <- sprintf('%s_%s_%+.8f.csv',entry$key,method,step)
      write.table(cbind(d,result$p),file.path(root,file),sep=',',row.names=FALSE,col.names=FALSE,na='nan')
      rows[[length(rows)+1L]] <- data.frame(key=entry$key,method=method,step=step,file=file,
        reason=paste(c(error_text,result$reason)[nzchar(c(error_text,result$reason))],collapse='; '),
        warnings=paste(c(warning_text,result$warnings)[nzchar(c(warning_text,result$warnings))],collapse='; '),
        scale=spec$c,marker_input_unchanged=TRUE)
    }
  }
}
write.csv(do.call(rbind,rows),file.path(root,'outputs.csv'),row.names=FALSE,na='')
