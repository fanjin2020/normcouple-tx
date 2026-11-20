# Frozen dtangle specifications for finite clinical probes and matched imaging.
run_probes <- function(mode,bank,work) {
stopifnot(mode %in% c("clinical","imaging"))
plan <- jsonlite::fromJSON(file.path(bank,"plan.json"),simplifyVector=FALSE)
vec <- function(name) unname(unlist(plan[[name]]))
readmat <- function(path,nr,nc) {
  stopifnot(file.info(path)$size==8*nr*nc)
  v <- readBin(path,"double",n=nr*nc,size=8,endian="little")
  stopifnot(all(is.finite(v)),all(v>=0))
  matrix(v,nr,nc)
}
classes <- vec("classes"); marks <- vec("marker_rows"); owner <- vec("marker_owner")
a <- readmat(file.path(bank,"A.bin"),length(marks),length(classes))
full_yield <- vec("full_yield")
tmm <- function(raw) {
  obj <- edgeR::calcNormFactors(edgeR::DGEList(counts=raw),method="TMM")
  as.numeric(obj$samples$lib.size*obj$samples$norm.factors)
}
predict_pair <- function(raw,den,samples) {
  if (min(tabulate(owner,nbins=length(classes)))<10) return(NULL)
  byclass <- lapply(seq_along(classes),function(i)which(owner==i))
  m <- t(raw[marks,,drop=FALSE])
  stopifnot(all(rowSums(m)>0),all(colSums(a)>0))
  run <- function(y,r,yield) {
    q <- dtangle::dtangle(Y=y,references=r,markers=byclass,
        n_markers=lengths(byclass),data_type="rna-seq")$estimates
    p <- sweep(q,2,yield,"/"); p <- p/rowSums(p)
    stopifnot(all(is.finite(p)),all(p>=0))
    colnames(p) <- classes
    data.frame(sample_id=samples,p,check.names=FALSE)
  }
  anchor <- 1e6*plan$anchor
  f <- run(log2(m/den*1e6+1),log2(t(a)*1e6+1),full_yield)
  g <- run(log2(m/rowSums(m)*anchor+1),log2(t(sweep(a,2,colSums(a),"/"))*anchor+1),
           if(mode=="clinical") vec("marker_yield") else full_yield*colSums(a))
  rbind(transform(f,method="dtangle_tmm"),transform(g,method="dtangle_marker"))
}
if(mode=="imaging"){
  spec <- jsonlite::fromJSON(file.path(work,"imaging_plan.json"))
  raw <- readmat(file.path(work,"imaging.bin"),spec$genes,length(spec$samples))
  a <- a[spec$panel_indices,,drop=FALSE]; owner <- owner[spec$panel_indices]
  marks <- spec$count_rows
  stopifnot(!anyDuplicated(marks),min(tabulate(owner,nbins=7))>=10)
  p <- predict_pair(raw,tmm(raw),spec$samples)
  write.csv(p,file.path(work,"imaging_fractions.csv"),row.names=FALSE)
}else{
  ids <- vec("sample_ids"); group <- vec("exposure")
  raw <- readmat(file.path(bank,"raw.bin"),length(vec("genes")),length(ids))
  probes <- jsonlite::fromJSON(file.path(work,"probes.json"),simplifyVector=FALSE)
  original <- colSums(raw); d0 <- tmm(raw)
  stored <- scan(file.path(bank,"outcome_denominator.csv"),quiet=TRUE)
  stopifnot(length(d0)==length(stored),max(abs(d0/stored-1))<1e-9)
  out <- list(); denom <- list()
  for (probe in probes){
    rows <- as.integer(unlist(probe$rows)); step <- probe$step
    observed <- raw
    if(length(rows)) observed[rows,] <- sweep(raw[rows,,drop=FALSE],2,exp(step*group),"*")
    d <- if(step==0) d0 else tmm(observed)
    share <- if(length(rows)) colSums(raw[rows,,drop=FALSE])/original else rep(0,length(ids))
    lchange <- colSums(observed)/original-1
    stopifnot(max(abs(lchange-share*expm1(step*group)))<1e-10)
    denom[[probe$id]] <- data.frame(probe=probe$id,sample_id=ids,group=group,
       tested_input_share=share,library_relative_change=lchange,effective_relative_change=d/d0-1)
    p <- predict_pair(observed,d,ids)
    if(!is.null(p)) out[[probe$id]] <- transform(p,probe=probe$id)
  }
  write.csv(do.call(rbind,denom),file.path(work,"dose.csv"),row.names=FALSE)
  if(length(out)) write.csv(do.call(rbind,out),file.path(work,"probe_fractions.csv"),row.names=FALSE)
}
}

if(sys.nframe()==0L) {
  args <- commandArgs(TRUE)
  if(length(args)>3) .libPaths(c(args[-(1:3)],.libPaths()))
  stopifnot(as.character(packageVersion("edgeR"))=="4.8.2",
            as.character(packageVersion("dtangle"))=="2.0.9")
  run_probes(args[1],args[2],args[3])
}
