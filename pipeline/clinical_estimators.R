args <- commandArgs(TRUE)
root <- normalizePath(args[1])
.libPaths(c(args[-1L],.libPaths()))
suppressPackageStartupMessages(library(SingleCellExperiment))
script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
source(file.path(dirname(normalizePath(script)),"estimators.R"))
plan <- jsonlite::fromJSON(file.path(root,"plan.json"),simplifyVector=FALSE)
stopifnot(identical(packageDescription("MuSiC")$RemoteSha,plan$music_commit),
          as.character(packageVersion("dtangle"))=="2.0.9",
          as.character(packageVersion("e1071"))=="1.7.17",
          as.character(packageVersion("edgeR"))=="4.8.2")
classes <- unlist(plan$classes)
genes <- unlist(plan$genes)
samples <- unlist(plan$sample_ids)
markers <- unlist(plan$marker_rows)
tests <- unlist(plan$test_rows)
selected <- unlist(plan$stress_rows)
group <- unlist(plan$exposure)
n <- length(samples)
readmat <- function(name,nr,nc) {
  v <- readBin(file.path(root,name),"double",n=nr*nc,size=8,endian="little")
  stopifnot(length(v)==nr*nc)
  matrix(v,nr,nc)
}
raw <- readmat("raw.bin",length(genes),n)
a <- readmat("A.bin",length(markers),7L)
rownames(raw) <- genes
colnames(raw) <- samples
marks <- lapply(seq_len(7L),function(k)which(unlist(plan$marker_owner)==k))
marker_ok <- isTRUE(plan$marker_methods_available)
stopifnot(!any(markers%in%tests),all(raw>=0),all(is.finite(raw)),
          identical(marker_ok,min(lengths(marks))>=10L))
music <- cached_music_engine(plan$reference,classes,plan$music_commit)
yield_marker <- unlist(plan$marker_yield)
yield_full <- unlist(plan$full_yield)
convert <- function(q,yield) {
  p <- sweep(q,2,yield,"/")
  p/rowSums(p)
}
output <- list()
parity <- list()
writep <- function(p,method,probe) {
  p <- as.matrix(p)
  stopifnot(identical(dim(p),c(as.integer(n),7L)),all(is.finite(p)),
            all(p>=0),max(abs(rowSums(p)-1))<1e-6)
  value <- data.frame(sample_id=samples,method=method,probe=probe,p,check.names=FALSE)
  names(value)[-(1:3)] <- classes
  output[[length(output)+1L]] <<- value
}
for(probe in c("baseline","test_counts","marker_counts")) {
  cat("CLINICAL PROBE",probe,"samples",n,"\n"); flush.console()
  observed <- raw
  rows <- if(probe=="test_counts") tests else if(probe=="marker_counts") selected else integer(0)
  if(length(rows)) observed[rows,] <- sweep(observed[rows,,drop=FALSE],2,exp(.1*group),"*")
  obj <- edgeR::calcNormFactors(edgeR::DGEList(counts=observed),method="TMM")
  denominator <- obj$samples$lib.size*obj$samples$norm.factors
  if(probe=="baseline") {
    y <- log2(t(raw[tests,,drop=FALSE])/denominator*1e6+1)
    writeBin(as.double(y),file.path(root,"outcomes.bin"),size=8,endian="little")
    write.table(denominator,file.path(root,"outcome_denominator.csv"),sep=",",row.names=FALSE,col.names=FALSE)
  }
  m <- t(observed[markers,,drop=FALSE])
  if(marker_ok) {
  fit <- dtangle::dtangle(Y=log2(m/denominator*1e6+1),references=log2(t(a)*1e6+1),
                         markers=marks,n_markers=lengths(marks),data_type="rna-seq")$estimates
  writep(convert(fit,yield_full),"dtangle_tmm",probe)
  scale <- 1e6*plan$anchor
  fit <- dtangle::dtangle(Y=log2(m/rowSums(m)*scale+1),
                         references=log2(t(sweep(a,2,colSums(a),"/"))*scale+1),
                         markers=marks,n_markers=lengths(marks),data_type="rna-seq")$estimates
  writep(convert(fit,yield_marker),"dtangle_marker",probe)
  writep(convert(linear_svr_fractions(a,t(m)),yield_marker),"linear_svr_marker",probe)
  }
  policies <- list(music_marker=genes[markers],music_broad=NULL,
                   music_non_test=genes[-tests])
  for(method in names(policies)) {
    if(method=="music_marker" && !marker_ok) next
    cat("CLINICAL METHOD",method,probe,"\n"); flush.console()
    selected_genes <- policies[[method]]
    result <- music(bulk.mtx=observed,sc.sce=NULL,markers=selected_genes,
        clusters="cell_type",samples="donor_id",select.ct=classes,verbose=FALSE,
        iter.max=1000,nu=1e-4,eps=.01,centered=FALSE,normalize=FALSE)
    props <- result$Est.prop.weighted[,classes,drop=FALSE]
    stopifnot(identical(rownames(props),samples))
    writep(props,method,probe)
    if(probe=="baseline") {
      sce <- readRDS(file.path(plan$reference,"reference.rds"))
      original <- MuSiC::music_prop(bulk.mtx=observed[,c(1,n),drop=FALSE],
          sc.sce=sce,markers=selected_genes,clusters="cell_type",samples="donor_id",
          select.ct=classes,verbose=FALSE,iter.max=1000,nu=1e-4,eps=.01,
          centered=FALSE,normalize=FALSE)$Est.prop.weighted[,classes,drop=FALSE]
      error <- max(abs(original-props[c(1,n),,drop=FALSE]))
      stopifnot(is.finite(error),error<1e-10)
      parity[[method]] <- list(samples=2L,maximum_absolute_difference=error)
      rm(sce)
      gc(FALSE)
    }
  }
}
write.csv(do.call(rbind,output),file.path(root,"fractions.csv"),row.names=FALSE)
jsonlite::write_json(list(official_parity=parity,methods=unique(do.call(rbind,output)$method),
    marker_methods_available=marker_ok,marker_counts_by_class=unlist(plan$marker_counts_by_class),
    R=R.version.string,music_commit=plan$music_commit,
    e1071=as.character(packageVersion("e1071")),
    interpretation="Sensitivity on a fixed outcome matrix and gene family; no clinical FDR truth"),
    file.path(root,"engine.json"),auto_unbox=TRUE,pretty=TRUE)
