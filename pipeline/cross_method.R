args <- commandArgs(TRUE)
if(args[1]=="--environment") {
  .libPaths(c(args[-1L], .libPaths()))
  packages <- c("MuSiC","TOAST","nnls","SingleCellExperiment","SummarizedExperiment",
                "Matrix","dtangle","edgeR","jsonlite")
  versions <- setNames(lapply(packages,function(p) as.character(packageVersion(p))),packages)
  stopifnot(identical(packageDescription("MuSiC")$RemoteSha,"f21fe67f5670d5e9fca0ad7550abaae3423eb59c"),
            versions$dtangle=="2.0.9",versions$edgeR=="4.8.2")
  cat(jsonlite::toJSON(list(R=R.version.string,platform=R.version$platform,packages=versions,
                          music_commit=packageDescription("MuSiC")$RemoteSha),auto_unbox=TRUE))
  quit(status=0)
}
if(args[1]=="--prepare") {
  .libPaths(c(args[-c(1L,2L)], .libPaths()))
  suppressPackageStartupMessages(library(SingleCellExperiment))
  folder <- args[2]
  destination <- file.path(folder,"reference.rds")
  checksum <- paste0(destination,".md5")
  input_file <- paste0(destination,".inputs.md5")
  inputs <- unname(tools::md5sum(file.path(folder,c("reference.mtx","genes.csv","cells.csv"))))
  if(file.exists(destination) && file.exists(checksum) && file.exists(input_file)) {
    stopifnot(identical(unname(tools::md5sum(destination)), readLines(checksum)),
              identical(inputs,readLines(input_file)))
  } else {
    x <- as(Matrix::readMM(file.path(folder,"reference.mtx")), "CsparseMatrix")
    rownames(x) <- read.csv(file.path(folder,"genes.csv"))$gene_id
    meta <- read.csv(file.path(folder,"cells.csv"))
    stopifnot(ncol(x)==nrow(meta), all(Matrix::colSums(x)>0))
    sce <- SingleCellExperiment(assays=list(counts=x), colData=meta)
    temporary <- paste0(destination,".partial")
    saveRDS(sce,temporary,compress=FALSE)
    stopifnot(file.rename(temporary,destination))
    writeLines(unname(tools::md5sum(destination)),checksum)
    writeLines(inputs,input_file)
  }
  basis_file <- file.path(folder,"music_basis.rds")
  if(file.exists(basis_file) && file.exists(paste0(basis_file,".md5"))) {
    stopifnot(identical(unname(tools::md5sum(basis_file)), readLines(paste0(basis_file,".md5"))))
    value <- readRDS(basis_file)
    stopifnot(identical(value$reference_md5,readLines(checksum)),
              identical(value$commit,packageDescription("MuSiC")$RemoteSha))
  } else {
    sce <- readRDS(destination)
    classes <- c("Astro","EndoMural","Excit","Inhib","Micro","OPC","Oligo")
    basis <- MuSiC::music_basis(sce, non.zero=TRUE, markers=NULL, clusters="cell_type", samples="donor_id",
                              select.ct=classes, verbose=FALSE)
    value <- list(basis=basis, classes=classes, reference_md5=readLines(checksum),
                  commit=packageDescription("MuSiC")$RemoteSha)
    temporary <- paste0(basis_file,".partial")
    saveRDS(value,temporary,compress=FALSE)
    stopifnot(file.rename(temporary,basis_file))
    writeLines(unname(tools::md5sum(basis_file)),paste0(basis_file,".md5"))
  }
  cat("REFERENCE READY",folder,"\n")
  quit(status=0)
}
root <- normalizePath(args[1])
.libPaths(c(args[-1L], .libPaths()))
plan <- jsonlite::fromJSON(file.path(root, "plan.json"), simplifyVector=FALSE)
.libPaths(c(unlist(plan$libraries), .libPaths()))
suppressPackageStartupMessages(library(SingleCellExperiment))
for (p in c("MuSiC", "dtangle", "edgeR", "jsonlite", "Matrix")) {
  if (!requireNamespace(p, quietly=TRUE)) stop(paste("Missing dependency", p))
}
description <- packageDescription("MuSiC")
stopifnot(identical(description$RemoteSha, plan$music_commit),
          as.character(packageVersion("dtangle")) == "2.0.9",
          as.character(packageVersion("edgeR")) == "4.8.2")
classes <- unlist(plan$classes)
genes <- unlist(plan$norm_genes)
markers <- unlist(plan$marker_rows)
tests <- unlist(plan$test_rows)
owner <- unlist(plan$marker_owner)
stopifnot(!any(markers %in% tests), length(tests)==100L)
cached <- readRDS(file.path(plan$reference,"music_basis.rds"))
stopifnot(identical(cached$commit, plan$music_commit), identical(cached$classes, classes),
          identical(cached$reference_md5,readLines(file.path(plan$reference,"reference.rds.md5"))))
# Only reference moments are memoized; the official fitting function body and
# its zero filtering, normalization, weighting, and solver remain unchanged.
cached_basis <- function(x, non.zero=TRUE, markers=NULL, clusters, samples, select.ct=NULL,
                         cell_size=NULL, ct.cov=FALSE, verbose=TRUE) {
  stopifnot(isTRUE(non.zero), !is.null(markers), identical(clusters,"cell_type"),
            identical(samples,"donor_id"), identical(select.ct,classes), is.null(cell_size), !ct.cov)
  value <- cached$basis
  chosen <- intersect(unlist(markers), rownames(value$Disgn.mtx))
  for(name in c("Disgn.mtx","M.theta","Sigma")) value[[name]] <- value[[name]][chosen,,drop=FALSE]
  value
}
music_cached <- MuSiC::music_prop
environment(music_cached) <- list2env(list(music_basis=cached_basis), parent=environment(MuSiC::music_prop))
stopifnot(identical(body(music_cached),body(MuSiC::music_prop)), identical(formals(music_cached),formals(MuSiC::music_prop)))
marks <- lapply(seq_len(7L), function(k) which(owner==k))
readmat <- function(name, nr, nc) {
  value <- readBin(file.path(root, name), "double", n=nr*nc, size=8, endian="little")
  stopifnot(length(value)==nr*nc)
  matrix(value, nr, nc)
}
check <- function(x, n) {
  x <- as.matrix(x)
  stopifnot(identical(dim(x), c(as.integer(n),7L)), all(is.finite(x)), all(x>=0), max(abs(rowSums(x)-1))<1e-6)
  x
}
savep <- function(value, key, probe, method) {
  write.table(value, file.path(root, paste(key, probe, paste0(method,".csv"), sep="_")),
              sep=",", col.names=FALSE, row.names=FALSE)
}
dtfit <- function(y, a, scale=1e6) {
  check(dtangle::dtangle(Y=log2(y*scale+1), references=log2(t(a)*scale+1),
                        markers=marks, n_markers=lengths(marks), data_type="rna-seq")$estimates, nrow(y))
}
bulks <- list()
index <- list()
diagnostics <- list()
for (entry in plan$entries) {
  key <- entry$key
  n <- as.integer(entry$n)
  raw <- readmat(entry$raw, length(genes), n)
  total <- as.numeric(readmat(paste0(key,"_total.bin"), n, 1L))
  a <- readmat(paste0(key,"_A.bin"), length(markers), 7L)
  added <- readmat(paste0(key,"_added.bin"), length(tests), n)
  extra <- readmat(paste0(key,"_extra.bin"), length(tests), n)
  stress <- readmat(paste0(key,"_stress.bin"), length(markers), n)
  for (probe in c("baseline", "test_counts", "marker_counts")) {
    observed <- raw
    observed[tests,] <- observed[tests,] + added + if(probe=="test_counts") extra else 0
    if(probe=="marker_counts") observed[markers,] <- observed[markers,] + stress
    library <- total + colSums(added) + if(probe=="test_counts") colSums(extra) else if(probe=="marker_counts") colSums(stress) else 0
    stopifnot(all(is.finite(observed)), all(observed>=0), all(library>0), all(colSums(observed[markers,,drop=FALSE])>0))
    obj <- edgeR::calcNormFactors(edgeR::DGEList(counts=observed, lib.size=library), method="TMM")
    denom <- library * obj$samples$norm.factors
    m <- t(observed[markers,,drop=FALSE])
    savep(dtfit(m/denom, a), key, probe, "dtangle_tmm")
    savep(dtfit(m/rowSums(m), sweep(a,2,colSums(a),"/"), 1e6*entry$anchor), key, probe, "dtangle_marker")
    rownames(observed) <- genes
    colnames(observed) <- paste(key, probe, seq_len(n), sep="_")
    bulks[[length(bulks)+1L]] <- observed
    index[[length(index)+1L]] <- list(key=key, probe=probe, n=n)
    diagnostics[[length(diagnostics)+1L]] <- data.frame(key=key, probe=probe,
       total_min=min(library), total_max=max(library), marker_min=min(rowSums(m)),
       all_zero_markers=sum(colSums(m)==0))
  }
}
bulk <- do.call(cbind, bulks)
result <- music_cached(bulk.mtx=bulk, sc.sce=NULL, markers=genes[markers],
    clusters="cell_type", samples="donor_id", select.ct=classes, verbose=FALSE,
    iter.max=1000, nu=1e-4, eps=.01, centered=FALSE, normalize=FALSE)
p <- result$Est.prop.weighted
stopifnot(all(classes %in% colnames(p)), identical(rownames(p), colnames(bulk)))
p <- check(p[,classes,drop=FALSE], ncol(bulk))
parity <- list(checked=FALSE)
if(isTRUE(plan$validate_official)) {
  n <- as.integer(plan$entries[[1]]$n)
  selected <- unlist(lapply(seq_along(index), function(i) (i-1L)*n+c(1L,n)))
  sce <- readRDS(file.path(plan$reference,"reference.rds"))
  original <- MuSiC::music_prop(bulk.mtx=bulk[,selected,drop=FALSE], sc.sce=sce, markers=genes[markers],
      clusters="cell_type", samples="donor_id", select.ct=classes, verbose=FALSE,
      iter.max=1000, nu=1e-4, eps=.01, centered=FALSE, normalize=FALSE)$Est.prop.weighted[,classes,drop=FALSE]
  error <- max(abs(original-p[selected,,drop=FALSE]))
  stopifnot(is.finite(error), error<1e-10)
  parity <- list(checked=TRUE, samples=length(selected), max_error=error,
                 comparison="Uncached official API on a separately fitted subset of samples")
}
offset <- 0L
for (entry in index) {
  rows <- offset + seq_len(entry$n)
  savep(p[rows,,drop=FALSE], entry$key, entry$probe, "music_marker")
  offset <- offset + entry$n
}
versions <- setNames(lapply(c("MuSiC","dtangle","edgeR","SingleCellExperiment","Matrix","jsonlite"),
                           function(x) as.character(packageVersion(x))),
                     c("MuSiC","dtangle","edgeR","SingleCellExperiment","Matrix","jsonlite"))
jsonlite::write_json(list(R=R.version.string, versions=versions, music_commit=description$RemoteSha, official_parity=parity,
                         samples=ncol(bulk), genes=nrow(bulk), markers=length(markers),
                         diagnostics=do.call(rbind, diagnostics)), file.path(root,"engine.json"),
                     auto_unbox=TRUE, pretty=TRUE)
