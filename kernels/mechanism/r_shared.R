args <- commandArgs(trailingOnly=TRUE)
.libPaths(c(args[2], .libPaths()))
needed <- c("jsonlite", "Matrix", "dtangle", "edgeR", "DESeq2", "SingleCellExperiment", "MEAD", "vashr")
for (p in needed) if (!requireNamespace(p, quietly=TRUE)) stop(paste("Missing R dependency:", p))
version_of <- function(p) as.character(getNamespaceVersion(p))
stopifnot(version_of("dtangle") == "2.0.9", version_of("MEAD") == "1.0.2")
suppressPackageStartupMessages(library(SingleCellExperiment))
versions <- setNames(lapply(needed, version_of), needed)
versions$R <- R.version.string
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
upstream <- new.env(parent=asNamespace("MEAD"))
for (file in list.files(file.path(dirname(normalizePath(script)),"vendor","MEAD","R"),pattern="[.]R$",full.names=TRUE)) {
  sys.source(file,envir=upstream)
}
for (fun in c("MEAD_getX","MEAD_est","MEAD_preprocessing")) {
  installed <- getExportedValue("MEAD",fun)
  stopifnot(identical(deparse(body(installed)),deparse(body(upstream[[fun]]))),
            identical(formals(installed),formals(upstream[[fun]])))
}
if (args[1] == "--check") {
  jsonlite::write_json(versions, args[3], auto_unbox=TRUE, pretty=TRUE)
  quit(status=0)
}
root <- args[1]
plan <- jsonlite::fromJSON(file.path(root, "jobs.json"), simplifyVector=FALSE)
ax <- jsonlite::fromJSON(file.path(plan$reference_folder, "axes.json"))
marks <- lapply(seq_len(7), function(k) which(ax$marker_owner == k))
readmat <- function(name, nr, nc) {
  x <- readBin(file.path(root, name), "double", n=nr*nc, size=8, endian="little")
  if (length(x) != nr*nc) stop(paste("Truncated binary", name))
  matrix(x, nrow=nr, ncol=nc)
}
checkp <- function(x, n) {
  x <- as.matrix(x)
  stopifnot(identical(dim(x), c(as.integer(n),7L)), all(is.finite(x)), all(x>=0), max(abs(rowSums(x)-1))<1e-5)
  x
}
dtfit <- function(Y, A, cpm=1e6, a=1) {
  if (a == 0 && (any(Y<=0) || any(A<=0))) stop("pure_log_requires_strictly_positive_input")
  checkp(dtangle::dtangle(Y=log2(Y*cpm+a), references=log2(t(A)*cpm+a),
     markers=marks, n_markers=lengths(marks), data_type="rna-seq")$estimates, nrow(Y))
}
safe <- function(fun, n) {
  warnings <- character()
  tryCatch({
    value <- withCallingHandlers(fun(), warning=function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    })
    list(p=checkp(value,n), reason="", warnings=paste(unique(warnings),collapse="; "))
  }, error=function(e) {
    if (grepl("^(INTEGRITY|INFRASTRUCTURE):",conditionMessage(e))) stop(e)
    list(p=matrix(NA_real_,n,7), reason=conditionMessage(e), warnings=paste(unique(warnings),collapse="; "))
  })
}
savep <- function(key, arm, result) {
  file <- paste0(key,"__",arm,".csv")
  write.table(result$p, file.path(root,file), sep=",",row.names=FALSE,col.names=FALSE,na="nan")
  list(key=key, arm=arm, file=file, reason=result$reason, warnings=result$warnings)
}

# These calls are upstream MEAD functions, without edits to the estimator.
# Repeated reference draws retain distinct draw IDs, not fictitious new donors.
ref_raw <- NULL
getref <- function(weights) {
  code <- paste(weights,collapse="_")
  dir.create(plan$mead_cache, recursive=TRUE, showWarnings=FALSE)
  dest <- file.path(plan$mead_cache,paste0(code,".rds"))
  receipt <- paste0(dest,".md5")
  if (file.exists(dest) && file.exists(receipt)) {
    if (!identical(unname(tools::md5sum(dest)),readLines(receipt,warn=FALSE))) stop("INTEGRITY: damaged MEAD reference cache")
    return(readRDS(dest))
  }
  if (is.null(ref_raw)) {
    x <- as(Matrix::readMM(file.path(plan$reference_folder,"reference.mtx")),"CsparseMatrix")
    meta <- read.csv(file.path(plan$reference_folder,"reference_cells.csv"),stringsAsFactors=FALSE)
    genes <- read.csv(file.path(plan$reference_folder,"mead_genes.csv"),stringsAsFactors=FALSE)$gene_id
    stopifnot(ncol(x)==nrow(meta),nrow(x)==length(genes))
    rownames(x) <- genes
    ref_raw <<- list(x=x,meta=meta,genes=genes)
  }
  draws <- rep(seq_along(weights),weights)
  ix <- unlist(lapply(draws,function(j) which(ref_raw$meta$ref_index==j-1L)))
  ids <- unlist(lapply(seq_along(draws),function(j) rep(paste0("draw_",j),sum(ref_raw$meta$ref_index==draws[j]-1L))))
  x <- ref_raw$x[,ix,drop=FALSE]
  keep <- Matrix::colSums(x)>0
  sce <- SingleCellExperiment(assays=list(counts=x[,keep,drop=FALSE]))
  sce$individual <- ids[keep]
  sce$cell_type <- ref_raw$meta$cell_type[ix][keep]
  tab <- table(sce$individual,factor(sce$cell_type,levels=ax$classes))
  if (!all(tab>0)) stop("MEAD_reference_class_missing_after_zero_library_filter")
  set.seed(1012L+sum(weights*c(1L,10L,100L,1000L)))
  value <- MEAD::MEAD_getX(sce,cell_types=ax$classes,individuals=paste0("draw_",seq_along(draws)))
  value$genes <- rownames(sce)
  stopifnot(identical(colnames(value$X),ax$classes))
  tmp <- paste0(dest,".",Sys.getpid(),".partial")
  saveRDS(value,tmp)
  if (!file.rename(tmp,dest)) stop("INFRASTRUCTURE: MEAD reference cache write failed")
  writeLines(unname(tools::md5sum(dest)),receipt)
  value
}

