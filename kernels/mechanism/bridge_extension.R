args <- commandArgs(TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
source(file.path(dirname(normalizePath(script)), "r_shared.R"))
errors <- list(); factors <- list(); controls <- list()
for (entry in plan$entries) {
  n <- as.integer(entry$n)
  raw <- readmat(entry$raw, ax$normalization_genes, n)
  A <- readmat(entry$A, length(ax$marker_rows), 7)
  original <- as.numeric(readmat(entry$total, n, 1))
  stopifnot(all(is.finite(raw)), all(raw >= 0), all(original > 0))
  marker_fraction <- colSums(A)
  if (any(!is.finite(marker_fraction)) || any(marker_fraction <= 0)) stop("INTEGRITY: nonpositive marker yield")
  AM <- sweep(A, 2, marker_fraction, "/")
  anchor <- as.numeric(entry$marker_anchor)
  stopifnot(is.finite(anchor), anchor > 0, anchor <= 1)
  if (!isTRUE(plan$omit_mead)) {
    mead <- safe(function() {
      ref <- getref(unlist(entry$weights))
      bulk <- raw[ax$mead_rows,,drop=FALSE]
      rownames(bulk) <- ref$genes
      colnames(bulk) <- paste0("mix_", seq_len(n))
      t(MEAD::MEAD_est(y=bulk, X=ref$X, Vg=ref$V, w=ref$w,
        calc_var=FALSE, groups=NULL, centeringXY=FALSE, use.QP=FALSE, lambda=1)$p_hat)
    }, n)
    errors[[length(errors)+1L]] <- savep(entry$key, "mead_native", mead)
  }
  stress <- if (is.null(entry$stress)) matrix(0, length(ax$marker_rows), n) else
    readmat(entry$stress, length(ax$marker_rows), n)
  for (cond in entry$conditions) {
    key <- paste(entry$key, cond$name, sep="_")
    added <- readmat(cond$additions, length(ax$test_rows), n)
    observed <- raw
    observed[ax$test_rows,] <- observed[ax$test_rows,] + added
    if (isTRUE(cond$marker_stress)) observed[ax$marker_rows,] <- observed[ax$marker_rows,] + stress
    total <- original + colSums(added) + if (isTRUE(cond$marker_stress)) colSums(stress) else 0
    marker_raw <- t(observed[ax$marker_rows,,drop=FALSE])
    denoms <- list(observed=total, fixed=original, marker=rowSums(marker_raw))
    tmm_error <- ""
    denoms$tmm <- tryCatch({
      obj <- edgeR::calcNormFactors(edgeR::DGEList(counts=observed, lib.size=total), method="TMM")
      total * obj$samples$norm.factors
    }, error=function(e) { tmm_error <<- conditionMessage(e); rep(NA_real_, n) })
    for (kind in names(denoms)) {
      factors[[length(factors)+1L]] <- data.frame(key=key, normalizer=kind, sample=seq_len(n),
        denominator=denoms[[kind]], marker_anchor=anchor)
    }
    specs <- list(observed=list(d=denoms$observed, A=A, c=1e6),
                  fixed=list(d=denoms$fixed, A=A, c=1e6),
                  tmm=list(d=denoms$tmm, A=A, c=1e6),
                  marker_matched=list(d=denoms$marker, A=AM, c=1e6*anchor),
                  marker_unmatched_scale=list(d=denoms$marker, A=AM, c=1e6),
                  marker_denominator_only=list(d=denoms$marker/anchor, A=A, c=1e6))
    if (isTRUE(plan$real_bulk)) specs$fixed <- NULL
    for (method in names(specs)) {
      spec <- specs[[method]]
      value <- safe(function() {
        if (any(!is.finite(spec$d)) || any(spec$d<=0)) stop(paste("invalid_denominator", tmm_error))
        dtfit(marker_raw/spec$d, spec$A, cpm=spec$c, a=1)
      }, n)
      errors[[length(errors)+1L]] <- savep(key, method, value)
    }
  }
}
set.seed(20261013)
Ap <- matrix(rexp(length(ax$marker_rows)*7)+.1, ncol=7)
Yp <- matrix(rexp(length(ax$marker_rows)*12)+.1, nrow=12)
controls$pure_log_scaling_error <- max(abs(dtfit(Yp,Ap,a=0)-dtfit(Yp/exp(seq(-2,2,length.out=12)),Ap,a=0)))
controls$joint_scale_pseudocount_error <- max(abs(dtfit(Yp,Ap,cpm=1e4,a=.01)-dtfit(Yp,Ap,cpm=1e6,a=1)))
stopifnot(controls$pure_log_scaling_error<1e-8, controls$joint_scale_pseudocount_error<1e-8)
write.csv(do.call(rbind,lapply(errors,as.data.frame)),file.path(root,"outputs.csv"),row.names=FALSE)
write.csv(do.call(rbind,factors),file.path(root,"normalization_factors.csv"),row.names=FALSE)
jsonlite::write_json(controls,file.path(root,"r_controls.json"),auto_unbox=TRUE,pretty=TRUE)
