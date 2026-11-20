args <- commandArgs(trailingOnly=TRUE)
.libPaths(c(args[2], .libPaths()))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(jsonlite))
if (args[1] == "--probe") {
  cat(toJSON(list(R=R.version.string, edgeR=as.character(packageVersion("edgeR")),
    limma=as.character(packageVersion("limma")), statmod=as.character(packageVersion("statmod")),
    jsonlite=as.character(packageVersion("jsonlite"))), auto_unbox=TRUE))
  quit(status=0)
}
folder <- args[1]
cfg <- fromJSON(file.path(folder, "job.json"))
con <- file(file.path(folder, "counts.bin"), "rb")
counts <- matrix(readBin(con, "double", n=cfg$genes*cfg$n, size=8, endian="little"), nrow=cfg$genes)
close(con)
genes <- read.csv(file.path(folder, "genes.csv"), stringsAsFactors=FALSE)
design <- as.matrix(read.csv(file.path(folder, "design.csv"), check.names=FALSE))
rownames(counts) <- genes$gene_id
stopifnot(ncol(counts)==nrow(design), all(is.finite(counts)), all(counts>=0),
          all(counts==floor(counts)), qr(design)$rank==ncol(design), colnames(design)[2]=="AD")
y <- DGEList(counts=counts)
keep <- filterByExpr(y, design=design)
stopifnot(sum(keep) >= 100, all(colSums(counts)>0))
# Retain whole-transcriptome library totals after the expression filter.
y <- y[keep, , keep.lib.sizes=TRUE]
normalizer_rows <- if (cfg$normalizer=="TMM_nonmarker") !genes$is_any_marker[keep] else rep(TRUE,sum(keep))
stopifnot(sum(normalizer_rows)>=100)
nf <- calcNormFactors(y[normalizer_rows,,keep.lib.sizes=TRUE], method="TMM")$samples$norm.factors
y$samples$norm.factors <- nf
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)
test <- glmQLFTest(fit, coef=2)
tab <- topTags(test, n=Inf, sort.by="none")$table
out <- data.frame(gene_id=rownames(y), log2FC=tab$logFC, logCPM=tab$logCPM,
                  qlf_p=tab$PValue, qlf_q_all_expressed=tab$FDR)
for (dose in c(25,50,100)) {
  lfc <- log2(1+dose/100)
  tr <- glmTreat(fit, coef=2, lfc=lfc)
  tt <- topTags(tr, n=Inf, sort.by="none")$table
  stopifnot(identical(rownames(tt), out$gene_id))
  out[[paste0("treat_p_",dose)]] <- tt$PValue
}
write.csv(out, file.path(folder,"gene_results.csv"),row.names=FALSE)
write.csv(data.frame(gene_id=genes$gene_id, expressed=keep), file.path(folder,"gene_filter.csv"),row.names=FALSE)
write.csv(data.frame(SampleID=cfg$sample_ids, library_size=y$samples$lib.size, TMM=nf,
                    effective_library=y$samples$lib.size*nf), file.path(folder,"normalization_factors.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(folder,"sessionInfo.txt"))
