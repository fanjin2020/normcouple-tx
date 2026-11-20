args <- commandArgs(TRUE)
dir.create(args[2], recursive=TRUE, showWarnings=FALSE)
inner <- tempfile(fileext=".rda")
con <- gzfile(args[1], "rb"); to <- file(inner, "wb")
repeat { b <- readBin(con,"raw",1048576); if (!length(b)) break; writeBin(b,to) }
close(con); close(to)
e <- new.env(); load(inner, envir=e)
stopifnot(exists("targets", e), exists("normExpr.reg",e))
x <- as.matrix(e$normExpr.reg); meta <- e$targets
stopifnot(!anyDuplicated(rownames(x)), !anyDuplicated(colnames(x)),
          !anyDuplicated(meta$SampleID), identical(colnames(x), as.character(meta$SampleID)))
write.csv(meta, file.path(args[2],"bulk_metadata.csv"), row.names=FALSE)
write.csv(data.frame(gene_id=rownames(x),x,check.names=FALSE), file.path(args[2],"bulk_processed_expression.csv"),row.names=FALSE)
writeLines(c("Objects: normExpr.reg, targets", "normExpr.reg is processed expression, NOT raw counts.",
  "Never exponentiate/round this matrix to manufacture a raw-count input."),file.path(args[2],"MATRIX_TYPE.txt"))
unlink(inner)
