args <- commandArgs(TRUE)
.libPaths(c(args[-c(1L,2L)], .libPaths()))
suppressPackageStartupMessages({library(SingleCellExperiment);library(HDF5Array)})
out <- args[2]
sce <- HDF5Array::loadHDF5SummarizedExperiment(args[1])
stopifnot(methods::is(sce,"SingleCellExperiment"), "counts" %in% assayNames(sce))
cd <- as.data.frame(colData(sce))
classes <- c("Astro","EndoMural","Excit","Inhib","Micro","OPC","Oligo")
meta <- data.frame(cell_id=colnames(sce),donor_id=as.character(cd$BrNum),
                   cell_type=as.character(cd$cellType_broad_hc),source_id="LIBD_spatialDLPFC_snRNAseq")
keep <- !is.na(meta$donor_id) & nzchar(meta$donor_id) & meta$cell_type %in% classes
meta <- meta[keep,,drop=FALSE]
sce <- sce[,keep]
valid <- function(x,n) length(x)==n && !anyNA(x) && all(nzchar(trimws(x))) && !anyDuplicated(x)
original <- rownames(sce)
genes <- original
field <- "rownames"
if(!valid(genes,nrow(sce))) {
  values <- lapply(as.data.frame(rowData(sce)),function(x) if(is.atomic(x) && is.null(dim(x))) as.character(x) else NULL)
  candidates <- names(values)[vapply(values,function(x) valid(x,nrow(sce)) && all(grepl("^ENSG[0-9]+(\\.[0-9]+)?$",x)),logical(1))]
  stopifnot(length(candidates)>0)
  canonical <- lapply(values[candidates],function(x) sub("\\.[0-9]+$","",x))
  stopifnot(!any(vapply(canonical,function(x) anyDuplicated(x)>0L,logical(1))),
            all(vapply(canonical,function(x) identical(x,canonical[[1]]),logical(1))))
  field <- candidates[1]
  genes <- values[[field]]
}
cells <- colnames(sce)
stopifnot(length(cells)==nrow(meta), !anyNA(cells), all(nzchar(cells)))
meta$source_cell_id <- meta$cell_id
if(anyDuplicated(cells)) cells <- paste0(nchar(meta$donor_id,type="bytes"),":",meta$donor_id,":",cells)
stopifnot(valid(cells,nrow(meta)))
meta$cell_id <- cells
blocks <- split(seq_len(ncol(sce)),ceiling(seq_len(ncol(sce))/1000))
mats <- lapply(blocks,function(j) as(as.matrix(assay(sce,"counts")[,j,drop=FALSE]),"dgCMatrix"))
x <- do.call(cbind,mats)
dimnames(x) <- list(genes,cells)
stopifnot(all(is.finite(x@x)), all(x@x>=0), all(abs(x@x-round(x@x))<1e-8))
Matrix::writeMM(x,file.path(out,"counts.mtx"))
write.csv(meta,file.path(out,"cells.csv"),row.names=FALSE)
write.csv(data.frame(gene_id=genes),file.path(out,"genes.csv"),row.names=FALSE)
write.csv(data.frame(source_row=seq_len(nrow(sce)),source_gene_id=original,
                    gene_id=genes,identity_field=field),file.path(out,"gene_identity_map.csv"),row.names=FALSE)
cat("LIBD RECONSTRUCTED\n")
