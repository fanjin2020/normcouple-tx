args <- commandArgs(TRUE)
.libPaths(c(args[1],.libPaths()))
suppressPackageStartupMessages(library(SingleCellExperiment))
set.seed(1012)
K <- 7L; G <- 350L; D <- 4L; cells <- 10L
classes <- c("Astro","EndoMural","Excit","Inhib","Micro","OPC","Oligo")
genes <- paste0("g",seq_len(G))
celltype <- rep(rep(classes,each=cells),D)
individual <- rep(paste0("d",seq_len(D)),each=K*cells)
mu <- matrix(1,G,K)
for (k in seq_len(K)) mu[((k-1)*50+1):(k*50),k] <- 20
x <- matrix(0,G,length(celltype))
for (i in seq_len(ncol(x))) x[,i] <- rpois(G,mu[,match(celltype[i],classes)]*runif(G,.85,1.15))
rownames(x) <- genes
ref <- SingleCellExperiment(assays=list(counts=Matrix::Matrix(x,sparse=TRUE)))
ref$cell_type <- celltype
ref$individual <- individual
bulk <- matrix(rpois(G*12,5),G,12,dimnames=list(genes,paste0("mix",seq_len(12))))
for (i in seq_len(12)) bulk[,i] <- bulk[,i]+round(mu%*%runif(K,1,15))
R01 <- Matrix::Diagonal(G)
dimnames(R01) <- list(genes,genes)
set.seed(111)
direct <- MEAD::MEAD(bulk,ref,cell_types=classes,
  preprocessing_control=list(filter.gene=FALSE),
  estimation_control=list(calc_var=FALSE,groups=NULL),R01=R01)$p_hat
set.seed(111)
mom <- MEAD::MEAD_getX(ref,classes,unique(individual))
modular <- MEAD::MEAD_est(bulk,mom$X,mom$V,mom$w,calc_var=FALSE,groups=NULL)$p_hat
delta <- max(abs(direct-modular))
stopifnot(is.finite(delta),delta<1e-10,identical(colnames(mom$X),classes),
          all(is.finite(modular)),max(abs(colSums(modular)-1))<1e-10)
lambda0 <- MEAD::MEAD_est(bulk,mom$X,mom$V,mom$w,calc_var=FALSE,lambda=0)$p_hat
w <- mom$w/sum(mom$w)*G
raw <- pmax(solve(crossprod(mom$X*sqrt(w)),crossprod(mom$X*sqrt(w),bulk*sqrt(w))),0)
expected <- sweep(raw,2,colSums(raw),"/")
lambda0_error <- max(abs(expected-lambda0))
stopifnot(lambda0_error<1e-9)
total <- colSums(bulk)
tmm <- edgeR::calcNormFactors(edgeR::DGEList(counts=bulk,lib.size=total),method="TMM")$samples$norm.factors
sf <- DESeq2::estimateSizeFactorsForMatrix(bulk,type="poscounts")
stopifnot(all(is.finite(tmm)),all(tmm>0),abs(exp(mean(log(tmm)))-1)<1e-8,
          all(is.finite(sf)),all(sf>0))
jsonlite::write_json(list(status="passed",MEAD_top_level_modular_max_error=delta,
   MEAD_lambda0_weighted_fit_max_error=lambda0_error,TMM_factors_positive=TRUE,
   DESeq2_poscounts_factors_positive=TRUE,fixture_only=TRUE),args[2],auto_unbox=TRUE,pretty=TRUE)
