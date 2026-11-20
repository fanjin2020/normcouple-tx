# Deterministic scientific-kernel check with a synthetic bank, never a clinical result.
args <- commandArgs(TRUE)
if(length(args)) .libPaths(c(args,.libPaths()))
this <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
source(file.path(dirname(dirname(normalizePath(this))),"pipeline","probe_engine.R"))
set.seed(125583)
work <- tempfile("normcouple-fixture-"); dir.create(work)
n <- 20L; g <- 110L; k <- 7L
owners <- rep(seq_len(k),each=10)
marks <- seq_along(owners)
ids <- paste0("S",seq_len(n)); genes <- paste0("G",seq_len(g))
classes <- c("Astro","EndoMural","Excit","Inhib","Micro","Oligo","OPC")
a <- matrix(runif(length(marks)*k,.0001,.002),length(marks),k)
a[cbind(marks,owners)] <- a[cbind(marks,owners)]+.03
mix <- matrix(rexp(n*k),k,n); mix <- sweep(mix,2,colSums(mix),"/")
raw <- rbind(round(1e5*a%*%mix)+1,matrix(rpois((g-length(marks))*n,250)+1,g-length(marks),n))
x <- rep(0:1,length.out=n)
fy <- seq(.8,1.4,length.out=k); my <- fy*colSums(a)
plan <- list(genes=genes,sample_ids=ids,classes=classes,marker_rows=marks,marker_owner=owners,
   exposure=x,anchor=.308,full_yield=fy,marker_yield=my)
savej <- function(obj,name) jsonlite::write_json(obj,file.path(work,name),auto_unbox=TRUE,digits=NA)
matbin <- function(m,nm) writeBin(as.double(m),file.path(work,nm),size=8,endian="little")
savej(plan,"plan.json"); matbin(raw,"raw.bin"); matbin(a,"A.bin")
p <- list(list(id="baseline",step=0,rows=list()),list(id="up",step=.1,rows=c(77,78)),
          list(id="down",step=-.1,rows=c(77,78)),list(id="marker",step=.1,rows=c(3,13,23)))
savej(p,"probes.json")
denom <- function(v) {z<-edgeR::calcNormFactors(edgeR::DGEList(counts=v)); z$samples$lib.size*z$samples$norm.factors}
d0 <- denom(raw)
write.table(d0,file.path(work,"outcome_denominator.csv"),row.names=FALSE,col.names=FALSE)
run_probes("clinical",work,work)
got <- read.csv(file.path(work,"probe_fractions.csv"))
dose <- read.csv(file.path(work,"dose.csv"))
rows <- function(probe,method) as.matrix(got[got$method==method & got$probe==probe,classes])
for (j in p) {
    v <- raw; rr <- unlist(j$rows)
    if(length(rr)) v[rr,] <- sweep(raw[rr,,drop=FALSE],2,exp(j$step*x),"*")
    di <- dose[dose$probe==j$id,]
    stopifnot(identical(di$sample_id,ids),max(abs(di$effective_relative_change-(denom(v)/d0-1)))<1e-12)
    z <- dtangle::dtangle(Y=log2(t(v[marks,])/denom(v)*1e6+1),references=log2(t(a)*1e6+1),
        markers=split(marks,owners),n_markers=rep(10,k),data_type="rna-seq")$estimates
    z <- sweep(z,2,fy,"/"); z<-z/rowSums(z)
    stopifnot(max(abs(z-rows(j$id,"dtangle_tmm")))<1e-10)
}
stopifnot(max(abs(rows("up","dtangle_marker")-rows("baseline","dtangle_marker")))<1e-12,
          max(abs(rows("down","dtangle_marker")-rows("baseline","dtangle_marker")))<1e-12,
          max(abs(rows("marker","dtangle_marker")-rows("baseline","dtangle_marker")))>1e-5)
savej(list(genes=g,samples=ids,panel_indices=marks,count_rows=marks),"imaging_plan.json")
matbin(raw,"imaging.bin"); run_probes("imaging",work,work)
image <- read.csv(file.path(work,"imaging_fractions.csv"))
for (m in c("dtangle_marker","dtangle_tmm")) {
  stopifnot(max(abs(rows("baseline",m)-as.matrix(image[image$method==m,classes])))<1e-11)
}
cat("PASS synthetic clinical and imaging fits; edgeR",as.character(packageVersion("edgeR")),
    "dtangle",as.character(packageVersion("dtangle")),"\n")
