# Offline tests: no package downloads and no changes to installed R libraries.
source("install_r.R")
stopifnot(bioc_for_r("4.4.2") == "3.20", bioc_for_r("4.5.3") == "3.22",
          bioc_for_r("4.6.0") == "3.23")
stopifnot(inherits(try(bioc_for_r("4.3.3"), silent=TRUE), "try-error"))
repos <- package_repos("3.22", "https://cloud.r-project.org/", "https://bioconductor.posit.co/")
stopifnot(repos[["BioCsoft"]] == "https://bioconductor.posit.co/packages/3.22/bioc",
          repos[["BioCann"]] == "https://bioconductor.posit.co/packages/3.22/data/annotation")
stopifnot(inherits(try(package_repos("3.22", cran="http://invalid.example"), silent=TRUE), "try-error"))

exercise_installer <- function(initial) {
  e <- new.env(parent=globalenv())
  sys.source("install_r.R", envir=e)
  e$packages <- initial
  e$calls <- list()
  e$getRversion <- function() package_version("4.5.3")
  e$.Platform <- list(OS.type="unix")
  e$Sys.which <- function(...) "/test/cmake"
  e$system2 <- function(...) "cmake version 3.31.10"
  e$installed.packages <- function(...) matrix("1", nrow=length(e$packages), ncol=1,
                                              dimnames=list(e$packages, "Version"))
  e$requireNamespace <- function(package, ...) package %in% e$packages
  e$find.package <- function(package, ...) if (package %in% e$packages) "/test/dtangle" else character()
  e$read.dcf <- function(...) matrix("2.0.9", 1, 1)
  e$getNamespaceVersion <- function(...) package_version("1.0.0")
  e$install.packages <- function(pkgs, ...) {
    e$calls[[length(e$calls)+1L]] <- list(pkgs=pkgs, ...)
    names <- sub("_2\\.0\\.9$", "", sub("\\.tar\\.gz$", "", basename(pkgs)))
    e$packages <- union(e$packages, names)
  }
  path <- tempfile("paper3_R_setup_test_")
  old_options <- options()
  old_paths <- .libPaths()
  on.exit({options(old_options); .libPaths(old_paths); unlink(path, recursive=TRUE)}, add=TRUE)
  e$install_main(c(path, "/test/vendor"))
  requested <- unlist(lapply(e$calls, `[[`, "pkgs"))
  stopifnot(!any(grepl("BiocManager", requested)),
            file.exists(file.path(path, "paper3_install_repositories.txt")))
  requested
}
requested <- exercise_installer(character())
stopifnot(all(c("Rfast", "edgeR", "DESeq2", "SingleCellExperiment", "qvalue") %in% requested))
requested <- exercise_installer(c("jsonlite", "DEoptimR", "Rfast", "quadprog", "ashr", "SQUAREM",
                                  "GenomeInfoDbData", "edgeR", "DESeq2", "SingleCellExperiment", "qvalue",
                                  "vashr", "dtangle", "MEAD"))
stopifnot(identical(requested, "/test/vendor/MEAD"))
cat("PASS: R release selection, HTTPS guard, install ordering, successful-package reuse\n")
