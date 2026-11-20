# Select a documented R/Bioconductor pair locally, without a live config.yaml call.
# Mapping: https://bioconductor.org/about/release-announcements/
bioc_for_r <- function(version=getRversion()) {
  key <- paste(as.integer(strsplit(as.character(version), "\\.")[[1]][1:2]), collapse=".")
  versions <- c("4.4"="3.20", "4.5"="3.22", "4.6"="3.23")
  if (!(key %in% names(versions))) {
    stop("Unsupported R version: ", version, ". Use R 4.4, 4.5 or 4.6; do not guess a Bioconductor release.")
  }
  unname(versions[[key]])
}

package_repos <- function(bioc,
                          cran=Sys.getenv("PAPER3_CRAN", "https://cloud.r-project.org"),
                          mirror=Sys.getenv("PAPER3_BIO_MIRROR", "https://bioconductor.posit.co")) {
  urls <- c(cran, mirror)
  if (any(!startsWith(urls, "https://"))) stop("Package mirrors must use HTTPS; TLS verification remains enabled.")
  root <- paste0(sub("/+$", "", mirror), "/packages/", bioc)
  c(BioCsoft=paste0(root, "/bioc"), BioCann=paste0(root, "/data/annotation"),
    BioCexp=paste0(root, "/data/experiment"), CRAN=sub("/+$", "", cran))
}

install_main <- function(args) {
  stopifnot(length(args) == 2L)
  dir.create(args[1], recursive=TRUE, showWarnings=FALSE)
  lib <- normalizePath(args[1], mustWork=TRUE)
  .libPaths(c(lib, .libPaths()))
  bioc <- bioc_for_r()
  repos <- package_repos(bioc)
  options(timeout=600, repos=repos)
  if (.Platform$OS.type == "windows") {
    options(pkgType="binary", install.packages.compile.from.source="never")
  } else {
    cmake <- Sys.which("cmake")
    if (!nzchar(cmake)) stop("CMake missing from PATH; rerun the launcher with --setup.")
    cmake_text <- system2(cmake, "--version", stdout=TRUE, stderr=TRUE)
    cat(paste(cmake_text, collapse="\n"), "\n")
    if (!is.null(attr(cmake_text, "status"))) stop("CMake --version failed.")
  }
  cat("R:", as.character(getRversion()), "Bioconductor repository:", bioc, "\n")
  print(repos)
  writeLines(c(paste("R", getRversion()), paste("Bioconductor repository", bioc),
               paste(names(repos), repos)), file.path(lib, "paper3_install_repositories.txt"))

  installed <- function() rownames(installed.packages(lib.loc=.libPaths()))
  check_packages <- function(packages) {
    ok <- vapply(packages, requireNamespace, logical(1), quietly=TRUE)
    if (!all(ok)) stop("Packages missing or unloadable: ", paste(packages[!ok], collapse=", "),
                      ". See the package-specific ERROR above; analysis has not started.")
  }
  # A failed staged install removes the package; successful installs are reused.
  cran <- c("jsonlite", "DEoptimR", "Rfast", "quadprog", "ashr", "SQUAREM")
  missing <- setdiff(cran, installed())
  if (length(missing)) install.packages(missing, lib=lib, Ncpus=2)
  check_packages(cran)
  # Annotation packages have source distributions even on Windows.
  if (!("GenomeInfoDbData" %in% installed())) {
    install.packages("GenomeInfoDbData", lib=lib, type="source")
  }
  bio <- c("edgeR", "DESeq2", "SingleCellExperiment", "qvalue")
  missing <- setdiff(bio, installed())
  if (length(missing)) install.packages(missing, lib=lib, Ncpus=2)
  check_packages(c("GenomeInfoDbData", bio))
  if (!requireNamespace("vashr", quietly=TRUE)) {
    install.packages(file.path(args[2], "vashr.tar.gz"), lib=lib, repos=NULL, type="source")
  }
  dtangle_path <- find.package("dtangle", quiet=TRUE)
  dtangle_version <- if (length(dtangle_path)) {
    read.dcf(file.path(dtangle_path, "DESCRIPTION"), fields="Version")[1, 1]
  } else ""
  if (dtangle_version != "2.0.9") {
    install.packages(file.path(args[2], "dtangle_2.0.9.tar.gz"), lib=lib, repos=NULL, type="source")
  }
  # Upstream omits dependency fields; install its imports explicitly above.
  install.packages(file.path(args[2], "MEAD"), lib=lib, repos=NULL, type="source")
  needed <- c(cran, bio, "dtangle", "MEAD", "vashr")
  check_packages(needed)
  print(sapply(needed, function(x) as.character(getNamespaceVersion(x))))
  writeLines(capture.output(sessionInfo()), file.path(lib, "paper3_install_session.txt"))
}

if (sys.nframe() == 0L) install_main(commandArgs(trailingOnly=TRUE))
