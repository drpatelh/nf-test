#!/usr/bin/env Rscript
packages <- c("devtools", "BiocManager")
suppressMessages(invisible(lapply(packages, library, character.only = TRUE, quietly=TRUE)))
Sys.setenv(CONDA_BUILD_SYSROOT="/")
devtools::install_github("GreenleafLab/ArchR", ref="release_1.0.1", repos = BiocManager::repositories(), upgrade="always", quiet=TRUE, build_opts=c("--no-build-vignettes"))
devtools::install_github("immunogenomics/presto",upgrade="always", quiet=TRUE, build_opts=c("--no-build-vignettes"))


