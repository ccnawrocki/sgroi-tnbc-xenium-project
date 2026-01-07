renv::init(bare = T)
.libPaths()

# [1] "/Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/xenium/analysis/sgroi-tnbc/renv/library/macos/R-4.4/aarch64-apple-darwin20"
# [2] "/Users/ccn22/Library/Caches/org.R-project.R/R/renv/sandbox/macos/R-4.4/aarch64-apple-darwin20/f7156815" 

Sys.which("clang") # "/usr/bin/clang" 

renv::install("bioc::rhdf5")
renv::install("Matrix")

renv::install(packages = list("bioc::SingleCellExperiment", "Seurat", "bioc::SpatialExperiment"))
renv::install("bnprks/BPCells/r")

renv::install("pak")
renv::install("remotes")
remotes::install_github("scverse/anndataR", ref="084a187") # This is an older version that works with R 4.4.2

renv::install("korsunskylab/spatula")
renv::install("immunogenomics/presto@glmm")
renv::install("immunogenomics/singlecellmethods")

renv::install(packages = list("Bioc::limma", "Bioc::edgeR", "Bioc::DESeq2", "Bioc::glmGamPoi"))

renv::install("bioc::sparseMatrixStats")
install.packages("irlba", type = "source")

renv::install(packages = list("lme4", "lmerTest"))

pak::pkg_install(pkg = "Rfast")

renv::snapshot()
.rs.restartR()

# This is not everything, but this covers most packages we will use in R and
# will enable us to install pretty much any other package that we may need 
# later on in this project.

