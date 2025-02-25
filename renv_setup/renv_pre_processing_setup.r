#run this the first time
renv::init()
#ignore warnings about restarting

#run this when you launch the docker instance again
renv::restore()

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

install.packages("renv")
install.packages(c("ggplot2"))
install.packages(c("dplyr"))
renv::install("Seurat@4.4.0")
renv::install("stringr")
renv::install("babelgene")
renv::install("devtools")
renv::install("HGNChelper")
renv::install("openxlsx")
renv::install("BiocManager")
BiocManager::install("SummarizedExperiment")
BiocManager::install("SingleCellExperiment")
BiocManager::install("zellkonverter")



#run this after launching a new docker image
renv::snapshot(type = "all", prompt = FALSE)


