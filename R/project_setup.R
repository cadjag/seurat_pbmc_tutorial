# project setup

renv::activate()
renv::install("markdown")
renv::install("BiocManager")
renv::install("tidyverse")
renv::install("Seurat")
renv::install("Matrix")
renv::install("shiny")
renv::install("shinythemes")
renv::install("plotly")
renv::install("periscope")
renv::install("rsconnect")
BiocManager::install("limma")

reticulate::py_install(packages = 'umap-learn')

renv::snapshot()
