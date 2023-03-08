# project setup

renv::activate()
renv::install("markdown")
renv::install("BiocManager")
renv::install("tidyverse")
renv::install("Seurat")
BiocManager::install("limma")

renv::snapshot()