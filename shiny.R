#shiny app to generate gene-wise visualizations for seurat pbmc tutorial data
library(shiny)
library(shinythemes)


# environment setup / Seurat analysis & processing
  # load data and initialize Seurat object
    pbmc.data <- Seurat::Read10X(data.dir = "./data/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/")
    pbmc <- Seurat::CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
  # extract gene list from Seurat object
    gene_list <- rownames(pbmc)
  # add mitochondrial RNA % column to QC stats dataframe
    pbmc[["percent.mt"]] <- Seurat::PercentageFeatureSet(pbmc, pattern = "^MT-")
  # subset cells for further analysis based on QC metrics
    pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  # normalize data
    pbmc <- Seurat::NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  # find top 2000 variable features
    pbmc <- Seurat::FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  # scale data
    pbmc <- Seurat::ScaleData(pbmc, features = gene_list)
  # run PCA
    pbmc <- Seurat::RunPCA(pbmc, features = Seurat::VariableFeatures(object = pbmc))
  # cluster cells
    pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:10)
    pbmc <- Seurat::FindClusters(pbmc, resolution = 0.5)
  # rename clusters
    new.cluster.ids <- c("Naive CD4+ T", "CD14+ Mono", "Memory CD4+ T", "B", "CD8+ T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
    names(new.cluster.ids) <- levels(pbmc)
    pbmc <- Seurat::RenameIdents(pbmc, new.cluster.ids)
  # run umap
    pbmc <- Seurat::RunUMAP(pbmc, dims = 1:10)
  

# app ui

ui <- navbarPage("PBMC scRNA-seq Data",
        theme = shinytheme("flatly"),
        navbarMenu("Gene Plots",
        tabPanel("Violin Plot",
          sidebarLayout(
            sidebarPanel(
              selectizeInput("gene_violin", #type to find gene symbol for plotting
                             label = "Gene Symbol:",
                             choices = gene_list,
                             selected = "IL32",
                             multiple = FALSE,
                             options = NULL)
                        ),
            
            mainPanel(
              plotOutput("violin_plot", height = "600px")
                      )
          
                    )
        ),
        tabPanel("Ridge Plot",
                 sidebarLayout(
                   sidebarPanel(
                     selectizeInput("gene_ridge", #type to find gene symbol for plotting
                                    label = "Gene Symbol:",
                                    choices = gene_list,
                                    selected = "IL32",
                                    multiple = FALSE,
                                    options = NULL)
                   ),
                   
                   mainPanel(
                     plotOutput("ridge_plot", height = "600px")
                   )
                   
                 )),
        tabPanel("Feature Plot",
                 sidebarLayout(
                   sidebarPanel(
                     selectizeInput("gene_feature", #type to find gene symbol for plotting
                                    label = "Gene Symbol:",
                                    choices = gene_list,
                                    selected = "IL32",
                                    multiple = FALSE,
                                    options = NULL)
                   ),
                   
                   mainPanel(
                     plotOutput("feature_plot", height = "600px")
                   )
                   
                 ))
        ),
        tabPanel("Interactive UMAP"),
        tabPanel("QC Plots")
)



# app server

server <- function(input, output) {

  output$violin_plot <- renderPlot({ Seurat::VlnPlot(pbmc, features = input$gene_violin) +
                                     ggplot2::theme(legend.position="none",
                                                    axis.title.x = element_blank()) })
  
  output$ridge_plot <- renderPlot({ Seurat::RidgePlot(pbmc, features = input$gene_ridge) +
                                    ggplot2::theme(legend.position="none",
                                                   axis.title.y = element_blank()) })
  
  output$feature_plot <- renderPlot({ Seurat::FeaturePlot(pbmc, features = input$gene_feature)})
  
}

# run app

shinyApp(ui = ui, server = server)