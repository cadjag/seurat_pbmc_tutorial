#shiny app to generate gene-wise visualizations for seurat pbmc tutorial data
library(shiny)
library(shinythemes)
library(plotly)
library(periscope)



# load rds containing pre-analyzed data    
pbmc <- readRDS(file = "pbmc.rds")


# app ui

ui <- navbarPage("PBMC scRNA-seq Data",
                 theme = shinytheme("flatly"),
                 navbarMenu("Gene Plots",
                            tabPanel("Violin Plot",
                                     sidebarLayout(
                                       sidebarPanel(
                                         textOutput("violin_text"),
                                         br(),
                                         selectizeInput("gene_violin", #type to find gene symbol for plotting
                                                        label = "Gene Symbol:",
                                                        choices = NULL,
                                                        multiple = FALSE),
                                         
                                         downloadButton("download_violin", "Save Plot")
                                       ),
                                       
                                       mainPanel(
                                         plotOutput("violin_plot", height = "600px")
                                       )
                                       
                                     )
                            ),
                            tabPanel("Ridge Plot",
                                     sidebarLayout(
                                       sidebarPanel(
                                         textOutput("ridge_text"),
                                         br(),
                                         selectizeInput("gene_ridge", #type to find gene symbol for plotting
                                                        label = "Gene Symbol:",
                                                        choices = NULL,
                                                        multiple = FALSE),
                                         
                                         downloadButton("download_ridge", "Save Plot")
                                       ),
                                       
                                       mainPanel(
                                         plotOutput("ridge_plot", height = "600px")
                                       )
                                       
                                     )),
                            tabPanel("Feature Plot",
                                     sidebarLayout(
                                       sidebarPanel(
                                         textOutput("feature_text"),
                                         br(),
                                         selectizeInput("gene_feature", #type to find gene symbol for plotting
                                                        label = "Gene Symbol:",
                                                        choices = NULL,
                                                        multiple = FALSE),
                                         downloadButton("download_feature", "Save Plot")
                                       ),
                                       
                                       mainPanel(
                                         plotOutput("feature_plot", height = "600px")
                                       )
                                     )),
                            tabPanel("Multi-Feature Plot",
                                     sidebarLayout(
                                       sidebarPanel(
                                         textOutput("multifeature_text"),
                                         br(),
                                         selectizeInput("gene_multifeature", #type to find gene symbol for plotting
                                                        label = "Gene Symbols:",
                                                        choices = NULL,
                                                        multiple = TRUE,
                                                        options = list(maxItems = 2)),
                                         downloadButton("download_multifeature", "Save Plot")
                                       ),
                                       
                                       mainPanel(
                                         plotOutput("multifeature_plot", height = "600px")
                                       )
                                     )),
                            tabPanel("Dot Plot",
                                     sidebarLayout(
                                       sidebarPanel(
                                         textOutput("dot_text"),
                                         br(),
                                         selectizeInput("gene_dot", #type to find gene symbol for plotting
                                                        label = "Gene Symbols:",
                                                        choices = NULL,
                                                        selected = c("PPBP","LYZ","S100A9","IGLL5","GNLY","FTL","PF4","FTH1","GNG11","S100A8", "IL32", "CCL5", "CD8A"),
                                                        multiple = TRUE,
                                                        options = list(maxItems = 15)),
                                         downloadButton("download_dot", "Save Plot")
                                       ),
                                       
                                       mainPanel(
                                         plotOutput("dot_plot", height = "600px")
                                       )
                                     )
                            )),
                 navbarMenu("UMAP",
                            tabPanel("Interactive",
                                     mainPanel(
                                       plotlyOutput("umap_plot", width = "800px", height = "600px")
                                     )),
                            tabPanel("Static",
                                     plotOutput("umap_plot2", width = "800px", height = "600px"))
                 ),
                 
                 tabPanel("QC Plots",
                          mainPanel(
                            plotOutput("qc_violin_plot", height = "600px"),
                            plotOutput("qc_scatter_plot", height = "600px"),
                            plotOutput("qc_variable_features_plot", width = "50%", height = "600px")
                          ))
)



# app server

server <- function(input, output, session) {
  
  ################## violin plot
  # render text to describe example plot
  output$violin_text <- renderText({ "Enter a gene symbol to generate a violin plot. Try IL32 to see an example." })
  # update selectizeInput with gene list
  updateSelectizeInput(session, 'gene_violin', choices = gene_list, server = TRUE)  
  # draw violin plot based on selected gene
  violin_plot <- reactive({Seurat::VlnPlot(pbmc, features = input$gene_violin) +
      ggplot2::theme(legend.position="none", axis.title.x = element_blank())  })
  # render violin plot for display in shiny app
  output$violin_plot <- renderPlot({ violin_plot() })
  # render violin plot for downloading to png
  output$download_violin <- downloadHandler(filename = function() { paste ("violin_", input$gene_violin, ".png", sep = "") },
                                            content = function(file) { ggsave(file, violin_plot(), device = "png") } )
  
  ################## ridge plot
  # render text to describe example plot
  output$ridge_text <- renderText({ "Enter a gene symbol to generate a ridge plot. Try IL32 to see an example." })
  # update selectizeInput with gene list
  updateSelectizeInput(session, 'gene_ridge', choices = gene_list, server = TRUE)  
  # draw ridge plot based on selected gene
  ridge_plot <- reactive({Seurat::RidgePlot(pbmc, features = input$gene_ridge) +
      ggplot2::theme(legend.position="none", axis.title.y = element_blank())  })
  # render ridge plot for display in shiny app
  output$ridge_plot <- renderPlot({ ridge_plot() })
  # render ridge plot for downloading to png
  output$download_ridge <- downloadHandler(filename = function() { paste ("ridge_", input$gene_ridge, ".png", sep = "") },
                                           content = function(file) { ggsave(file, ridge_plot(), device = "png") } )
  
  ################## single-gene feature plot
  # render text to describe example plot
  output$feature_text <- renderText({ "Enter a gene symbol to generate a UMAP plot overlaid with gene expression. Try IL32 to see an example." })
  # update selectizeInput with gene list
  updateSelectizeInput(session, 'gene_feature', choices = gene_list, server = TRUE)  
  # draw feature plot based on selected gene
  feature_plot <- reactive({ Seurat::FeaturePlot(pbmc, features = input$gene_feature) })
  # render feature plot for display in shiny app
  output$feature_plot <- renderPlot({ feature_plot() })
  # render feature plot for downloading to png
  output$download_feature <- downloadHandler(filename = function() { paste ("feature_", input$gene_feature, ".png", sep = "") },
                                             content = function(file) { ggsave(file, feature_plot(), device = "png") } )
  
  ################## multi-gene feature plot
  # render text to describe example plot
  output$multifeature_text <- renderText({ "Enter two gene symbols to visualize colocalization. Try MS4A1 & CD79A to see an example." })
  # update selectizeInput with gene list
  updateSelectizeInput(session, 'gene_multifeature', choices = gene_list, server = TRUE)  
  # draw multifeature plot based on selected gene
  multifeature_plot <- reactive({ Seurat::FeaturePlot(pbmc, features = input$gene_multifeature, blend = TRUE) })
  # render multifeature plot for display in shiny app
  output$multifeature_plot <- renderPlot({ multifeature_plot() })
  # render multifeature plot for downloading to png
  output$download_multifeature <- downloadHandler(filename = function() { paste ("multifeature_", input$gene_multifeature, ".png", sep = "") },
                                                  content = function(file) { ggsave(file, multifeature_plot(), device = "png",
                                                                                    width = 16, height = 6, units = "in") } )  
  ################## dot plot  
  # render text to describe example plot
  output$dot_text <- renderText({ "Enter up to 15 gene symbols to generate a dot plot. 
                                   Try PPBP, LYZ, S100A9, IGLL5, GNLY, FTL, PF4, FTH1, GNG11, and S100A8 to see the top 10 most variable genes." })
  # update selectizeInput with gene list
  updateSelectizeInput(session, 'gene_dot', choices = gene_list, server = TRUE)  
  # draw dot plot based on selected gene
  dot_plot <- reactive({ Seurat::DotPlot(pbmc, features = input$gene_dot) + 
      ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) })
  # render dot plot for display in shiny app
  output$dot_plot <- renderPlot({ dot_plot() })
  # render dot plot for downloading to png
  output$download_dot <- downloadHandler(filename = function() { paste ("dot_", input$gene_dot, ".png", sep = "") },
                                         content = function(file) { ggsave(file, dot_plot(), device = "png", bg = "white",
                                                                           width = 12, height = 6, units = "in") } ) 
  
  ################## umap
  output$umap_plot <- renderPlotly({ baseplot <- Seurat::DimPlot(pbmc, reduction = "umap")
  baseplot <- Seurat::HoverLocator(plot = baseplot, 
                                   information = Seurat::FetchData(pbmc, vars = c("ident", "percent.mt", "nFeature_RNA", "nCount_RNA")))
  baseplot })
  output$umap_plot2 <- renderPlot({ Seurat::DimPlot(pbmc, reduction = "umap") })
  
  ################## qc plots  
  output$qc_violin_plot <- renderPlot({ qc_violin_plot })
  
  output$qc_scatter_plot <- renderPlot({ qc_scatter_plot1 + qc_scatter_plot2 })
  
  output$qc_variable_features_plot <- renderPlot({
    # Identify the 10 most highly variable genes
    top10 <- head(Seurat::VariableFeatures(pbmc), 10)
    
    # plot variable features with labels
    plot1 <- Seurat::VariableFeaturePlot(pbmc)
    plot1 <- Seurat::LabelPoints(plot = plot1, points = top10, repel = TRUE)
    plot1
  })
  
}
 
# run app

shinyApp(ui = ui, server = server)