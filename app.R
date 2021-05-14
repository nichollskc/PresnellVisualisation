library(shiny)
library(ggplot2)
library(plotly)
library(dplyr)
library(heatmaply)

setwd("~/docs/PhD/PresnellVisualisation/")
source("biclust_utils.R")
source("presnell.R")

SSLB_result <- read_thresholded_factor_matrices(main_run)
rownames(SSLB_result$B) <- gene_symbols
rownames(SSLB_result$X) <- sample_names
colnames(SSLB_result$X) <- paste("factor", 1:SSLB_result$K, sep="_")
colnames(SSLB_result$B) <- paste("factor", 1:SSLB_result$K, sep="_")

sample_info_with_fac <- cbind(sample_info, SSLB_result$X)

ui <- fluidPage(
  titlePanel("Biclustering on Presnell sorted blood cell dataset"),
  sidebarLayout(
    sidebarPanel(
      numericInput(inputId="factor", label="Factor",
                   value=11, min=1, max=SSLB_result$K, step=1),
    ),
    mainPanel(
      plotlyOutput("factorcontribution_heatmap"),
      tableOutput("nz_samples_table"),
      tableOutput("nz_genes_table"),
    )
  )
)

server <- function(input, output) {
  factor_contribution <- reactive({
    contribution <- calc_factor_contribution(SSLB_result, input$factor)
    contribution
  })
  nz_samples <- reactive({
    (SSLB_result$X[, input$factor] != 0)
  })
  nz_genes <- reactive({
    (SSLB_result$B[, input$factor] != 0)
  })
  nz_factor_contribution <- reactive({
    factor_contribution()#[nz_samples(), nz_genes()]
  })
  output$factorcontribution_heatmap <- renderPlotly({
    heatmaply(nz_factor_contribution())
  })
  
  output$nz_samples_table <- renderTable({
    sample_info[nz_samples()]
  })
  
  output$nz_genes_table <- renderTable({
    nz_genes()
  })
}

shinyApp(ui = ui, server = server)
