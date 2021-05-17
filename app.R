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

abs_X <- abs(SSLB_result$X)
colnames(abs_X) <- paste0("abs_", colnames(abs_X))

sample_info_with_fac <- cbind(sample_info, SSLB_result$X)
sample_info_with_fac <- cbind(sample_info_with_fac, abs_X)

ui <- fluidPage(
  titlePanel("Biclustering on Presnell sorted blood cell dataset"),
  sidebarLayout(
    sidebarPanel(
      numericInput(inputId="factor", label="Factor",
                   value=11, min=1, max=SSLB_result$K, step=1),
    ),
    mainPanel(
      plotlyOutput("factorcontribution_heatmap"),
      dataTableOutput("nz_samples_table"),
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
  sorted_samples <- reactive({
    factor_column <- paste0("factor_", input$factor)
    abs_factor_column <- paste0("abs_factor_", input$factor)
    sample_info_with_fac %>%
      arrange(desc(abs_factor_column)) %>%
      select("sample_description", "sex", "disease", "cell", "treatment", "age",
             factor_column, abs_factor_column,
             starts_with("factor_"), starts_with("abs_factor_"))
  })
  sorted_samples_nz <- reactive({
    abs_factor_column <- paste0("abs_factor_", input$factor)
    sorted_samples() %>% filter(abs(!!as.name(abs_factor_column)) > 0)
  })
  nz_factor_contribution_sorted <- reactive({
    factor_cont <- factor_contribution()
    factor_cont[rownames(sorted_samples_nz()), nz_genes()]
  })
  output$factorcontribution_heatmap <- renderPlotly({
    heatmaply(nz_factor_contribution_sorted(), Rowv=FALSE, Colv=FALSE)
  })
  
  output$nz_samples_table <- renderDataTable({
    sorted_samples_nz()
  })
  
  output$nz_genes_table <- renderTable({
    nz_genes()
  })
}

shinyApp(ui = ui, server = server)
