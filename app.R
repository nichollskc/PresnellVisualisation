library(shiny)
library(ggplot2)
library(plotly)
library(dplyr)
library(heatmaply)
library(limma)

setwd("~/docs/PhD/PresnellVisualisation/")
source("biclust_utils.R")
source("presnell.R")

SSLB_result <- read_thresholded_factor_matrices(main_run)
add_row_col_names(SSLB_result, gene_symbols, sample_names)
sample_info_with_fac <- generate_sample_info_with_fac(SSLB_result, sample_info)
gene_info_with_fac <- generate_gene_info_with_fac(SSLB_result, gene_info)

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
      dataTableOutput("nz_genes_table"),
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

  sorted_genes <- reactive({
    enriched_pathway_names <- enriched_pathways() %>%
      filter(pvalue < 0.05) %>%
      colnames()

    factor_column <- paste0("factor_", input$factor)
    abs_factor_column <- paste0("abs_factor_", input$factor)
    gene_info_with_fac %>%
      arrange(desc(abs_factor_column)) %>%
      select("ensembl",
             factor_column, abs_factor_column,
             enriched_pathway_names,
             starts_with("factor_"), starts_with("abs_factor_"))
  })
  sorted_genes_nz <- reactive({
    abs_factor_column <- paste0("abs_factor_", input$factor)
    sorted_genes() %>% filter(abs(!!as.name(abs_factor_column)) > 0)
  })

  enriched_pathways <- rownames(gene_info_with_fac %>% filter(HALLMARK_INTERFERON_ALPHA_RESPONSE != 0))


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
  output$nz_genes_table <- renderDataTable({
    sorted_genes_nz()
  })
}

shinyApp(ui = ui, server = server)

