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
SSLB_result <- add_row_col_names(SSLB_result, gene_symbols, sample_names)
sample_info_with_fac <- generate_sample_info_with_fac(SSLB_result, sample_info)
gene_info_with_fac <- generate_gene_info_with_fac(SSLB_result, gene_info)

pathway_info <- gene_info %>% select(-one_of("Ensembl", "ensembl", "GeneSymbol", "annotated"))
pathway_enrichment <- calculate_pathway_enrichment(SSLB_result, pathway_info)
pathway_enrichment <- add_odds_ratio_and_counts(pathway_enrichment, SSLB_result, pathway_info)


ui <- fluidPage(
  titlePanel("Biclustering on Presnell sorted blood cell dataset"),
  sidebarLayout(
    sidebarPanel(
      numericInput(inputId="factor", label="Factor",
                   value=11, min=1, max=SSLB_result$K, step=1),
    ),
    mainPanel(
      dataTableOutput("pathways_table"),
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
             factor_column, abs_factor_column)
  })
  sorted_samples_nz <- reactive({
    abs_factor_column <- paste0("abs_factor_", input$factor)
    sorted_samples() %>% filter(abs(!!as.name(abs_factor_column)) > 0)
  })

  sorted_genes <- reactive({
    enriched_pathway_names <- enriched_pathways()$PathwayName
    print(enriched_pathway_names)

    factor_column <- paste0("factor_", input$factor)
    abs_factor_column <- paste0("abs_factor_", input$factor)
    gene_info_with_fac %>%
      arrange(desc(abs_factor_column)) %>%
      select(GeneSymbol, ensembl,
             factor_column, abs_factor_column,
             enriched_pathway_names)
  })
  sorted_genes_nz <- reactive({
    sorted_genes()
    abs_factor_column <- paste0("abs_factor_", input$factor)
    sorted_genes() %>% filter(abs(!!as.name(abs_factor_column)) > 0)
  })

  enriched_pathways <- reactive({
    print(colnames(pathway_enrichment))
    pathway_enrichment %>%
      arrange(CameraPR_qvalue) %>%
      filter(CameraPR_qvalue < 0.05, FactorIndex == input$factor) %>%
      select(PathwayName, GenesInPathway, GenesInIntersection,
             OddsRatio, CameraPR_qvalue, GeneSetTest_qvalue)
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
  output$nz_genes_table <- renderDataTable({
    sorted_genes_nz()
  })
  output$pathways_table <- renderDataTable({
    enriched_pathways()
  })
}

shinyApp(ui = ui, server = server)

