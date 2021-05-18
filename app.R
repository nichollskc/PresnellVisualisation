library(BiocManager)
options(repos = BiocManager::repositories())

library(shiny)
library(ggplot2)
library(plotly)
library(dplyr)
library(tidyr)
library(heatmaply)
library(limma)
library(shinydashboard)
library(shinyjs)

source("biclust_utils.R")
source("presnell.R")

SSLB_result <- read_thresholded_factor_matrices(main_run)
SSLB_result <- add_row_col_names(SSLB_result, gene_symbols, sample_names)
sample_info_with_fac <- generate_sample_info_with_fac(SSLB_result, sample_info)
gene_info_with_fac <- generate_gene_info_with_fac(SSLB_result, gene_info)

pathway_info <- gene_info %>% select(-one_of("Ensembl", "ensembl", "GeneSymbol", "annotated"))
load("pathway_enrichment.Rda")
#pathway_enrichment <- calculate_pathway_enrichment(SSLB_result, pathway_info)
#pathway_enrichment <- add_odds_ratio_and_counts(pathway_enrichment, SSLB_result, pathway_info)

total_sample_counts <- count_samples_by_type(sample_info)


sample_heatmap_output <- box(title="Sample types",
    collapsible=TRUE, solidHeader = TRUE, status="primary", width=NULL,
    div(paste0("Types of samples included in the factor, broken down by disease, ",
               "cell type and sex. The numbers in the cells give the percentage of samples",
               " of that type which are contained in this factor. The darker the colour ", 
               "(red for sex, or blue for cell type), the higher this percentage.",
               "Grey represents sample types not present in the dataset.")),
    plotlyOutput("sample_heatmap"),
)

pathways_table_output <- box(title="Enriched pathways",
    collapsible=TRUE, solidHeader = TRUE, status="primary", width=NULL,
    div(paste0("Pathways with strongest enrichment (measured by the q-value returned by ",
               "CameraPR) in this factor.")),
    div(style = 'overflow-x: scroll', dataTableOutput("pathways_table"))
)

factorcontribution_output <- box(title="Factor contribution",
    div(paste0("Heatmap showing the values this factor contributes to the overall ",
               "reconstruction of the original matrix. Constructed by $x_k b_k^T$ ",
               "where $x_k$ is the length n vector giving sample loadings and ",
               "$b_k$ is the length p vector giving gene loadings.")),
    plotlyOutput("factorcontribution_heatmap"),
    collapsible=TRUE, solidHeader = TRUE, status="primary", width=NULL,
)

sample_table_output <- box(title="Samples table",
    div(style = 'overflow-x: scroll', dataTableOutput("nz_samples_table")),
    collapsible=TRUE, solidHeader = TRUE, status="primary", width=NULL,
)

gene_table_output <- box(title="Genes table",
    div(style = 'overflow-x: scroll', dataTableOutput("nz_genes_table")),
    collapsible=TRUE, solidHeader = TRUE, status="primary", width=NULL,
)

ui <- dashboardPage(
  dashboardHeader(title="Biclustering on Presnell sorted blood cell dataset"),
  dashboardSidebar(
    numericInput(inputId="factor", label="Factor",
                 value=12, min=1, max=SSLB_result$K, step=1),
    htmlOutput("summary_text")
  ),
  dashboardBody(
    # Allow custom JS script to collapse box when you click on header
    useShinyjs(),
    # Darken header on hover to make it more obvious you can click
    tags$style(HTML("
                    .box-header {
                    padding: 0 10px 0 0;
                    }
                    .box-header:hover {
                    filter: brightness(80%);
                    }
                    .box-header h3 {
                    width: 100%;
                    padding: 10px;
                    }")),
    fluidRow(
      column(width=6,
        sample_heatmap_output,
        sample_table_output,
      ),
      column(width=6,
        pathways_table_output,
        factorcontribution_output,
        gene_table_output,
      )
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
      arrange(desc(!!as.name(factor_column))) %>%
      select("sample_description", "sex", "short_disease", "cell", "treatment", "age",
             factor_column, abs_factor_column)
  })
  sorted_samples_nz <- reactive({
    abs_factor_column <- paste0("abs_factor_", input$factor)
    sorted_samples() %>% filter(abs(!!as.name(abs_factor_column)) > 0)
  })

  sorted_genes <- reactive({
    enriched_pathway_names <- enriched_pathways()$PathwayName

    factor_column <- paste0("factor_", input$factor)
    abs_factor_column <- paste0("abs_factor_", input$factor)
    gene_info_with_fac %>%
      arrange(desc(!!as.name(factor_column))) %>%
      select(GeneSymbol, ensembl,
             factor_column, abs_factor_column,
             enriched_pathway_names)
  })
  sorted_genes_nz <- reactive({
    abs_factor_column <- paste0("abs_factor_", input$factor)
    sorted_genes() %>%
      filter(abs(!!as.name(abs_factor_column)) > 0)
  })

  enriched_pathways <- reactive({
    pathway_enrichment %>%
      arrange(CameraPR_qvalue) %>%
      filter(CameraPR_qvalue < 0.05, FactorIndex == input$factor) %>%
      select(PathwayName, GenesInPathway, GenesInIntersection,
             OddsRatio, CameraPR_qvalue, GeneSetTest_qvalue)
  })

  nz_factor_contribution_sorted <- reactive({
    factor_cont <- factor_contribution()
    factor_cont[rownames(sorted_samples_nz()), rownames(sorted_genes_nz())]
  })
  
  
  output$sample_heatmap <- renderPlotly({
    factor_sample_counts <- count_samples_by_type(sorted_samples_nz())
    factor_sample_proportions <- calculate_proportions_sample_types(sample_info_with_fac,
                                                                    input$factor)
    sex_hm <- sample_heatmap(factor_sample_proportions$by_sex_disease, "firebrick")
    cell_hm <- sample_heatmap(factor_sample_proportions$by_cell_disease, "dodgerblue")
    
    subplot(cell_hm, sex_hm, nrows=2, heights=c(5/7, 2/7))
  })
  output$factorcontribution_heatmap <- renderPlotly({
    heatmaply(nz_factor_contribution_sorted(),
              # Don't apply clustering to rows and columns since we've already sorted
              Rowv=FALSE, Colv=FALSE,
              scale_fill_gradient_fun = scale_fill_gradient2(low="blue", high="red", midpoint=0))
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
  output$summary_text <- renderUI({
    HTML(
      sprintf("Total samples: %d<br>Total genes: %d",
              sum(nz_samples()),
              sum(nz_genes()))
    )
  })
  
  # JS code to collapse box when you click on header
  runjs("
          $('.box').on('click', '.box-header h3', function() {
        $(this).closest('.box')
        .find('[data-widget=collapse]')
        .click();
        });")
  
}

shinyApp(ui = ui, server = server)

