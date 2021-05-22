library(BiocManager)
options(repos = BiocManager::repositories())

library(shiny)
library(ggplot2)
library(plotly)
library(dplyr)
library(tidyr)
library(heatmaply)
library(limma)

#setwd("~/docs/PhD/PresnellVisualisation")
source("biclust_utils.R")
source("plot_utils.R")
source("presnell.R")

#################################################################################################
# Loading data                                                                                  #
#################################################################################################
SSLB_result <- read_thresholded_factor_matrices(main_run)
SSLB_result <- add_row_col_names(SSLB_result, gene_symbols, sample_names)
sample_info_with_fac <- generate_sample_info_with_fac(SSLB_result, sample_info)
gene_info_with_fac <- generate_gene_info_with_fac(SSLB_result, gene_info)

pathway_info <- gene_info %>% select(-one_of("Ensembl", "ensembl", "GeneSymbol", "annotated"))
load("pathway_enrichment.Rda")
#pathway_enrichment <- calculate_pathway_enrichment(SSLB_result, pathway_info)
#pathway_enrichment <- add_odds_ratio_and_counts(pathway_enrichment, SSLB_result, pathway_info)

total_sample_counts <- count_samples_by_type(sample_info)
factor_contribution_maxes <- apply(SSLB_result$X, 2, max) * apply(SSLB_result$B, 2, max)
total_counts_sex_disease <- counts_by_variables(sample_info_with_fac, "sex", "disease")
total_counts_cell_disease <- counts_by_variables(sample_info_with_fac, "cell", "disease")

#################################################################################################
# Server - processing data and generating plots                                                 #
#################################################################################################
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
  
  #################################################################################################
  # Defining plots                                                                                #
  #################################################################################################
  output$sample_heatmap <- renderPlotly({
    factor_counts_sex_disease <- counts_by_variables(sample_info_with_fac, "sex", "disease", input$factor)
    factor_counts_cell_disease <- counts_by_variables(sample_info_with_fac, "cell", "disease", input$factor)
    
    sex_hm <- sample_heatmap(total_counts_sex_disease, factor_counts_sex_disease, "Reds")
    cell_hm <- sample_heatmap(total_counts_cell_disease, factor_counts_cell_disease, "Blues")
    
    subplot(cell_hm, sex_hm, nrows=2, heights=c(5/7, 2/7))
  })
  output$factorcontribution_heatmap <- renderPlotly({

    if (input$factorcontribution_scale == "This factor") {
      max_value <- max(abs(nz_factor_contribution_sorted()))
    } else {
      max_value <- max(factor_contribution_maxes)
    }
    limits <- c(-max_value, max_value)
    
    if (input$factorcontribution_trans == "No transformation") {
      transformation = "identity"
      breaks = neat_symmetric_cbar_breaks(max_value)
    } else {
      transformation = "pseudo_log"
      breaks = neat_symmetric_cbar_breaks_log(max_value)
    }
    
    heatmaply(nz_factor_contribution_sorted(),
              # Don't apply clustering to rows and columns since we've already sorted
              Rowv=FALSE, Colv=FALSE,
              scale_fill_gradient_fun = scale_fill_gradient2(low="blue", high="red", midpoint=0,
                                                             trans=transformation,
                                                             limits=limits,
                                                             breaks=breaks,
              ))
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

server

