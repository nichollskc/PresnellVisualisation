library(BiocManager)
options(repos = BiocManager::repositories())

library(shiny)
library(ggplot2)
library(plotly)
library(dplyr)
library(tidyr)
library(heatmaply)
library(limma)
library(withr)

#setwd("~/docs/PhD/PresnellVisualisation")
source("biclust_utils.R")
source("plot_utils.R")

#################################################################################################
# Loading/generating data                                                                       #
#################################################################################################
load_data_from_file <- TRUE
if (load_data_from_file) {
  load("presnellVisualisation.Rda")  
} else {
  source("presnell.R")
  SSLB_result <- read_thresholded_factor_matrices(main_run)
  SSLB_result <- add_row_col_names(SSLB_result, gene_symbols, sample_names)
  sample_info_with_fac <- generate_sample_info_with_fac(SSLB_result, sample_info)
  gene_info_with_fac <- generate_gene_info_with_fac(SSLB_result, gene_info)
  gene_info_with_fac <- add_variance_from_factors(gene_info_with_fac, SSLB_result, Y)
  
  pathway_info <- gene_info %>% select(-one_of("Ensembl", "ensembl", "GeneSymbol", "annotated"))
  load("pathway_enrichment.Rda")
  #pathway_enrichment <- calculate_pathway_enrichment(SSLB_result, pathway_info)
  #pathway_enrichment <- add_odds_ratio_and_counts(pathway_enrichment, SSLB_result, pathway_info)
  
  factor_contribution_maxes <- apply(SSLB_result$X, 2, max) * apply(SSLB_result$B, 2, max)
  total_counts_sex_disease <- counts_by_variables(sample_info_with_fac, "sex", "short_disease")
  total_counts_cell_disease <- counts_by_variables(sample_info_with_fac, "cell", "short_disease")
  
  factor_info <- generate_factor_info(gene_info_with_fac, sample_info_with_fac, SSLB_result)
  
  #save.image(file="presnellVisualisation.Rda")
}

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
    prop_var_column <- paste0("prop_var_factor_", input$factor)
    var_column <- paste0("var_factor_", input$factor)
    gene_info_with_fac %>%
      arrange(desc(!!as.name(factor_column))) %>%
      select(GeneSymbol, ensembl,
             factor_column, abs_factor_column, prop_var_column, var_column, total_variance,
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
    factor_counts_sex_disease <- counts_by_variables(sample_info_with_fac, "sex", "short_disease", input$factor)
    factor_counts_cell_disease <- counts_by_variables(sample_info_with_fac, "cell", "short_disease", input$factor)
    
    sex_hm <- sample_heatmap(total_counts_sex_disease, factor_counts_sex_disease, "firebrick")
    cell_hm <- sample_heatmap(total_counts_cell_disease, factor_counts_cell_disease, "dodgerblue", show_x_labels=FALSE)
    
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
              row_side_colors = sorted_samples_nz() %>% select(cell, sex, short_disease, treatment),
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
  output$gene_importance_plot <- renderPlotly({
    factor_gene_info <- sorted_genes_nz()
    factor_gene_info$factor_loading <- factor_gene_info[[paste0("factor_",
                                                                input$factor)]]
    factor_gene_info$prop_var_explained <- factor_gene_info[[paste0("prop_var_factor_",
                                                                    input$factor)]]
    
    # Check if pathway_to_highlight is either not yet defined, or None
    # In this case, just colour all genes
    if (is.null(input$pathway_to_highlight) || input$pathway_to_highlight == "None") {
      factor_gene_info$gene_highlighted <- TRUE
      showlegend <- FALSE
    } else {
      # Only colour genes belonging to this pathway
      pathway_name <- gsub(" \\([0-9]+ / [0-9]+ genes\\)", "", input$pathway_to_highlight)
      factor_gene_info$gene_highlighted <- factor_gene_info[[pathway_name]] == 1
      showlegend <- TRUE
    }
    
    p <- ggplot(factor_gene_info, aes(x=factor_loading,
                                      y=prop_var_explained,
                                      colour=gene_highlighted,
                                      text=paste("Gene:", GeneSymbol))) +
      scale_colour_manual(values=c("TRUE"="goldenrod", "FALSE"="black")) +
      geom_point(alpha=0.8, size=1.5) +
      labs(x="Gene's loading in factor",
           y="Proportion of variance of gene explained by factor")
    with_options(list(digits=4), # Only print 4 digits in hovertext
                 ggplotly(p, tooltip=c("text", "x", "y")) %>%
                   layout(legend=list(bgcolor = "cornsilk",
                                      bordercolor = "#FFFFFF",
                                      borderwidth = 2,
                                      title=list(text="In selected pathway"),
                                      orientation = "h",
                                      y = -0.3),
                          showlegend=showlegend))
  })
  output$factor_scatterplot <- renderPlotly({
    factor_info$highlighted <- (factor_info$factor == paste0("factor_", input$factor))
    factor_info$hovertext <- paste0(factor_info$factor,
                               "<br>Genes: ", factor_info$num_genes,
                               "<br>Samples: ", factor_info$num_samples,
                               "<br>Prop. variance explained (mean in factor): ", signif(factor_info$local_mean_prop_var_explained, 3),
                               "<br>Prop. variance explained (mean all genes): ", signif(factor_info$mean_prop_var_explained, 3))
    
    variable_names <- list("num_genes"="Number of genes",
                           "num_samples"="Number of samples",
                           "local_mean_prop_var_explained"="Proportion of variance explained (mean of genes in factor)",
                           "mean_prop_var_explained"="Proportion of variance explained (mean of all genes)")
    p <- ggplot(factor_info, aes_string(x=input$factor_scatterplot_x,
                                        y=input$factor_scatterplot_y,
                                        text="hovertext",
                                        colour="highlighted")) +
      geom_point() +
      labs(x=variable_names[[input$factor_scatterplot_x]],
           y=variable_names[[input$factor_scatterplot_y]]) +
      scale_colour_manual(values=c("TRUE"="goldenrod", "FALSE"="black"))
    
    with_options(list(digits=3), # Only print 4 digits in hovertext,
                 ggplotly(p,
                          tooltip=c("text")) %>%
                   layout(showlegend=FALSE))
  })
  
  output$pathway_dropdown <- renderUI({
    pathways <- as.list(paste0(enriched_pathways()$PathwayName,
                               " (", enriched_pathways()$GenesInIntersection,
                               " / ", enriched_pathways()$GenesInPathway, " genes)"))
    selectInput("pathway_to_highlight", "Pathway to highlight", choices=c("None", pathways))
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

