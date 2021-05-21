library(shiny)
library(plotly)
library(heatmaply)
library(shinydashboard)
library(shinyjs)

numericSpinnerInput <- function(inputId, value = 0) {
  tagList(
    singleton(tags$head(tags$script(src = "js/numeric_spinner.js"))),
    tags$input(id = inputId,
               class = "spinner",
               value = as.character(value))
  )
}

#################################################################################################
# Defining plot HTML objects                                                                    #
#################################################################################################
sample_heatmap_output <- box(title="Sample types",
                             div(paste0("Types of samples included in the factor, broken down by disease, ",
                                        "cell type and sex. The numbers in the cells give the percentage of samples",
                                        " of that type which are contained in this factor. The darker the colour ", 
                                        "(red for sex, or blue for cell type), the higher this percentage.",
                                        "Grey represents sample types not present in the dataset.")),
                             plotlyOutput("sample_heatmap"),
                             collapsible=TRUE, solidHeader = TRUE, status="primary", width=NULL,
)

pathways_table_output <- box(title="Enriched pathways",
                             div(paste0("Pathways with strongest enrichment (measured by the q-value returned by ",
                                        "CameraPR) in this factor.")),
                             div(style = 'overflow-x: scroll', dataTableOutput("pathways_table")),
                             collapsible=TRUE, solidHeader = TRUE, status="primary", width=NULL,
)

factorcontribution_output <- box(title="Factor contribution",
                                 id="factorcontribution_box",
                                 div(paste0("Heatmap showing the values this factor contributes to the overall ",
                                            "reconstruction of the original matrix. Constructed by $x_k b_k^T$ ",
                                            "where $x_k$ is the length n vector giving sample loadings and ",
                                            "$b_k$ is the length p vector giving gene loadings.")),
                                 plotlyOutput("factorcontribution_heatmap"),
                                 selectInput("factorcontribution_trans", "Transformation for colour scale:",
                                             c("No transformation", "Log transformation"), "No transformation"),
                                 selectInput("factorcontribution_scale", "Max/min defined by:",
                                             c("This factor", "All factors"), "This factor"),
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

#################################################################################################
# UI - Specifying HTML layout                                                                   #
#################################################################################################
ui <- dashboardPage(
  dashboardHeader(title="Biclustering on Presnell sorted blood cell dataset"),
  dashboardSidebar(
    div(
      tags$label("Factor:"),
      numericSpinnerInput("factor", 12),
    ),
    div(
      htmlOutput("summary_text")
    )
  ),
  dashboardBody(
    # Allow custom JS script to collapse box when you click on header
    useShinyjs(),
    
    # Include jquery core
    tags$head(
#        tags$script(src="js/jquery-ui-1.12.1.custom/external/jquery/jquery.js"),
        tags$script(src="js/jquery-ui-1.12.1.custom/jquery-ui.js"),
        tags$link(rel = "stylesheet", type = "text/css", href = "js/jquery-ui-1.12.1.custom/jquery-ui.css"),
        tags$link(rel = "stylesheet", type = "text/css", href = "css/biclust.css"),
    ),
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

ui

