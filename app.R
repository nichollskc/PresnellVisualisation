library(shiny)
library(ggplot2)
library(plotly)
library(dplyr)
library(heatmaply)

setwd("~/docs/PhD/PresnellVisualisation/")
source("biclust_utils.R")
source("presnell.R")

SSLB_result <- read_thresholded_factor_matrices(main_run)

ui <- fluidPage(
  titlePanel("Biclustering on Presnell sorted blood cell dataset"),
  sidebarLayout(
    sidebarPanel(
      numericInput(inputId="factor", label="Factor",
                   value=11, min=1, max=SSLB_result$K, step=1),
    ),
    mainPanel(
      plotlyOutput("factorcontribution_heatmap"),
    )
  )
)

server <- function(input, output) {
  factor_contribution <- reactive({
    calc_factor_contribution(SSLB_result, input$factor)
  })
  nz_factor_contribution <- reactive({
    calc_nz_factor_contribution(SSLB_result, input$factor)
  })
  output$factorcontribution_heatmap <- renderPlotly({
    heatmaply(nz_factor_contribution())
  })
}

shinyApp(ui = ui, server = server)
