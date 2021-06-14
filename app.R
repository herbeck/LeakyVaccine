
library(shiny)
library(ggplot2)
library(plotly)

library(deSolve)
library(tidyverse)
library(EpiModel)
library(survival)
library(EasyABC)
library(shinythemes)

source("ve_sim.R")
source("tabContent.R")
source("sim_fns.R")
source("sim_plots.R")


server <- function(input, output, session) {
  
  updateTabsetPanel(session, "page-nav", "Introduction")
  #--------------
  # for old model
  #--------------
  reacOld <- reactiveValues()
  
  observe({reacOld$beta = input$betaOld})
  observe({reacOld$contactRate = input$contactRateOld})
  observe({reacOld$prev = input$prevOld})
  observe({reacOld$epsilon = input$epsilonOld})
  observe({reacOld$size = input$sizeOld})
  observe({reacOld$inc = input$incOld})
  observe({reacOld$nsteps = input$nstepsOld})
  observe({reacOld$sampleSize = input$sampleSizeOld})
  
  output$plotOld1  <-  renderPlotly(runSimByPropHigh(reacOld))
  output$plotOld2  <-  renderPlotly(runSimByInc(reacOld))
  output$plotOld3  <-  renderPlotly(runSimByEpsilon(reacOld))
  
  #------------------
  # end for old model
  #------------------
  
  reac <- reactiveValues()
  
  observe({reac$beta = input$beta})
  observe({reac$contactRate = input$contactRate})
  observe({reac$prev = input$prev})
  observe({reac$epsilon = input$epsilon})
  observe({reac$risk = input$risk}) 
  

  createCInfectionPlot(output, reac)
  createPlaceboRiskPlot(output,reac)
  createCIncidencePlot(output, reac)
  createPlaceboVaccinePlot(output,reac)
  createPlaceboVaccineRiskPlot(output,reac)

}

#------------------------------------------------------------------------------
# for creating UI
#------------------------------------------------------------------------------
ui <- navbarPage(
  
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "stylesContent.css"),
  ),
  
  title ="LeakyVaccine",
  id = "page-nav",
  theme = shinytheme("cerulean"),
  
  #tabs
  getIntroContent(),
  getModelSetupContent(),
  getCalibrationContent(),
  getJoshPlotContent(),
  getParameterSweepContent()
  #titlePanel(htmlTemplate("template.html"))

)


#Run ----
shinyApp(ui = ui, server=server)