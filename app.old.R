
library(shiny)
library(ggplot2)
library(plotly)

source("ve_sim.R")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Server code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
server <- function(input,output) {
  reac <- reactiveValues()
  
  observe({reac$beta = input$beta})
  observe({reac$contactRate = input$contactRate})
  observe({reac$prev = input$prev})
  observe({reac$epsilon = input$epsilon})
  #observe({reac$size = input$size})
  observe({reac$inc = input$inc})
  observe({reac$nsteps = input$nsteps})
  observe({reac$sampleSize = input$sampleSize})
  
  #browser()
  
  output$plot1  <-  renderPlotly(runSimByPropHigh(reac))
  output$plot2  <-  renderPlotly(runSimByInc(reac))
  output$plot3  <-  renderPlotly(runSimByEpsilon(reac))
  
  
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# UI layout ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ui <- fluidPage(
 
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
  ),
  
  titlePanel(htmlTemplate("header.html")),

  sidebarLayout(
    sidebarPanel(  
      sliderInput('beta', 'beta (per contact transmission rate):', min=0, max=0.01,
                  value=0.004, step=0.001, round=-4),
      sliderInput('contactRate', 'contact rate (contacts per day):', min=0, max=1,
                  value=25/365, step=0.01, round=FALSE),
      sliderInput('prev', 'prev (prevalence):', min=0, max=1,
                  value=0.1, step=0.1, round=FALSE),
      sliderInput('epsilon', 'epsilon (per contact vaccine efficacy):', min=0, max=1,
                  value=0.3, step=0.1, round=FALSE),
      #sliderInput('size', 'sample size (population size):', min=0, max=10000,
       #           value=5000, step=500, round=FALSE),
      sliderInput('inc', 'inc (incidence, per 100 person years:', min=0, max=0.05,
                  value=0.04, step=0.001, round=-3),
      sliderInput('sampleSize', 'Sample size (population size):', min=0, max=10000,
                  value=5000, step=500, round=FALSE),
      sliderInput('nsteps', 'nsteps (days):', min=0, max=3650,
                  value=365*3, step=100, round=FALSE),
    ),
    
    mainPanel(
      
      plotlyOutput("plot1"),
      plotlyOutput("plot2"),
      plotlyOutput("plot3")
      
    ),
  ),
  titlePanel(htmlTemplate("template.html"))
             
  
)

#Run ----
shinyApp(ui = ui, server=server)
