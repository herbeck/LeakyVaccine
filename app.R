
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

#------------------------------------------------------------------------------
# Right now this includes just two vaccine trial populations, each with a 
# vaccine arm and a placebo arm. One population has homogeneous exposure / risk
# of infection; the other population includes exposure heterogeneity, and this 
# heterogeneity is the same in both trial arms.
#------------------------------------------------------------------------------
si_ode <- function(times, init, param){
  with(as.list(c(init, param)), {
    
    # Flows
    # the number of people moving from S to I at each time step
    #Susceptible, Infected, placebo
    SIp.flow <- lambda*Sp
    SIv.flow <- lambda*(1-epsilon)*Sv
    
    #Susceptible, Infected, placebo, high, medium, low
    SIph.flow <- risk*lambda*Sph
    SIpm.flow <- lambda*Spl
    #SIpl.flow <- 0*lambda*Spl  #zero risk in the low risk group
    SIpl.flow <- (1/risk)*lambda*Spl  #inverse risk multiplier for low risk group
    
    #Susceptible, Infected, vaccine, high, medium, low
    SIvh.flow <- risk*lambda*(1-epsilon)*Svh
    SIvm.flow <- lambda*(1-epsilon)*Spl
    #SIvl.flow <- 0*lambda*(1-epsilon)*Svl  #zero risk in the low risk group 
    SIvl.flow <- (1/risk)*lambda*(1-epsilon)*Svl  #inverse risk 
    
    # ODEs
    # placebo; homogeneous risk
    dSp <- -SIp.flow
    dIp <- SIp.flow  #lambda*Sp
    
    # vaccine; homogeneous risk
    dSv <- -SIv.flow
    dIv <- SIv.flow  #lambda*epsilon*Sv
    
    # placebo; heterogeneous risk
    dSph <- -SIph.flow
    dIph <- SIph.flow  #risk*lambda*Sph
    dSpm <- -SIpm.flow
    dIpm <- SIpm.flow  #lambda*Spm
    dSpl <- -SIpl.flow
    dIpl <- SIpl.flow  #0*lambda*Spl
    
    # vaccine; heterogeneous risk
    dSvh <- -SIvh.flow
    dIvh <- SIvh.flow  #risk*lambda*(1-epsilon)*Svh
    dSvm <- -SIvm.flow
    dIvm <- SIvm.flow  #lambda*Svm
    dSvl <- -SIvl.flow
    dIvl <- SIvl.flow  #0*lambda*(1-epsilon)*Svl
    
    #Output
    list(c(dSp,dIp,
           dSv,dIv,
           dSph,dIph,
           dSpm,dIpm,
           dSpl,dIpl,
           dSvh,dIvh,
           dSvm,dIvm,
           dSvl,dIvl,
           SIp.flow,SIv.flow,
           SIph.flow,SIpm.flow,SIpl.flow,
           SIvh.flow,SIvm.flow,SIvl.flow))
  })
}

#------------------------------------------------------------------------------
# This function just takes the model output (a `mod` file) and uses the data to 
# create other data (e.g. incidence and VE estimates) for plotting.
#------------------------------------------------------------------------------
mod.manipulate <- function(mod){
  
  mod <- mutate_epi(mod, total.Svh.Svm.Svl = Svh + Svm + Svl) #all susceptible in heterogeneous risk vaccine pop
  mod <- mutate_epi(mod, total.Sph.Spm.Spl = Sph + Spm + Spl) #all susceptible in heterogeneous risk placebo pop
  mod <- mutate_epi(mod, total.Ivh.Ivm.Ivl = Ivh + Ivm + Ivl) #all infected in heterogeneous risk vaccine pop
  mod <- mutate_epi(mod, total.Iph.Ipm.Ipl = Iph + Ipm + Ipl) #all infected in heterogeneous risk placebo pop
  mod <- mutate_epi(mod, total.SIvh.SIvm.SIvl.flow = SIvh.flow + SIvm.flow + SIvl.flow) #all infections per day in heterogeneous risk vaccine pop
  mod <- mutate_epi(mod, total.SIph.SIpm.SIpl.flow = SIph.flow + SIpm.flow + SIpl.flow) #all infections in heterogeneous risk placebo pop
  
  #Instantaneous ncidence (hazard) estimates, per 100 person years
  #Instantaneous incidence / hazard
  mod <- mutate_epi(mod, rate.Vaccine = (SIv.flow/Sv)*365*100)
  mod <- mutate_epi(mod, rate.Placebo = (SIp.flow/Sp)*365*100)
  mod <- mutate_epi(mod, rate.Vaccine.het = (total.SIvh.SIvm.SIvl.flow/total.Svh.Svm.Svl)*365*100)
  mod <- mutate_epi(mod, rate.Placebo.het = (total.SIph.SIpm.SIpl.flow/total.Sph.Spm.Spl)*365*100)
  
  #Cumulative incidence
  mod <- mutate_epi(mod, cumul.Sv = cumsum(Sv))
  mod <- mutate_epi(mod, cumul.Sp = cumsum(Sp))
  mod <- mutate_epi(mod, cumul.Svh.Svm.Svl = cumsum(total.Svh.Svm.Svl))
  mod <- mutate_epi(mod, cumul.Sph.Spm.Spl = cumsum(total.Sph.Spm.Spl))
  mod <- mutate_epi(mod, cumul.rate.Vaccine = (Iv/cumul.Sv)*365*100)
  mod <- mutate_epi(mod, cumul.rate.Placebo = (Ip/cumul.Sp)*365*100)
  mod <- mutate_epi(mod, cumul.rate.Vaccine.het = (total.Ivh.Ivm.Ivl/cumul.Svh.Svm.Svl)*365*100)
  mod <- mutate_epi(mod, cumul.rate.Placebo.het = (total.Iph.Ipm.Ipl/cumul.Sph.Spm.Spl)*365*100)
  
  #Vaccine efficacy (VE) estimates
  #VE <- 1 - Relative Risk; this is VE for hazard
  mod <- mutate_epi(mod, VE1.inst = 1 - rate.Vaccine/rate.Placebo)
  mod <- mutate_epi(mod, VE2.inst = 1 - rate.Vaccine.het/rate.Placebo.het)
  
  #VE <- 1 - Relative Risk; this is VE from cumulative incidence
  mod <- mutate_epi(mod, VE1.cumul = 1 - cumul.rate.Vaccine/cumul.rate.Placebo)
  mod <- mutate_epi(mod, VE2.cumul = 1 - cumul.rate.Vaccine.het/cumul.rate.Placebo.het)
  
  return(mod)
}

#------------------------------------------------------------------------------
# Initial parameter settings
#------------------------------------------------------------------------------
# beta <- 0.004   #transmission rate (per contact)
# c <- 90/365    #contact rate (contacts per day)
# prev <- 0.10    #needs some consideration
# #prev <- 0.01
# lambda <- beta*c*prev
# #lambda <- 0.000008398975
# epsilon <- 0.30  #per contact vaccine efficacy
# risk <- 10.0   #risk multiplier


runSim <- function(param) {
  
  param <- param.dcm(lambda = param$beta * param$contactRate * param$prev, 
                     epsilon = param$epsilon, 
                     risk = param$risk)
  init <- init.dcm(Sp = 10000, Ip = 0,
                   Sv = 10000, Iv = 0,
                   Sph = 1000, Iph = 0,    #placebo, high risk
                   Spm = 7000, Ipm = 0,    #placebo, medium risk
                   Spl = 2500, Ipl = 0,    #placebo, low risk
                   Svh = 1000, Ivh = 0,    #vaccine, high risk
                   Svm = 7000, Ivm = 0,    #vaccine, medium risk
                   Svl = 2500, Ivl = 0,    #vaccine, low risk
                   SIp.flow = 0, SIv.flow = 0, 
                   SIph.flow = 0, SIpm.flow = 0, SIpl.flow = 0,
                   SIvh.flow = 0, SIvm.flow = 0, SIvl.flow = 0)
  
  control <- control.dcm(nsteps = 365*3, new.mod = si_ode)
  mod <- dcm(param, init, control)
  #mod
  
  mod <- mod.manipulate(mod)
  return (mod)
}

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
  #observe({reacOld$size = input$sizeOld})
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
  

  output$plot1 <- renderPlot({
    
      mod <- runSim(reac)
    
      plot(mod, y = c("Ip", "total.Iph.Ipm.Ipl"), 
       alpha = 0.8, 
       main = "Cumulative infections",
       legend = FALSE,
       ylab = "infected",
       xlab = "days",
       col = c("blue", "red"))
    
  legend("bottomright", legend = c("homogeneous risk", "heterogeneous risk"), 
         col = c("blue", "red"), lwd = 2, cex = 0.9, bty = "n")
  })
  
  output$plot2 <- renderPlot({
    
    mod <- runSim(reac)
    
    plot(mod, y=c("rate.Placebo", "rate.Placebo.het"),
         alpha = 0.8,
         ylim = c(0, 4.5),
         main = "Hazard",
         xlab = "days",
         ylab = "infections per 100 person yrs",
         legend = FALSE,
         col = c("blue", "red"))
    legend("bottomright", legend = c("homogen. risk", "heterogen. risk"), col = c("blue", "red"), lwd = 2, cex = 0.9, bty = "n")
    
    par(mfrow=c(1,1))
  })
  
  #Vaccine now in the heterogeneous risk population, added to the above plot
  
  output$plot3 <- renderPlot({
    
    mod <- runSim(reac)
    
    #mod <- mod.manipulate(mod)
      
    plot(mod, y=c("cumul.rate.Placebo", "cumul.rate.Vaccine"),
         alpha = 0.8,
         ylim = c(0, 4.5),
         main = "Cumulative incidence",
         xlab = "days",
         ylab = "infections per 100 person yrs",
         legend = FALSE,
         col = c("blue", "green"))
    legend("bottomright", legend = c("placebo", "vaccine"), col = c("blue", "green"), lwd = 2, cex = 0.9, bty = "n")
  })
  
  #Instantaneous incidence / hazard, for the same comparison
  
  output$plot4 <- renderPlot({
    mod <- runSim(reac)

    plot(mod, y=c("rate.Placebo", "rate.Vaccine"),
         alpha = 0.8,
         ylim = c(0, 4.5),
         main = "Hazard",
         xlab = "days",
         ylab = "infections per 100 person yrs",
         legend = FALSE,
         col = c("blue", "green"))
    legend("bottomright", legend = c("placebo", "vaccine"), col = c("blue", "green"), lwd = 2, cex = 0.9, bty = "n")
  })
  
  #Vaccine now in the heterogeneous risk
  
  output$plot5 <- renderPlot({
    mod <- runSim(reac)
    plot(mod, y=c("rate.Placebo", "rate.Vaccine", "rate.Placebo.het", "rate.Vaccine.het"),
         alpha = 0.8,
         ylim = c(0, 4.5),
         main = "Hazard",
         xlab = "days",
         ylab = "infections per 100 person yrs",
         legend = FALSE,
         col = c("blue", "green", "red", "orange"))
    legend("bottomright", legend = c("placebo, homogeneous risk", "vaccine, homogeneous risk", 
          "placebo, heterogeneous risk", "vaccine, heterogeneous risk"), 
            col = c("blue", "green", "red", "orange"), lwd = 2, cex = 0.9, bty = "n")
  })
  
}

ui <- navbarPage(

        tags$head(
          tags$link(rel = "stylesheet", type = "text/css", href = "stylesContent.css"),
        ),
        title ="Leaky vaccines and exposure heterogeneity",
        id = "page-nav",
        theme = shinytheme("cerulean"),
        
        tabPanel("Introduction",
          HTML("<div class='content'>"),
          h3("Using models to assess the impact of HIV exposure heterogeneity on trial vaccine efficacy measures"),
          p("It is hypothesized that exposure heterogeneity (i.e. infection risk heterogeneity) can affect efficacy estimation for leaky vaccines (e.g. Halloran et al., 1992; White et al., 2010; O'Hagan et al.,2013; Edlefsen, 2014; Coley et al., 2016; Gomes et al., 2016; Kahn et al., 2018). Our goal is to make a simple deterministic compartmental model to facilitate straightforward simulation-based evaluation of this process within and across populations, in the context of HIV prevention trials or longitudinal studies."),
          HTML("<ol type='1'>"),
          HTML("<li>In acute infection studies it seems like many participants get infected early. What is the magnitude of this effect that might be due to frailty bias?</li>"),
          HTML("<li>Assess if this effect might contribute to the differences between the RV 144 and HVTN 702 vaccine trial outcomes. (There have been a couple of analyses of this, and we can build on this and make future analyses of other trial results more straightforward to evaluate.)</li>"),
          HTML("<li>Assess if this effect might contribute to the waning efficacies seen in HIV prevention trials (specifically the AMP VRC01 bnAb trial).</li>"),
          HTML("<li>In the context of the AMP Trial and the different results seen in the sub-studies (703 vs 704); is this due to different forces of infection between the populations?</li>"),
          HTML("<li>Continue to raise awareness of this issue to HIV prevention trials, with the ultimate goal of better design and interpretation of efficacy outcomes.</li>"),
          HTML("</ol>"),
          p("From Gomes et al., 2016:  \"This effect is more pronounced in the control group as individuals within it experience higher rates of infection overall. Consequently, the ratio of disease rates in vaccinated over control groups increases, and vaccine efficacy, as measured by simple rate ratios, decreases as the trial progresses. Finally, the magnitude of this effect increases with the intensity of transmission.\"  "),
          HTML("</div>"),
        ),
        tabPanel("Model setup",
          HTML("<div class='content'>"),
          p(paste("We are modeling a vaccine trial using an SI deterministic compartmental model.", 
                  "We are not modeling infections from the I to S compartments but rather only infections from the outside (non-trial) population, ", 
                  "with the infection rate based on the population prevalence `prev` (of viremic individuals), the exposure rate ",
                  "(serodiscordant sexual contacts per time) `c`, and the transmission rate (per contact) `p`. The per contact ",
                  "effect of vaccination is `epsilon`, and with this iteration of the model `epsilon` is: ")),
          HTML("<ol type='1'><li> not time-varying (the per contact vaccine effect does not decay over time) and</li>  <li>assumes a homogeneous effect (does not vary by mark / viral genotype).</li></ol>"),
          p("This model structure also removes the possibility of indirect effects from vaccination.  "),
          p("We are, with this early iteration, including just three subgroups in the heterogeneous exposure population: high, medium, and low exposure. Right now we do not know the correct size of these subgroups (i.e. fraction of the population) or their relative contribution to overall incidence. First pass is 10% high risk, 80% medium risk, 10% low risk (and low risk is set at zero risk), as this, in combination with a 10% `risk` multiplier, results in 3.5% incidence in a putative population with no vaccine (placebo arm) and no exposure heterogeneity. "),
          
          HTML(paste("<div class='code'>", 
                     "<div class='flex'><div class='definition'>beta</div><div>transmission rate (per contact)</div></div>",
                     "<div class='flex'><div class='definition'>c</div><div>exposure rate (serodiscordant sexual contacts per time)</div></div>",
                     "<div class='flex'><div class='definition'>prev</div><div>prevalence  (prevalence of viremic individuals)</div></div>",
                     "<div class='flex'><div class='definition'>lambda</div><div>lambda = beta * c * prev</div></div>",
                     "<div class='flex'><div class='definition'>risk</div><div>risk multiplier</div></div>",
                     "<div class='flex'><div class='definition'>epsilon</div><div>per contact vaccine efficacy; vaccine-induced reduction in the risk of HIV infection from a single exposure</div></div>",
                     "</div><br/>")),
          HTML("<p>The risk multiplier is an amalgam of increases in transmission risk that could be due to:</p>"),
          HTML("<ol type='1'>"),
          HTML("<li>increased per contact transmission risk</li>"),
          HTML("<li>increased exposure rate (number of contacts)</li>"),
          HTML("<li>increased prevalence of HIV viremia in partners.</li></ol>"),
          HTML("<p>Individual risk of infection can vary for these separately or in combination.</p>"),
          h4("Key for our model functions:"),
          HTML(paste("<div class='code'>", 
                     "<div class='flex'><div class='definition'>sp</div><div>susceptible placebo</div></div>",
                     "<div class='flex'><div class='definition'>Ip</div><div>infected placebo</div></div>",
                     "<div class='flex'><div class='definition'>Sv</div><div>susceptible vaccinated</div></div>",
                     "<div class='flex'><div class='definition'>Iv</div><div>infected vaccinated</div></div>",
                     "<div class='flex'><div class='definition'>Svh</div><div>susceptible vaccinated high exposure</div></div>",
                     "<div class='flex'><div class='definition'>Svm</div><div>susceptible vaccinated medium exposure</div></div>",
                     "<div class='flex'><div class='definition'>SvL</div><div>susceptible vaccinated low exposure (zero in this instance)</div></div>",
                     "<div class='flex'><div class='definition'>Ivh</div><div>infected vaccinated high exposure</div></div>",
                     "<div class='flex'><div class='definition'>Ivm</div><div>infected vaccinated medium exposure</div></div>",
                     "<div class='flex'><div class='definition'>Ivl</div><div>infected vaccinated low exposure (zero in this instance)</div></div>",
                     "</div><br/>")),
          HTML("</div>"),
          
   
        ),
        tabPanel("Calibration",
          HTML("<div class='content'>"),
          p("We use calibration in the following steps: 1. calibrate placebo incidence to RV144 (0.035%) At this stage we just eyeball-calibrated the incidence to ~3.5% per 100 person years, to be reasonably consistent with HVTN 702 in South Africa."),
          p("We used an initial set of transmission parameters for sub-Saharan Africa borrowing from an SIR model from Alain Vandormael (2018): 'We used realistic parameter values for the SIR model, based on earlier HIV studies that have been undertaken in the sub-Saharan Africa context. To this extent, we varied `c` within the range of 50 to 120 sexual acts per year based on data collected from serodiscordant couples across eastern and southern African sites. Previous research has shown considerable heterogeneity in the probability of HIV transmission per sexual contact, largely due to factors associated with the viral load level, genital ulcer disease, stage of HIV progression, condom use, circumcision and use of ART. Following a systematic review of this topic by Boily et al., we selected values for `beta` within the range of 0.003–0.008. ... Here, we chose values for `v` within the range of 0.15–0.35, which are slightly conservative, but supported by population-based estimates from the sub-Saharan African context."),
          HTML(paste("<div class='code'>", 
                     "<div class='flex'><div class='definition'>c</div><div>varies from 50 to 120 per year</div></div>",
                     "<div class='flex'><div class='definition'>beta</div><div>varies from 0.003 to 0.008</div></div>",
                     "<div class='flex'><div class='definition'>prev</div><div>which here is population prevalence of unsuppressed VL, varies from 0.15 to 0.35</div></div>",
                     "<div class='flex'><div class='definition'>Sv</div><div>could be parameterized using the RV144 Thai Trial results: VE = 61% at 12 months, 31% at 42 months, but below we start with 30% and no waning efficacy. A vaccine duration parameter is not needed because we are only modeling a 3 year trial without boosters.</div></div>",
                     "</div>")),
          HTML("</div>")
        ),
        
        tabPanel("Simple plots", 
                 sidebarLayout(
                   sidebarPanel(  
                     sliderInput('beta', 'Beta:', min=0, max=0.01,
                                 value=0.004, step=0.001, round=-4),
                     sliderInput('contactRate', 'Contacts rate:', min=0, max=1,
                                 value=25/365, step=0.01, round=FALSE),
                     sliderInput('prev', 'prev:', min=0, max=1,
                                 value=0.1, step=0.1, round=FALSE),
                     sliderInput('epsilon', 'epsilon:', min=0, max=1,
                                 value=0.3, step=0.1, round=FALSE),
                     sliderInput('risk', 'risk:', min=0, max=20,
                                 value=10, step=1, round=FALSE)
                     
                   ),
                   mainPanel(
                      fluidRow(
                        column(6, plotOutput("plot1")),
                        column(6, plotOutput("plot2"))),
                      fluidRow(
                        column(6, plotOutput("plot3")),
                        column(6, plotOutput("plot4"))),
                      plotOutput("plot5")
                             )
                 )
        ),

        tabPanel("Parameter sweeps", 
                  sidebarLayout(
                   sidebarPanel(  
                     sliderInput('betaOld', 'beta (per contact transmission rate):', min=0, max=0.01,
                                 value=0.004, step=0.001, round=-4),
                     sliderInput('contactRateOld', 'Contact rate:', min=0, max=1,
                                 value=25/365, step=0.01, round=FALSE),
                     sliderInput('prevOld', 'Prevalence of viremia:', min=0, max=1,
                                 value=0.1, step=0.1, round=FALSE),
                     sliderInput('epsilonOld', 'Epsilon (per-contact vaccine efficacy):', min=0, max=1,
                                 value=0.3, step=0.1, round=FALSE),
                     #sliderInput('sizeOld', 'sample size:', min=0, max=10000,
                                 #value=5000, step=500, round=FALSE),
                     sliderInput('incOld', 'inc (cumulative incidence, per 100 person years):', min=0, max=1,
                                 value=0.04, step=0.01, round=-3),
                     sliderInput('sampleSizeOld', 'sample size (population size):', min=0, max=10000,
                                 value=5000, step=500, round=FALSE),
                     sliderInput('nstepsOld', 'nsteps:', min=0, max=3650,
                                 value=365*3, step=100, round=FALSE),
                   ),
                   
                   mainPanel(
                     
                     plotlyOutput("plotOld1"),
                     plotlyOutput("plotOld2"),
                     plotlyOutput("plotOld3")
                     
                   )
                 )
                        
        )
  #titlePanel(htmlTemplate("template.html"))
  
)


#Run ----
shinyApp(ui = ui, server=server)