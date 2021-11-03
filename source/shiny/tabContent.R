
#------------------------------------------------------------------------------
# for creating Intro tab content
#------------------------------------------------------------------------------
getAboutContent <- function() {
  return (tabPanel("About this tool",
                   
    HTML("<div class='mainPanel main'>"),
    h3("Using models to assess the impact of HIV exposure heterogeneity on trial vaccine efficacy measures"),
    p("It is hypothesized that exposure heterogeneity (i.e. infection risk heterogeneity) can affect estimates of vaccine efficacy for leaky vaccines (e.g. Halloran et al., 1992; White et al., 2010; O'Hagan et al.,2013; Edlefsen, 2014; Coley et al., 2016; Gomes et al., 2016; Kahn et al., 2018."), 
    p("A potential outcome of exposure heterogeneity (i.e. risk structure in a population) is that vaccine efficacy measured from a trial (i.e. the clinical efficacy) is lower than the biological vaccine efficacy (i.e. the per-exposure or per-contact vaccine efficacy). This distinction is important: the per-exposure vaccine efficacy is not necessarily equal to the clinical efficacy or the population effectiveness of the same vaccine."),
    p("From Gomes et al., 2016:  \"This effect is more pronounced in the control group as individuals within it experience higher rates of infection overall. Consequently, the ratio of disease rates in vaccinated over control groups increases, and vaccine efficacy, as measured by simple rate ratios, decreases as the trial progresses. Finally, the magnitude of this effect increases with the intensity of transmission.\"  "),
    p("Here we use epidemic models to simulate this process, within and across populations, in the context of HIV prevention trials or longitudinal studies. Our goals are to:"),
    HTML("<ol type='1'>"),
    HTML("<li>Raise awareness of the distinction between per-exposure vaccine efficacy, clinical vaccine efficacy, and population vaccine effectiveness.</li>"),
    HTML("<li>Assess if this effect might contribute to the difference between the RV144 and HVTN 702 vaccine trial outcomes.</li>"),
    HTML("<li>Assess if this effect might contribute to the waning efficacies seen in HIV prevention trials (for example the AMP VRC01 bnAb trials).</li>"),
    HTML("<li>In acute infection studies it seems like many participants get infected early. What is the magnitude of this effect that might be due to frailty bias?</li>"),
    HTML("<li>Assist in the design or interpretation of HIV prevention trials, from this exposure heterogeneity framework.</li>"),
    HTML("</ol>"),
    p("Thes separate tabs in this R Shiny app include:  
    Model Description; showing the structure of the model and the parameters included
    Initial Example Plots; showing how the model works and what simulated epidemic and trial outputs we focus on
    Parameter Sweeps; which allows the user to compare the impact of multiple parameter values in the same plots
    Calibration; which describes our approach at using this model for HIV populations, and
    Model Fitting; in which we use the model to examine specific trial results."),
    HTML("</div>"),
    titlePanel(htmlTemplate("template.html"))
  ))
}


#------------------------------------------------------------------------------
# for creating Model Description tab content
#------------------------------------------------------------------------------
getModelDescriptionContent <- function() {
  return(tabPanel("Model description",

    HTML("<div class='mainPanel main'>"),
    h3("Model structure"),
    p("We are modeling a vaccine trial using an SI deterministic compartmental model; this is a simple epidemic model that has two populations, the Susceptible (S) and the Infected (I). 
    We start the epidemic model with all trial individuals in the S group; over time S individuals move into the I group as they become infected.", 
    p("We are not modeling infections from the I to S compartments but rather only infections from the outside (non-trial) population."),  
    p("This model structure also removes the possibility of indirect effects from vaccination.  "),
    h3("Parameters"),
    p("The infection rate per time step is a combination of population prevalence `prev` (of viremic individuals), the exposure rate
    (serodiscordant sexual contacts per time) `c`, and the transmission rate (per contact) `p`. The per contact
    effect of vaccination is `epsilon`."),
    HTML("<ol type='1'><li> not time-varying (the per contact vaccine effect does not decay over time) and</li>  <li>assumes a homogeneous effect (does not vary by mark / viral genotype).</li></ol>"),
    p("We are including three subgroups in the heterogeneous exposure population: high, medium, and low exposure. 
    In reality we do not know the correct size of risk subgroups (i.e. fraction of the population) or their relative contribution to overall incidence. 
    But for our first pass we set the size of these groups at 10% high risk, 80% medium risk, 10% low risk (and low risk is set at zero risk)."),
           
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
     HTML("<li>higher per contact transmission risk</li>"),
     HTML("<li>higher exposure rate (number of contacts)</li>"),
     HTML("<li>higher prevalence of HIV viremia in partners.</li></ol>"),
     HTML("<p>Individual risk of infection can vary for these separately or in combination.</p>"),
     titlePanel(htmlTemplate("template.html"))
           
  ))
  
}


#------------------------------------------------------------------------------
# for creating Initial Example Plots tab content
#------------------------------------------------------------------------------

getInitialExamplePlotsContent <- function() {
  tabPanel("Initial example plots", 
           HTML("<div class='mainPanel'>"),
             sidebarPanel(  
               sliderInput('beta', 'beta (per-contact transmission probability):', min=0, max=0.01,
                           value=0.005, step=0.001, round=-4),
               sliderInput('contactRate', 'c (sexual contacts per day):', min=0, max=1,
                           value=90/365, step=0.01, round=FALSE),
               sliderInput('prev', 'prev (population prevalence of viremic individuals):', min=0, max=1,
                           value=0.10, step=0.1, round=FALSE),
               sliderInput('epsilon', 'epsilon (per-exposure vaccine efficacy):', min=0, max=1,
                           value=0.5, step=0.1, round=FALSE),
               sliderInput('risk', 'risk (risk multiplier; relative force of infection for high risk group):', min=0, max=25,
                           value=15, step=1, round=FALSE) 
             ),
             mainPanel(
              plotOutput("CumulativeInfectionsPlot") %>% withSpinner(color="#0dc5c1"),
              p("Figure 1. Cumulative infections in the placebo arms of a vaccine trial, for populations with homogeneous risk and heterogeneous risk. Note that the infections in the heterogeneous risk population accumulate faster early in the trial, as the high risk individuals are infected."),
              plotOutput("PlaceboRiskPlot") %>% withSpinner(color="#0dc5c1"),
              p("Figure 2. Incidence in the placebo arm of a vaccine trial. As expected from the cumulative infections plot above, the incidence in the heterogeneous risk population decreases over the course of the trial."),
              plotOutput("PlaceboVaccineRiskPlot") %>% withSpinner(color="#0dc5c1"),
              p("Figure 3. Incidence in the placebo and vaccine arms of a trial, for populations with homogeneous risk and heterogeneous risk."),
              plotOutput("VEPlot") %>% withSpinner(color="#0dc5c1"),
              p("Figure 4. Clinical vaccine efficacy over time in two vaccine trial population trials, one with homogeneous risk and the other with heterogeneous risk."),
              class = "plotPanel"
             ),
           HTML("</div>"),
           titlePanel(htmlTemplate("template.html"))
           
  )
}



#------------------------------------------------------------------------------
# for creating Parameter Sweeps content
#------------------------------------------------------------------------------
getParameterSweepContent <- function() {
  tabPanel("Parameter sweeps", 
           HTML("<div class='mainPanel'>"),
             sidebarPanel(  
               sliderInput('betaOld', 'beta (per-contact transmission probability):', min=0, max=0.01,
                           value=0.004, step=0.001, round=-4),
               sliderInput('contactRateOld', 'c (sexual contacts per day):', min=0, max=1,
                           value=90/365, step=0.01, round=FALSE),
               sliderInput('prevOld', 'prev (population prevalence of viremic individuals):', min=0, max=1,
                           value=0.1, step=0.1, round=FALSE),
               sliderInput('epsilonOld', 'epsilon (per-exposure vaccine efficacy):', min=0, max=1,
                           value=0.3, step=0.1, round=FALSE),
               sliderInput('incOld', 'inc (cumulative incidence, per 100 person years):', min=0, max=1,
                           value=0.04, step=0.01, round=-3),
               sliderInput('sampleSizeOld', 'sample size (population size):', min=0, max=10000,
                           value=5000, step=500, round=FALSE),
               sliderInput('nstepsOld', 'nsteps (number of times steps):', min=0, max=3650,
                           value=365*3, step=100, round=FALSE),
               class = "slider"
             ),
             
             mainPanel(
               p("Lorem ipsum dolor sit amet, consectetur adipiscing elit. Pellentesque fringilla aliquam ex, vitae scelerisque felis semper quis. Aenean nec pharetra ligula. Mauris vulputate purus ante, id faucibus leo facilisis sit amet. Fusce vestibulum justo eu enim commodo consectetur. Fusce nisi urna, ultrices at purus at, imperdiet efficitur ligula. Curabitur quis sapien ligula. Ut non orci ullamcorper, pulvinar nibh vel, molestie velit. Morbi vulputate hendrerit mi, a mollis risus blandit eget. Cras lacinia eget massa condimentum finibus. Morbi porta lorem augue, in sagittis orci vehicula vel. Sed ipsum nisi, scelerisque quis luctus at, efficitur sed erat. Nam aliquet hendrerit laoreet. Nulla rutrum, nisi pulvinar placerat eleifend, lacus metus ornare dolor, at iaculis augue mi sit amet lorem. Aliquam sit amet turpis nec quam aliquet pretium. Maecenas leo lectus, efficitur in magna ut, gravida iaculis felis."),
               plotlyOutput("plotOld1") %>% withSpinner(color="#0dc5c1"),
               p("Lorem ipsum dolor sit amet, consectetur adipiscing elit. Pellentesque fringilla aliquam ex, vitae scelerisque felis semper quis. Aenean nec pharetra ligula. Mauris vulputate purus ante, id faucibus leo facilisis sit amet. Fusce vestibulum justo eu enim commodo consectetur. Fusce nisi urna, ultrices at purus at, imperdiet efficitur ligula. Curabitur quis sapien ligula. Ut non orci ullamcorper, pulvinar nibh vel, molestie velit. Morbi vulputate hendrerit mi, a mollis risus blandit eget. Cras lacinia eget massa condimentum finibus. Morbi porta lorem augue, in sagittis orci vehicula vel. Sed ipsum nisi, scelerisque quis luctus at, efficitur sed erat. Nam aliquet hendrerit laoreet. Nulla rutrum, nisi pulvinar placerat eleifend, lacus metus ornare dolor, at iaculis augue mi sit amet lorem. Aliquam sit amet turpis nec quam aliquet pretium. Maecenas leo lectus, efficitur in magna ut, gravida iaculis felis."),
               plotlyOutput("plotOld2") %>% withSpinner(color="#0dc5c1"),
               p("Lorem ipsum dolor sit amet, consectetur adipiscing elit. Pellentesque fringilla aliquam ex, vitae scelerisque felis semper quis. Aenean nec pharetra ligula. Mauris vulputate purus ante, id faucibus leo facilisis sit amet. Fusce vestibulum justo eu enim commodo consectetur. Fusce nisi urna, ultrices at purus at, imperdiet efficitur ligula. Curabitur quis sapien ligula. Ut non orci ullamcorper, pulvinar nibh vel, molestie velit. Morbi vulputate hendrerit mi, a mollis risus blandit eget. Cras lacinia eget massa condimentum finibus. Morbi porta lorem augue, in sagittis orci vehicula vel. Sed ipsum nisi, scelerisque quis luctus at, efficitur sed erat. Nam aliquet hendrerit laoreet. Nulla rutrum, nisi pulvinar placerat eleifend, lacus metus ornare dolor, at iaculis augue mi sit amet lorem. Aliquam sit amet turpis nec quam aliquet pretium. Maecenas leo lectus, efficitur in magna ut, gravida iaculis felis."),
               plotlyOutput("plotOld3") %>% withSpinner(color="#0dc5c1"),
               class = "plotPanel"
             ),
           HTML("</div>"),
           titlePanel(htmlTemplate("template.html"))
           
  )
}

#------------------------------------------------------------------------------
# for creating Calibration tab content
#------------------------------------------------------------------------------
getCalibrationContent <- function() {
  tabPanel("Calibration",
           HTML("<div class='mainPanel main'>"),
           p("We use calibration to find model parameter settings that produce model outputs that are consistent with some target values."),
           p("We used an initial set of transmission parameters for sub-Saharan Africa borrowing from an SIR model from Alain Vandormael (2018): 
             'We used realistic parameter values for the SIR model, based on earlier HIV studies that have been undertaken in the sub-Saharan Africa context. 
             To this extent, we varied `c` within the range of 50 to 120 sexual acts per year based on data collected from serodiscordant couples across 
             eastern and southern African sites. Previous research has shown considerable heterogeneity in the probability of HIV transmission per sexual contact, 
             largely due to factors associated with the viral load level, genital ulcer disease, stage of HIV progression, condom use, circumcision and use of ART. 
             Following a systematic review of this topic by Boily et al., we selected values for `beta` within the range of 0.003–0.008. 
             Here, we chose values for `v` within the range of 0.15–0.35, which are slightly conservative, but supported by population-based estimates from the 
             sub-Saharan African context."),
           HTML(paste("<div class='code'>", 
                      "<div class='flex'><div class='definition'>c</div><div>varies from 50 to 120 per year</div></div>",
                      "<div class='flex'><div class='definition'>beta</div><div>varies from 0.003 to 0.008</div></div>",
                      "<div class='flex'><div class='definition'>prev</div><div>which here is population prevalence of unsuppressed VL, varies from 0.15 to 0.35</div></div>",
                      "<div class='flex'><div class='definition'>Sv</div><div>could be parameterized using the RV144 Thai Trial results: VE = 61% at 12 months, 31% at 42 months, but below we start with 30% and no waning efficacy. A vaccine duration parameter is not needed because we are only modeling a 3 year trial without boosters.</div></div>",
                      "</div>")),
           HTML("</div>"),
           titlePanel(htmlTemplate("template.html"))
           
  )
}

#------------------------------------------------------------------------------
# for creating Model Fitting tab content (Fitting the model to specific trial data)
#------------------------------------------------------------------------------

getTestTab <- function() {
  tabPanel("Model fitting", 
           
           HTML("<div class='mainPanel'>"),
           sidebarPanel(  
             sliderInput('lambdaTest', 'lambda:', min=0.000005, max=0.0001,
                         value=0.000028, step=0.000001, round=FALSE),
             sliderInput('epsilonTest', 'epilson:', min=0.0, max=1.0,
                         value=0.40, step=0.05, round=FALSE),
             sliderInput('riskTest', 'risk:', min=0, max=30,
                         value=10.0, step=1, round=FALSE),
             sliderInput('numExecution', '# of execution:', min=0, max=1000,
                         value=100, step=50, round=FALSE),
             class = "slider"
           ),
           
           mainPanel(
             p("test Lorem ipsum dolor sit amet, consectetur adipiscing elit. Pellentesque fringilla aliquam ex, vitae scelerisque felis semper quis. Aenean nec pharetra ligula. Mauris vulputate purus ante, id faucibus leo facilisis sit amet. Fusce vestibulum justo eu enim commodo consectetur. Fusce nisi urna, ultrices at purus at, imperdiet efficitur ligula. Curabitur quis sapien ligula. Ut non orci ullamcorper, pulvinar nibh vel, molestie velit. Morbi vulputate hendrerit mi, a mollis risus blandit eget. Cras lacinia eget massa condimentum finibus. Morbi porta lorem augue, in sagittis orci vehicula vel. Sed ipsum nisi, scelerisque quis luctus at, efficitur sed erat. Nam aliquet hendrerit laoreet. Nulla rutrum, nisi pulvinar placerat eleifend, lacus metus ornare dolor, at iaculis augue mi sit amet lorem. Aliquam sit amet turpis nec quam aliquet pretium. Maecenas leo lectus, efficitur in magna ut, gravida iaculis felis."),
             plotOutput("plotTestLambdaRisk")  %>% withSpinner(color="#0dc5c1"),
             p("Lorem ipsum dolor sit amet, consectetur adipiscing elit. Pellentesque fringilla aliquam ex, vitae scelerisque felis semper quis. Aenean nec pharetra ligula. Mauris vulputate purus ante, id faucibus leo facilisis sit amet. Fusce vestibulum justo eu enim commodo consectetur. Fusce nisi urna, ultrices at purus at, imperdiet efficitur ligula. Curabitur quis sapien ligula. Ut non orci ullamcorper, pulvinar nibh vel, molestie velit. Morbi vulputate hendrerit mi, a mollis risus blandit eget. Cras lacinia eget massa condimentum finibus. Morbi porta lorem augue, in sagittis orci vehicula vel. Sed ipsum nisi, scelerisque quis luctus at, efficitur sed erat. Nam aliquet hendrerit laoreet. Nulla rutrum, nisi pulvinar placerat eleifend, lacus metus ornare dolor, at iaculis augue mi sit amet lorem. Aliquam sit amet turpis nec quam aliquet pretium. Maecenas leo lectus, efficitur in magna ut, gravida iaculis felis."),
             plotOutput("plotTestEspilonRisk")  %>% withSpinner(color="#0dc5c1"),
             p("Lorem ipsum dolor sit amet, consectetur adipiscing elit. Pellentesque fringilla aliquam ex, vitae scelerisque felis semper quis. Aenean nec pharetra ligula. Mauris vulputate purus ante, id faucibus leo facilisis sit amet. Fusce vestibulum justo eu enim commodo consectetur. Fusce nisi urna, ultrices at purus at, imperdiet efficitur ligula. Curabitur quis sapien ligula. Ut non orci ullamcorper, pulvinar nibh vel, molestie velit. Morbi vulputate hendrerit mi, a mollis risus blandit eget. Cras lacinia eget massa condimentum finibus. Morbi porta lorem augue, in sagittis orci vehicula vel. Sed ipsum nisi, scelerisque quis luctus at, efficitur sed erat. Nam aliquet hendrerit laoreet. Nulla rutrum, nisi pulvinar placerat eleifend, lacus metus ornare dolor, at iaculis augue mi sit amet lorem. Aliquam sit amet turpis nec quam aliquet pretium. Maecenas leo lectus, efficitur in magna ut, gravida iaculis felis."),
             plotOutput("plotTestEpsilonLambda")  %>% withSpinner(color="#0dc5c1"),
             class = "plotPanel"
             
           ),
           HTML("</div>"),
           titlePanel(htmlTemplate("template.html"))
           
  )
}



