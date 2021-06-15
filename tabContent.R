


#------------------------------------------------------------------------------
# for creating Intro tab content
#------------------------------------------------------------------------------
getIntroContent <- function() {
  return (tabPanel("Introduction",
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
    HTML("</div>")
  ))
}


#------------------------------------------------------------------------------
# for creating Model setup tab content
#------------------------------------------------------------------------------
getModelSetupContent <- function() {
  
  return(tabPanel("Model setup",
           HTML("<div class='content'>"),
           p(paste("We are modeling a vaccine trial using an SI deterministic compartmental model.", 
                   "We are not modeling infections from the I to S compartments but rather only infections from the outside (non-trial) population, ", 
                   "with the infection rate based on the population prevalence `prev` (of viremic individuals), the exposure rate ",         "with the infection rate based on the population prevalence `prev` (of viremic individuals), the exposure rate ",
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
           
           
  ))
  
}


#------------------------------------------------------------------------------
# for creating Calibration tab content
#------------------------------------------------------------------------------
getCalibrationContent <- function() {
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
  )
}


#------------------------------------------------------------------------------
# for creating simple Plot tab content
#------------------------------------------------------------------------------
getSimplePlotContent <- function() {
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
                 column(6, plotOutput("plotCInfect")),
                 column(6, plotOutput("plot2"))),
               fluidRow(
                 column(6, plotOutput("plot3")),
                 column(6, plotOutput("plot4"))),
               plotOutput("plot5")
             )
           )
  )
}



#------------------------------------------------------------------------------
# for creating parameter sweep content
#------------------------------------------------------------------------------
getParameterSweepContent <- function() {
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
}

