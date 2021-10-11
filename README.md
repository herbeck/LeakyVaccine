Leaky vaccine 
=============

Modeling the effects of exposure heterogeneity on vaccine clinical efficacy  

IDM:  
Josh Herbeck (jherbeck@idmod.org)  
Adam Akullian (aakullian@idmod.org)    
Allen Roberts  
David Kong  
Minerva Enriquez  

FHCRC:  
Paul Edlefsen  (pedlefse@fredhutch.org)  

---  

### Background  

It has been hypothesized that exposure heterogeneity can affect estimates of clinical vaccine efficacy for leaky vaccines (e.g. Halloran et al., 1992; White et al., 2010; O'Hagan et al.,2013; Edlefsen, 2014; Coley et al., 2016; Gomes et al., 2016; Kahn et al., 2018; Langwig et al., 2019). Exposure heterogeneity can be broadly characterized as within- or across-population heterogeneity in infection risk.  

Within-population heterogeneity is the variation in risk of HIV infection within a population:  some individuals are at higher risk of infection, due to some combination of higher contact rate (e.g. number of sexual partners), higher per-exposure probability of transmission, or higher HIV prevalence in sexual partners. If this pattern exists within the vaccine and placebo arms of clinical trial, it can result in decreasing clinical vaccine efficacy over the course of the trial. This happens as high risk individuals in both arms are infected (and effectively removed from the susceptible population) at a higher rate than the low risk individuals; incidence declines over the course of this depletion, as the high-risk individuals get infected and only the lower risk individuals remain. If the vaccine at trial has some effect, this incidence decline occurs faster in the placebo arm, resulting in vaccine and placebo arms with unbalanced risk structure.

Across-population heterogeneity describes a situation where two or more populations have different forces of infection (e.g. there is variation in the background incidence or exposure rate). For leaky vaccines, which in theory partially protect all individuals on a per-exposure basis, repeated exposures will lead to declining vaccine efficacy; in populations with high HIV risk, the cumulative effect of multiple exposures can end up lower than the per-exposure vaccine efficacy (how much protection the vaccine provides for a single exposure). This situation may describe HIV vaccine trials in South Africa and Thailand; even if the per-exposure efficacy of the vaccine is the same, would we expect substantial differences in clinical VE that are due to the different incidences in the trial settings?

To quote from Gomes et al., 2016:  "This effect is more pronounced in the control group as individuals within it experience higher rates of infection overall. Consequently, the  ratio of disease rates in vaccinated over control groups increases, and vaccine efficacy, as measured by simple rate ratios, decreases as the trial progresses. Finally, the magnitude of this effect increases with the intensity of transmission." 

### Our goal  

Here we use epidemic models to simulate this process, within and across populations, in the context of HIV prevention trials or longitudinal studies. A key distinction for following our model and analyses is between the per-exposure (or per-contact) efficacy and the clinical efficacy (i.e. the trial outcome). Some initial questions that we address include:  

1. Did exposure heterogeneity contribute to the differences between the RV144 and HVTN 702 HIV vaccine trial outcomes?

2. Can exposure heterogeneity explain waning efficacies seen in other HIV prevention trials (e.g. the AMP VRC01 bnAb trials and the different results seen in the sub-studies, 703 and 704).  

3. In HIV cohort studies incidence often declines over the course of the study. How much of this effect may be due to frailty bias (i.e. individuals with high risk exposure or high exposure rates becoming infected early in the observation period, while individuals with lower risk become infected later)?  

### Model setup  

To simulate an HIV vaccine trial we use a simple deterministic compartmental model. The model includes two compartments:  S, susceptible individuals; and I, infected individuals. Individuals start as S and move to I over the course of the trial if they get infected. We do not model infections back from I to S; we assume that changes in the size of I do not affect the infection rate of S.  

The infection rate of individuals in S is based on: `prev`, the population prevalence (of viremic individuals); `c`, the exposure rate (serodiscordant sexual contacts per time); and `p`, the per-exposure transmission probability.  

The per-exposure (i.e. per-contact) effect of vaccination is `epsilon`, and with this iteration of the model `epsilon` is:  1) not time-varying (the per contact vaccine effect does not decay over time); and 2) assumes a homogeneous effect (does not vary by mark / viral genotype). This model structure also removes the possibility of indirect effects from vaccination.  

The risk structure is controlled by the size of the high-, medium-, and low-risk subgroups, and by `risk`, the risk multiplier, which is used to increase or decrease these risk subgroups.  

`beta` = transmission rate (per contact)   
`c` = exposure rate (serodiscordant sexual contacts per time)  
`prev` = prevalence  (prevalence of viremic individuals)  
`lambda = beta * c * prev`  
`risk` = risk multiplier  
`epsilon` = per-exposure vaccine efficacy; the vaccine-induced reduction in the risk of HIV infection from a single exposure  

The model's basic equations are:  

`dS/dt = -lambda*S`   
`dI/dt = lambda*S`  

The basic compartments are:  

Sp = susceptible placebo  
Ip = infected placebo  
Sv = susceptible vaccinated  
Iv = infected vaccinated  

Svh = susceptible vaccinated high exposure  
Svm = susceptible vaccinated medium exposure  
SvL = susceptible vaccinated low exposure    
Ivh = infected vaccinated high exposure  
Ivm = infected vaccinated medium exposure  
Ivl = infected vaccinated low exposure  

We use the EpiModel framework, http://www.epimodel.org/, from Sam Jenness (Emory University) to build the model.  

### First model demonstration  

Our first pass at the size of the high-, medium-, and low-risk subgroups are: 10% high risk, 80% medium risk, and 10% no (zero) risk. (This parameterization is tough:  Dimitrov et al 2015 even suggest that the MAJORITY of individuals in trials are NOT exposed; https://pubmed.ncbi.nlm.nih.gov/25569838/)  

Sp = susceptible placebo  
Ip = infected placebo  
Sv = susceptible vaccinated  
Iv = infected vaccinated  

Svh = susceptible vaccinated high exposure  
Svm = susceptible vaccinated medium exposure  
SvL = susceptible vaccinated low exposure    
Ivh = infected vaccinated high exposure  
Ivm = infected vaccinated medium exposure  
Ivl = infected vaccinated low exposure  

We use the EpiModel framework, http://www.epimodel.org/, from Sam Jenness (Emory University) to build the model.
