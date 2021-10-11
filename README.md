Leaky vaccine 
=============

Modeling the effects of exposure heterogeneity on vaccine clinical efficacy

IDM:  
Josh Herbeck  
Adam Akullian  
Allen Roberts  
David Kong  
Minerva Enriquez  

FHCRC:  
Paul Edlefsen  

---

It is hypothesized that exposure heterogeneity can affect estimates of clinical vaccine efficacy for leaky vaccines (e.g. Halloran et al., 1992; White et al., 2010; O'Hagan et al.,2013; Edlefsen, 2014; Coley et al., 2016; Gomes et al., 2016; Kahn et al., 2018; Langwig et al., 2019). 

Exposure heterogeneity can be broadly characterized as within- or across-population. Within-population heterogeneity is the variation in risk of HIV infection within a population:  some individuals are at higher risk of infection, due to some combination of higher contact rate (e.g. number of sexual partners), higher per-exposure probability of transmission (e.g. due to coinfections or individual traits that effect this probability), or higher HIV prevalence in sexual partners (e.g. a transmission network with higher rates of infection). If this pattern exists within the vaccine and placebo arms of clinical trial, this heterogeneity can lead to decreases in clinical vaccine efficacy over the course of the trial. This happens as high risk individuals in both arms are infected (and effectively removed from the trial susceptible population) at a higher rate than the low risk individuals; if the vaccine (or other prevention method) has some effect, this depletion of high risk individuals occurs faster in the placebo arm, resulting in vaccine and placebo arms with unbalanced risk structure. In short, with exposure heterogeneity the placebo arm incidence declines faster than the vaccine arm incidence, and the estimated clinical vaccine efficacy declines as well (clinical vaccine efficacy = 1 - vaccine incidence/placebo incidence).

; the **other** is between populations with different forces of infection (background incidence / exposure rates). The latter is analogous to comparing vaccine trials between South Africa and Thailand; even if the per-exposure efficacy of the vaccine is the same, would we expect substantial differences in clinical VE that are due to the different incidences in the trial settings? Our models can be used to look at both types of EH and VE.


Here we use epidemic models to simulate this process, within and across populations, in the context of HIV prevention trials or longitudinal studies. Some initial questions that we hoped to address include:  

1. Did exposure heterogeneity contribute to the differences between the RV144 and HVTN 702 HIV vaccine trial outcomes?

2. Can exposure heterogeneity explain waning efficacies seen in other HIV prevention trials (e.g. the AMP VRC01 bnAb trials and the different results seen in the sub-studies, 703 and 704).  

3. In HIV cohort studies incidence often declines over the course of the study. How much of this effect may be due to frailty bias (i.e. individuals with high risk exposure or high exposure rates becoming infected early in the observation period, while individuals with lower risk become infected later)? 

From Gomes et al., 2016:  "This effect is more pronounced in the control group as individuals within it experience higher rates of infection overall. Consequently, the  ratio of disease rates in vaccinated over control groups increases, and vaccine efficacy, as measured by simple rate ratios, decreases as the trial progresses. Finally, the magnitude of this effect increases with the intensity of transmission." 

---

Our goal to assess the impact of exposure heterogeneity on vaccine efficacy measurements. A key distinction is between the per-exposure (or per-contact) efficacy and the clinical efficacy (i.e. the trial outcome). Theory that has been developed over the last 20+ years suggests that exposure heterogeneity (EH) can lead to waning clinical efficacy and diminished population-level effectiveness (e.g. Halloran et al., 1992; White et al., 2010; O'Hagan et al.,2013; Edlefsen, 2014; Coley et al., 2016; Gomes et al., 2016; Kahn et al., 2018; Langwig et al., 2019). But we don’t know how much impact EH can have in HIV prevention trials. 

There are different mechanisms of EH:  **one** is the EH within a population with 2 trial arms, which leads to waning clinical vaccine efficacy (VE); the **other** is between populations with different forces of infection (background incidence / exposure rates). The latter is analogous to comparing vaccine trials between South Africa and Thailand; even if the per-exposure efficacy of the vaccine is the same, would we expect substantial differences in clinical VE that are due to the different incidences in the trial settings? Our models can be used to look at both types of EH and VE.

