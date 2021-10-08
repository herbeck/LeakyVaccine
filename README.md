Leaky vaccine 
=============

Modeling the effects of exposure heterogeneity on vaccine efficacy

IDM:  
Josh Herbeck  
Adam Akullian  
Allen Roberts  
David Kong  
Minerva Enriquez  

FHCRC:  
Paul Edlefsen  

---

It is hypothesized that exposure heterogeneity (i.e. variation in infection risk) can affect efficacy estimation for leaky vaccines (e.g. Halloran et al., 1992; White et al., 2010; O'Hagan et al.,2013; Edlefsen, 2014; Coley et al., 2016; Gomes et al., 2016; Kahn et al., 2018; Langwig et al., 2019). Our goal is to make a simple deterministic compartmental model to facilitate straightforward simulation-based evaluation of this process within and across populations, in the context of HIV prevention trials or longitudinal studies.  

1. Did exposure heterogeneity contribute to the differences between the RV144 and HVTN 702 HIV vaccine trial outcomes?

2. Can exposure heterogeneity explain waning efficacies seen in other HIV prevention trials (e.g. the AMP VRC01 bnAb trials and the different results seen in the sub-studies, 703 and 704).  

3. In HIV cohort studies incidence often declines over the course of the study. How much of this effect may be due to frailty bias (i.e. individuals with high risk exposure or high exposure rates becoming infected early in the observation period, while individuals with lower risk become infected later)? 

From Gomes et al., 2016:  "This effect is more pronounced in the control group as individuals within it experience higher rates of infection overall. Consequently, the  ratio of disease rates in vaccinated over control groups increases, and vaccine efficacy, as measured by simple rate ratios, decreases as the trial progresses. Finally, the magnitude of this effect increases with the intensity of transmission." 

---

Our goal to assess the impact of exposure heterogeneity on vaccine efficacy measurements. A key distinction is between the per-exposure (or per-contact) efficacy and the clinical efficacy (i.e. the trial outcome). Theory that has been developed over the last 20+ years suggests that exposure heterogeneity (EH) can lead to waning clinical efficacy and diminished population-level effectiveness (e.g. Halloran et al., 1992; White et al., 2010; O'Hagan et al.,2013; Edlefsen, 2014; Coley et al., 2016; Gomes et al., 2016; Kahn et al., 2018; Langwig et al., 2019). But we don’t know how much impact EH can have in HIV prevention trials. 

There are different mechanisms of EH:  **one** is the EH within a population with 2 trial arms, which leads to waning clinical vaccine efficacy (VE); the **other** is between populations with different forces of infection (background incidence / exposure rates). The latter is analogous to comparing vaccine trials between South Africa and Thailand; even if the per-exposure efficacy of the vaccine is the same, would we expect substantial differences in clinical VE that are due to the different incidences in the trial settings? Our models can be used to look at both types of EH and VE.

