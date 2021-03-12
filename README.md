Leaky vaccine 
=============

An attempt to use a compartmental model to evaluate the effects of exposure risk heterogeneity on vaccine efficacy estimates in HIV vaccine trials.

Josh Herbeck  
Paul Edlefsen  
Molly Reid   
Sam Jenness  

The model is here:  Leaky.vaccine.SI.trial.Rmd

It is hypothesized that exposure heterogeneity (i.e. infection risk heterogeneity) can affect efficacy estimation for leaky vaccines (e.g. Halloran et al., 1992; White et al., 2010; O'Hagan et al.,2013; Edlefsen, 2014; Coley et al., 2016; Gomes et al., 2016; Kahn et al., 2018; Langwig et al., 2019). Our goal is to make a simple deterministic compartmental model to facilitate straightforward simulation-based evaluation of this process within and across populations, in the context of HIV prevention trials or longitudinal studies.  

1. In acute infection studies it seems like many participants get infected early. What is the magnitude of this effect that might be due to frailty bias? 

2. Assess if this effect might contribute to the differences between the RV 144 and HVTN 702 vaccine trial outcomes.  (There has been a couple of analyses of this, and we can build on this and make future analyses of other trial results more straightforward to evaluate.)

3. Assess if this effect might contribute to the waning efficacies seen in HIV prevention trials (specifically the AMP VRC01 bnAb trial).  

4. In the context of the AMP Trial and the different results seen in the sub-studies (703 vs 704); is this due to different forces of infection between the populations?    

5. Continue to raise awareness of this issue to HIV prevention trials, with the ultimate goal of better design and interpretation of efficacy outcomes.   

From Gomes et al., 2016:  "This effect is more pronounced in the control group as individuals within it experience higher rates of infection overall. Consequently, the ratio of disease rates in vaccinated over control groups increases, and vaccine efficacy, as measured by simple rate ratios, decreases as the trial progresses. Finally, the magnitude of this effect increases with the intensity of transmission." 

---

The goal of this little project is to assess the impact of exposure heterogeneity on vaccine efficacy measurements. We know that exposure heterogeneity (EH) can lead to waning efficacy… but we don’t know how much EH has an appreciable (noticeable?) impact in HIV prevention trials. 

We also are talking about two different kinds of EH:  **one** is the EH within a population with 2 trial arms, and right now my model shows (as others have done previously, e.g. Halloran et al., 1992; White et al., 2010; O'Hagan et al.,2013; Edlefsen, 2014; Coley et al., 2016; Gomes et al., 2016; Kahn et al., 2018; Langwig et al., 2019) that EH leads to waning vaccine efficacy (VE) estimates; the **other** is between different populations that have different forces of infection (background incidence / exposure rates). The latter is analogous to comparing vaccine trials between South Africa and Thailand; even if the per-contact efficacy of the vaccine is the same, would we expect substantial differences in population/trial estimates of VE that are due to the different incidence (exposure rates) in the two trial settings? This model can be used to look at both types of EH and VE.

(This SI model is a proof of concept. Maybe to be followed by a Shiny app that anybody can play with.)

The current set of figures in the document https://github.com/herbeck/LeakyVaccine/blob/main/Leaky.vaccine.SI.trial.md show just that that VE will vary between 2 populations if the populations are with and w/out exposure heterogeneity. It also shows that EH causes the VE to wane. This just confirms theory, but it also shows that the model is working at that level.

The next steps are to use the model to ask questions about specific trials, in order to ask questions about HVTN702, RV144, AMP trials, other. The goal is to identify parameter sets (e.g. how much exposure heterogeneity) that would result in substantial differences in measured VE. An example would be, given the AMP trial results that show waning efficacy by year 3… what level of EH could explain this? 

I thought that we would first need to calibrate to the observed incidence in each arm, for each trial. I accept that the models should reflect reality (w/r/t incidence) at some level. And I thought that these initial calibrations would calibrate on `lambda` (`lambda` is a proxy for overall force of infection, as it combines contact rate, per contact transmission rate, and prevalence of viremic individuals). (Notice that in the model the rate of Susceptible -> Infected is not dependent on the number of Infected individuals, because we are modeling just the vaccine trial population.)

I thought that we would first we needed to calibrate to the observed incidence in each arm, for each trial. I accept that the models should reflect reality (w/r/t incidence) at some level. And I thought that these initial calibrations would calibrate on `lambda` (`lambda` is a proxy for overall force of infection, as it combines contact rate, per contact transmission rate, and prevalence of viremic individuals). (Notice that in the model the rate of Susceptible -> Infected is not dependent on the number of Infected individuals, because we are modeling just the vaccine trial population.)

And then we would calibrate on exposure heterogeneity, by fixing the `lambda` in both vaccine and placebo arms but setting the target stats as the incidence for the vaccine and placebo arms (separately). The free parameters would be the `risk` multiplier and the size of the `high risk subgroup`. 

But, now I am thinking that we could skip the calibration of `lambda` altogether and go straight to the calibration on exposure heterogeneity. And also, that the placebo and vaccine arms of the trials will have the same level of exposure heterogeneity (in each population, South Africa or Thailand, etc.). So do we need to calibrate with VE as the target stat? 

Or, we just calibrate the placebo arm, assuming homogeneous exposure, and then run the vaccine trial (add a vaccine and estimate VE using incidence from the placebo and vaccine arms). Then we add in more and more exposure heterogeneity until we see an impact on VE? Also, I think to “calibrate” is maybe the wrong word, because we basically want to explore the parameter space (of EH) that is consistent with observed waning VE (over the course of a trial) or observed lower VE (between two trials with different background incidence). Also, there are multiple mechanisms here:  with waning VE, it is depletion of high risk placebo individuals and violation of proportional hazard assumptions; with lower VE it is depletion of high risk vaccinees due to vaccine leakiness. (Also, cumulative incidence VE will wane even without EH, so...) 
 
Also, the EH parameters that I am using: risk multiplier and/or the size of the high risk population are just two different ways to define exposure heterogeneity. In this simple SI model we are just setting the size of the high risk compartment (i.e. the fraction of people who have accumulate the most exposures) OR set the number of exposures that the high risk people actually have (this would be the risk multiplier).
