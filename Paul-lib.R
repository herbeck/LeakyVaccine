library(deSolve)
library(tidyverse)
library(EpiModel)
library(survival)
library(EasyABC)
library(ggplot2)
library(pdfCluster)

si.ode.threegroup.fn <- function ( times, init, param ) {
  with(as.list(c(init, param)), {
    
    # Flows
    # the number of people moving from the S to I compartment at each time step
    
    # The one trial
    
    #PLACEBO arm
    #Susceptible, Infected, placebo, high, medium, low
    #SIph.flow <- risk*lambda*Sph # line from original model, FYI
    SIph.flow <- high.risk.multiplier*(10^log10lambda)*Sph
    SIpm.flow <- (10^log10lambda)*Spm
    SIpl.flow <- 0*(10^log10lambda)*Spl  #0 to give this group zero exposures
       # Could also use "1/high.risk.multiplier" if we don't want ZERO exposures
    
    #VACCINE arm
    #Susceptible, Infected, vaccine, high, medium, low
    SIvh.flow <- high.risk.multiplier*(10^log10lambda)*(1-epsilon)*Svh
    SIvm.flow <- (10^log10lambda)*(1-epsilon)*Svm
    SIvl.flow <- 0*(10^log10lambda)*(1-epsilon)*Svl  #0 to give this group zero exposures
       # Could also use "1/high.risk.multiplier" if we don't want ZERO exposures
    
    # ODEs
    # placebo; heterogeneous high.risk.multiplier
    # original ODE:  dSph <- -SIph.flow
    dSph <- -SIph.flow
    dIph <- SIph.flow  #high.risk.multiplier*lambda*Sph
    dSpm <- -SIpm.flow
    dIpm <- SIpm.flow  #lambda*Spm
    dSpl <- -SIpl.flow
    dIpl <- SIpl.flow  #0*lambda*Spl
    
    # vaccine; heterogeneous high.risk.multiplier
    dSvh <- -SIvh.flow
    dIvh <- SIvh.flow  #high.risk.multiplier*lambda*(1-epsilon)*Svh
    dSvm <- -SIvm.flow
    dIvm <- SIvm.flow  #lambda*Svm
    dSvl <- -SIvl.flow
    dIvl <- SIvl.flow  #0*lambda*(1-epsilon)*Svl


    #Output
    list(c(
           dSph,dIph,
           dSpm,dIpm,
           dSpl,dIpl,
           dSvh,dIvh,
           dSvm,dIvm,
           dSvl,dIvl,
           SIph.flow,SIpm.flow,SIpl.flow,
           SIvh.flow,SIvm.flow,SIvl.flow
           ))
  })
} # si.ode.threegroup.fn (..)
mod.manipulate.threegroup <- function( mod ) {
  #browser()

  # THREEGROUP
  mod <- mutate_epi(mod, total.Svh.Svm.Svl = Svh + Svm + Svl) #all susceptible in heterogeneous risk vaccine pop
  mod <- mutate_epi(mod, total.Sph.Spm.Spl = Sph + Spm + Spl) #all susceptible in heterogeneous risk placebo pop
  mod <- mutate_epi(mod, total.Ivh.Ivm.Ivl = Ivh + Ivm + Ivl) #all infected in heterogeneous risk vaccine pop
  mod <- mutate_epi(mod, total.Iph.Ipm.Ipl = Iph + Ipm + Ipl) #all infected in heterogeneous risk placebo pop
  mod <- mutate_epi(mod, total.SIvh.SIvm.SIvl.flow = SIvh.flow + SIvm.flow + SIvl.flow) #all infections per day in heterogeneous risk vaccine pop
  mod <- mutate_epi(mod, total.SIph.SIpm.SIpl.flow = SIph.flow + SIpm.flow + SIpl.flow) #all infections in heterogeneous risk placebo pop
  
  #Incidence estimates, per 100 person years
  #Instantaneous incidence / hazard
  mod <- mutate_epi(mod, rate.Vaccine.het = (total.SIvh.SIvm.SIvl.flow/total.Svh.Svm.Svl)*365*100)
  mod <- mutate_epi(mod, rate.Placebo.het = (total.SIph.SIpm.SIpl.flow/total.Sph.Spm.Spl)*365*100)
  
  #Cumulative incidence
  mod <- mutate_epi(mod, cumul.Svh.Svm.Svl = cumsum(total.Svh.Svm.Svl))
  mod <- mutate_epi(mod, cumul.Sph.Spm.Spl = cumsum(total.Sph.Spm.Spl))
  mod <- mutate_epi(mod, cumul.rate.Vaccine.het = (total.Ivh.Ivm.Ivl/cumul.Svh.Svm.Svl)*365*100)
  mod <- mutate_epi(mod, cumul.rate.Placebo.het = (total.Iph.Ipm.Ipl/cumul.Sph.Spm.Spl)*365*100)
  
  #Vaccine efficacy (VE) estimates
  #VE <- 1 - Relative Risk; this is VE for instantaneous incidence / hazard
  mod <- mutate_epi(mod, VE.inst = 100 * ( 1 - rate.Vaccine.het/rate.Placebo.het ) )
  
  #VE <- 1 - Relative Risk; this is VE from cumulative incidence
  mod <- mutate_epi(mod, VE.cumul = 100 * ( 1 - cumul.rate.Vaccine.het/cumul.rate.Placebo.het ) )

  return( mod );
} # mod.manipulate.threegroup (..)

run.and.compute.run.stats.threegroup <- function (
      epsilon,   #per contact vaccine efficacy
      log10lambda,     #log10( beta*c*prev ),
      high.risk.multiplier,          # Risk multiplier for high risk group
      highRiskProportion,
      lowRiskProportion,             # This is a proportion among those not high risk
      vaccinatedProportion = 0.5,  # In lieu of naming vaccine and placebo arms separately (and their N)
      trialSize = 10000,  # Now we just have to add this magic number for size
      trial.evaluation.time = 3*365
) {
      # Paul added the other params to this (risk here, others below):
      param <- param.dcm(epsilon = epsilon, log10lambda = log10lambda, high.risk.multiplier = high.risk.multiplier );

      # initial values
      Svh <- floor( highRiskProportion * vaccinatedProportion * trialSize );
      Sph <- floor( highRiskProportion * ( 1.0 - vaccinatedProportion ) * trialSize );
      Svl <- floor( ( 1.0 - highRiskProportion ) * lowRiskProportion * vaccinatedProportion * trialSize );
      Spl <- floor( ( 1.0 - highRiskProportion ) * lowRiskProportion * ( 1.0 - vaccinatedProportion ) * trialSize );
      Svm <- floor( ( 1.0 - highRiskProportion ) * ( 1.0 - lowRiskProportion ) * vaccinatedProportion * trialSize );
      Spm <- floor( ( 1.0 - highRiskProportion ) * ( 1.0 - lowRiskProportion ) * ( 1.0 - vaccinatedProportion ) * trialSize );

      Sp <- Spl + Spm + Sph;
      Sv <- Svl + Svm + Svh;
      if( vaccinatedProportion == 0.5 ) {
          stopifnot( Sp == Sv );
      }

      init <- init.dcm(Sph = Sph, Iph = 0,    #placebo, high risk
                       Spm = Spm, Ipm = 0,   #placebo, medium risk
                       Spl = Spl, Ipl = 0,    #placebo, low risk
                       Svh = Svh, Ivh = 0,    #vaccine
                       Svm = Svm, Ivm = 0,   #vaccine
                       Svl = Svl, Ivl = 0,    #vaccine
                       SIph.flow = 0, SIpm.flow = 0, SIpl.flow = 0,
                       SIvh.flow = 0, SIvm.flow = 0, SIvl.flow = 0
                       );
      
      control <- control.dcm( nsteps = trial.evaluation.time, new.mod = si.ode.threegroup.fn );
      mod <- dcm( param, init, control );
      #print( mod )
      
      mod.with.stats <- mod.manipulate.threegroup( mod );
      #print( mod.with.stats )
      mod.with.stats.df <- as.data.frame( mod.with.stats );
      
      # heterogeneous risk using cumulative VE:
      VE <- mod.with.stats.df$VE.cumul[ trial.evaluation.time ];
      
      ## The placebo incidence out stat vector is the _cumulative_ incidence at time trial.evaluation.time.
      placeboIncidence <- mod.with.stats.df$cumul.rate.Placebo.het[ trial.evaluation.time ];

      c( VE = VE, placeboIncidence = placeboIncidence );
} # run.and.compute.run.stats.threegroup (..)

common.parameters <- c( "epsilon" );
trial.parameters <- c( "log10lambda", "high.risk.multiplier", "highRiskProportion", "lowRiskProportion" );
