library(deSolve)
library(tidyverse)
library(EpiModel)
library(survival)
library(EasyABC)
library(ggplot2)
library(pdfCluster)

si.ode.onetrial.fn <- function ( times, init, param ) {
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
} # si.ode.onetrial.fn (..)
mod.manipulate.onetrial <- function( mod ) {
  #browser()

  # ONETRIAL
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
} # mod.manipulate.onetrial (..)

run.and.compute.run.stats.onetrial <- function (
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
      
      control <- control.dcm( nsteps = trial.evaluation.time, new.mod = si.ode.onetrial.fn );
      mod <- dcm( param, init, control );
      #print( mod )
      
      mod.with.stats <- mod.manipulate.onetrial( mod );
      #print( mod.with.stats )
      mod.with.stats.df <- as.data.frame( mod.with.stats );
      
      # heterogeneous risk using cumulative VE:
      VE <- mod.with.stats.df$VE.cumul[ trial.evaluation.time ];
      
      ## The placebo incidence out stat vector is the _cumulative_ incidence at time trial.evaluation.time.
      placeboIncidence <- mod.with.stats.df$cumul.rate.Placebo.het[ trial.evaluation.time ];

      c( VE = VE, placeboIncidence = placeboIncidence );
} # run.and.compute.run.stats.onetrial (..)

common.parameters <- c( "epsilon" );
trial.parameters <- c( "log10lambda", "high.risk.multiplier", "highRiskProportion", "lowRiskProportion" );
rv144.parameters <- paste( "rv144", trial.parameters, sep = "." );
hvtn702.parameters <- paste( "hvtn702", trial.parameters, sep = "." );
all.parameters <- c( common.parameters, rv144.parameters, hvtn702.parameters );

#------------------------------------------------------------------------------
# sim execution
#------------------------------------------------------------------------------
run.and.compute.run.stats.rv144.hvtn702 <- function (
      epsilon,   #per contact vaccine efficacy
      rv144.log10lambda,     #log10( beta*c*prev ),
      rv144.high.risk.multiplier,          # Risk multiplier for high risk group
      rv144.highRiskProportion,
      rv144.lowRiskProportion,             # This is a proportion among those not high risk
      hvtn702.log10lambda,
      hvtn702.high.risk.multiplier,
      hvtn702.highRiskProportion,
      hvtn702.lowRiskProportion,
      vaccinatedProportion = 0.5,  # In lieu of naming vaccine and placebo arms separately (and their N)
      trialSize = 10000,  # Now we just have to add this magic number for size
      trial.evaluation.time = 3*365
) {
    rv144.results <- run.and.compute.run.stats.onetrial(
      epsilon = epsilon,   # common
      log10lambda = rv144.log10lambda,
      high.risk.multiplier = rv144.high.risk.multiplier,
      highRiskProportion = rv144.highRiskProportion,
      lowRiskProportion = rv144.lowRiskProportion,
      vaccinatedProportion = vaccinatedProportion,
      trialSize = trialSize,
      trial.evaluation.time = trial.evaluation.time
    );
    hvtn702.results <- run.and.compute.run.stats.onetrial(
      epsilon = epsilon,   # common
      log10lambda = hvtn702.log10lambda,
      high.risk.multiplier = hvtn702.high.risk.multiplier,
      highRiskProportion = hvtn702.highRiskProportion,
      lowRiskProportion = hvtn702.lowRiskProportion,
      vaccinatedProportion = vaccinatedProportion,
      trialSize = trialSize,
      trial.evaluation.time = trial.evaluation.time
    );
    c( rv144.VE = rv144.results[[ "VE" ]], rv144.placeboIncidence = rv144.results[[ "placeboIncidence" ]], hvtn702.VE = hvtn702.results[[ "VE" ]], hvtn702.placeboIncidence = hvtn702.results[[ "placeboIncidence" ]] );
} # run.and.compute.run.stats.rv144.hvtn702 (..)

## OLD, Before refactoring:
si.ode.rv144.hvtn702.fn.old <- function ( times, init, param ) {
  with(as.list(c(init, param)), {
    
    # Flows
    # the number of people moving from the S to I compartment at each time step
    
    # RV144
    
    #PLACEBO arm
    #Susceptible, Infected, placebo, high, medium, low
    #SIph.flow <- risk*lambda*Sph # line from original model, FYI
    rv144.SIph.flow <- rv144.high.risk.multiplier*(10^rv144.log10lambda)*rv144.Sph
    #rv144.SIpm.flow <- rv144.lambda*rv144.Spl ### BUG HERE: Spl instead of Spm!!!
    rv144.SIpm.flow <- (10^rv144.log10lambda)*rv144.Spm ### BUG FIXED
    rv144.SIpl.flow <- 0*(10^rv144.log10lambda)*rv144.Spl  #0 to give this group zero exposures
       # Could also use "1/high.risk.multiplier" if we don't want ZERO exposures
    
    #VACCINE arm
    #Susceptible, Infected, vaccine, high, medium, low
    rv144.SIvh.flow <- rv144.high.risk.multiplier*(10^rv144.log10lambda)*(1-epsilon)*rv144.Svh
    ### Paul found this line has a bug:
    ### rv144.SIvm.flow <- rv144.lambda*(1-epsilon)*rv144.Spl
### DARN it actually had TWO BUGS! Also the l at the end! 
    #rv144.SIvm.flow <- rv144.lambda*(1-epsilon)*rv144.Svl ### BUG HERE: Svl instead of Svm!!!
    rv144.SIvm.flow <- (10^rv144.log10lambda)*(1-epsilon)*rv144.Svm ### BUG FIXED
    rv144.SIvl.flow <- 0*(10^rv144.log10lambda)*(1-epsilon)*rv144.Svl  #0 to give this group zero exposures
       # Could also use "1/high.risk.multiplier" if we don't want ZERO exposures
    
    # ODEs
    # placebo; heterogeneous rv144.high.risk.multiplier
    # original ODE:  dSph <- -SIph.flow
    drv144.Sph <- -rv144.SIph.flow
    drv144.Iph <- rv144.SIph.flow  #rv144.high.risk.multiplier*rv144.lambda*rv144.Sph
    drv144.Spm <- -rv144.SIpm.flow
    drv144.Ipm <- rv144.SIpm.flow  #rv144.lambda*rv144.Spm
    drv144.Spl <- -rv144.SIpl.flow
    drv144.Ipl <- rv144.SIpl.flow  #0*rv144.lambda*rv144.Spl
    
    # vaccine; heterogeneous rv144.high.risk.multiplier
    drv144.Svh <- -rv144.SIvh.flow
    drv144.Ivh <- rv144.SIvh.flow  #rv144.high.risk.multiplier*rv144.lambda*(1-epsilon)*rv144.Svh
    drv144.Svm <- -rv144.SIvm.flow
    drv144.Ivm <- rv144.SIvm.flow  #rv144.lambda*rv144.Svm
    drv144.Svl <- -rv144.SIvl.flow
    drv144.Ivl <- rv144.SIvl.flow  #0*rv144.lambda*(1-epsilon)*rv144.Svl

    # HVTN702
    # the number of people moving from S to I at each time step
    #Susceptible, Infected, placebo, high, medium, low
    hvtn702.SIph.flow <- hvtn702.high.risk.multiplier*(10^hvtn702.log10lambda)*hvtn702.Sph
    hvtn702.SIpm.flow <- (10^hvtn702.log10lambda)*hvtn702.Spm
    hvtn702.SIpl.flow <- 0*(10^hvtn702.log10lambda)*hvtn702.Spl  #0 to give this group zero exposures
    
    #Susceptible, Infected, vaccine, high, medium, low
    hvtn702.SIvh.flow <- hvtn702.high.risk.multiplier*(10^hvtn702.log10lambda)*(1-epsilon)*hvtn702.Svh
    hvtn702.SIvm.flow <- (10^hvtn702.log10lambda)*(1-epsilon)*hvtn702.Svm
    hvtn702.SIvl.flow <- 0*(10^hvtn702.log10lambda)*(1-epsilon)*hvtn702.Svl  #0 to give this group zero exposures
    
    # ODEs
    # placebo; heterogeneous hvtn702.high.risk.multiplier
    dhvtn702.Sph <- -hvtn702.SIph.flow
    dhvtn702.Iph <- hvtn702.SIph.flow  #hvtn702.high.risk.multiplier*hvtn702.lambda*hvtn702.Sph
    dhvtn702.Spm <- -hvtn702.SIpm.flow
    dhvtn702.Ipm <- hvtn702.SIpm.flow  #hvtn702.lambda*hvtn702.Spm
    dhvtn702.Spl <- -hvtn702.SIpl.flow
    dhvtn702.Ipl <- hvtn702.SIpl.flow  #0*hvtn702.lambda*hvtn702.Spl
    
    # vaccine; heterogeneous hvtn702.high.risk.multiplier
    dhvtn702.Svh <- -hvtn702.SIvh.flow
    dhvtn702.Ivh <- hvtn702.SIvh.flow  #hvtn702.high.risk.multiplier*hvtn702.lambda*(1-epsilon)*hvtn702.Svh
    dhvtn702.Svm <- -hvtn702.SIvm.flow
    dhvtn702.Ivm <- hvtn702.SIvm.flow  #hvtn702.lambda*hvtn702.Svm
    dhvtn702.Svl <- -hvtn702.SIvl.flow
    dhvtn702.Ivl <- hvtn702.SIvl.flow  #0*hvtn702.lambda*(1-epsilon)*hvtn702.Svl

    #Output
    list(c(
           drv144.Sph,drv144.Iph,
           drv144.Spm,drv144.Ipm,
           drv144.Spl,drv144.Ipl,
           drv144.Svh,drv144.Ivh,
           drv144.Svm,drv144.Ivm,
           drv144.Svl,drv144.Ivl,
           rv144.SIph.flow,rv144.SIpm.flow,rv144.SIpl.flow,
           rv144.SIvh.flow,rv144.SIvm.flow,rv144.SIvl.flow,

           dhvtn702.Sph,dhvtn702.Iph,
           dhvtn702.Spm,dhvtn702.Ipm,
           dhvtn702.Spl,dhvtn702.Ipl,
           dhvtn702.Svh,dhvtn702.Ivh,
           dhvtn702.Svm,dhvtn702.Ivm,
           dhvtn702.Svl,dhvtn702.Ivl,
           hvtn702.SIph.flow,hvtn702.SIpm.flow,hvtn702.SIpl.flow,
           hvtn702.SIvh.flow,hvtn702.SIvm.flow,hvtn702.SIvl.flow

           ))
  })
} # si.ode.rv144.hvtn702.fn.old (..)

mod.manipulate.rv144.hvtn702.old <- function( mod ) {
  #browser()

  # RV144
  mod <- mutate_epi(mod, total.rv144.Svh.rv144.Svm.rv144.Svl = rv144.Svh + rv144.Svm + rv144.Svl) #all susceptible in heterogeneous risk vaccine pop
  mod <- mutate_epi(mod, total.rv144.Sph.rv144.Spm.rv144.Spl = rv144.Sph + rv144.Spm + rv144.Spl) #all susceptible in heterogeneous risk placebo pop
  mod <- mutate_epi(mod, total.rv144.Ivh.rv144.Ivm.rv144.Ivl = rv144.Ivh + rv144.Ivm + rv144.Ivl) #all infected in heterogeneous risk vaccine pop
  mod <- mutate_epi(mod, total.rv144.Iph.rv144.Ipm.rv144.Ipl = rv144.Iph + rv144.Ipm + rv144.Ipl) #all infected in heterogeneous risk placebo pop
  mod <- mutate_epi(mod, total.rv144.SIvh.rv144.SIvm.rv144.SIvl.flow = rv144.SIvh.flow + rv144.SIvm.flow + rv144.SIvl.flow) #all infections per day in heterogeneous risk vaccine pop
  mod <- mutate_epi(mod, total.rv144.SIph.rv144.SIpm.rv144.SIpl.flow = rv144.SIph.flow + rv144.SIpm.flow + rv144.SIpl.flow) #all infections in heterogeneous risk placebo pop
  
  #Incidence estimates, per 100 person years
  #Instantaneous incidence / hazard
  mod <- mutate_epi(mod, rv144.rate.Vaccine.het = (total.rv144.SIvh.rv144.SIvm.rv144.SIvl.flow/total.rv144.Svh.rv144.Svm.rv144.Svl)*365*100)
  mod <- mutate_epi(mod, rv144.rate.Placebo.het = (total.rv144.SIph.rv144.SIpm.rv144.SIpl.flow/total.rv144.Sph.rv144.Spm.rv144.Spl)*365*100)
  
  #Cumulative incidence
  mod <- mutate_epi(mod, cumul.rv144.Svh.rv144.Svm.rv144.Svl = cumsum(total.rv144.Svh.rv144.Svm.rv144.Svl))
  mod <- mutate_epi(mod, cumul.rv144.Sph.rv144.Spm.rv144.Spl = cumsum(total.rv144.Sph.rv144.Spm.rv144.Spl))
  mod <- mutate_epi(mod, cumul.rv144.rate.Vaccine.het = (total.rv144.Ivh.rv144.Ivm.rv144.Ivl/cumul.rv144.Svh.rv144.Svm.rv144.Svl)*365*100)
  mod <- mutate_epi(mod, cumul.rv144.rate.Placebo.het = (total.rv144.Iph.rv144.Ipm.rv144.Ipl/cumul.rv144.Sph.rv144.Spm.rv144.Spl)*365*100)
  
  #Vaccine efficacy (VE) estimates
  #VE <- 1 - Relative Risk; this is VE for instantaneous incidence / hazard
  mod <- mutate_epi(mod, rv144.VE.inst = 100 * ( 1 - rv144.rate.Vaccine.het/rv144.rate.Placebo.het ) )
  
  #VE <- 1 - Relative Risk; this is VE from cumulative incidence
  mod <- mutate_epi(mod, rv144.VE.cumul = 100 * ( 1 - cumul.rv144.rate.Vaccine.het/cumul.rv144.rate.Placebo.het ) )


  
  # HVTN702
  mod <- mutate_epi(mod, total.hvtn702.Svh.hvtn702.Svm.hvtn702.Svl = hvtn702.Svh + hvtn702.Svm + hvtn702.Svl) #all susceptible in heterogeneous risk vaccine pop
  mod <- mutate_epi(mod, total.hvtn702.Sph.hvtn702.Spm.hvtn702.Spl = hvtn702.Sph + hvtn702.Spm + hvtn702.Spl) #all susceptible in heterogeneous risk placebo pop
  mod <- mutate_epi(mod, total.hvtn702.Ivh.hvtn702.Ivm.hvtn702.Ivl = hvtn702.Ivh + hvtn702.Ivm + hvtn702.Ivl) #all infected in heterogeneous risk vaccine pop
  mod <- mutate_epi(mod, total.hvtn702.Iph.hvtn702.Ipm.hvtn702.Ipl = hvtn702.Iph + hvtn702.Ipm + hvtn702.Ipl) #all infected in heterogeneous risk placebo pop
  mod <- mutate_epi(mod, total.hvtn702.SIvh.hvtn702.SIvm.hvtn702.SIvl.flow = hvtn702.SIvh.flow + hvtn702.SIvm.flow + hvtn702.SIvl.flow) #all infections per day in heterogeneous risk vaccine pop
  mod <- mutate_epi(mod, total.hvtn702.SIph.hvtn702.SIpm.hvtn702.SIpl.flow = hvtn702.SIph.flow + hvtn702.SIpm.flow + hvtn702.SIpl.flow) #all infections in heterogeneous risk placebo pop
  
  #Incidence estimates, per 100 person years
  #Instantaneous incidence / hazard
  mod <- mutate_epi(mod, hvtn702.rate.Vaccine.het = (total.hvtn702.SIvh.hvtn702.SIvm.hvtn702.SIvl.flow/total.hvtn702.Svh.hvtn702.Svm.hvtn702.Svl)*365*100)
  mod <- mutate_epi(mod, hvtn702.rate.Placebo.het = (total.hvtn702.SIph.hvtn702.SIpm.hvtn702.SIpl.flow/total.hvtn702.Sph.hvtn702.Spm.hvtn702.Spl)*365*100)
  
  #Cumulative incidence
  mod <- mutate_epi(mod, cumul.hvtn702.Svh.hvtn702.Svm.hvtn702.Svl = cumsum(total.hvtn702.Svh.hvtn702.Svm.hvtn702.Svl))
  mod <- mutate_epi(mod, cumul.hvtn702.Sph.hvtn702.Spm.hvtn702.Spl = cumsum(total.hvtn702.Sph.hvtn702.Spm.hvtn702.Spl))
  mod <- mutate_epi(mod, cumul.hvtn702.rate.Vaccine.het = (total.hvtn702.Ivh.hvtn702.Ivm.hvtn702.Ivl/cumul.hvtn702.Svh.hvtn702.Svm.hvtn702.Svl)*365*100)
  mod <- mutate_epi(mod, cumul.hvtn702.rate.Placebo.het = (total.hvtn702.Iph.hvtn702.Ipm.hvtn702.Ipl/cumul.hvtn702.Sph.hvtn702.Spm.hvtn702.Spl)*365*100)
  
  #Vaccine efficacy (VE) estimates
  #VE <- 1 - Relative Risk; this is VE for instantaneous incidence / hazard
  mod <- mutate_epi(mod, hvtn702.VE.inst = 100 * ( 1 - hvtn702.rate.Vaccine.het/hvtn702.rate.Placebo.het ) )
  
  #VE <- 1 - Relative Risk; this is VE from cumulative incidence
  mod <- mutate_epi(mod, hvtn702.VE.cumul = 100 * ( 1 - cumul.hvtn702.rate.Vaccine.het/cumul.hvtn702.rate.Placebo.het ) )

  #return(mod)
} # mod.manipulate.rv144.hvtn702.old (..)

#------------------------------------------------------------------------------
# sim execution
#------------------------------------------------------------------------------
run.and.compute.run.stats.rv144.hvtn702.old <- function (
      epsilon,   #per contact vaccine efficacy
      rv144.log10lambda,     #log10( beta*c*prev ),
      rv144.high.risk.multiplier,          # Risk multiplier for high risk group
      rv144.highRiskProportion,
      rv144.lowRiskProportion,             # This is a proportion among those not high risk
      hvtn702.log10lambda,
      hvtn702.high.risk.multiplier,
      hvtn702.highRiskProportion,
      hvtn702.lowRiskProportion,
      vaccinatedProportion = 0.5,  # In lieu of naming vaccine and placebo arms separately (and their N)
      trialSize = 10000,  # Now we just have to add this magic number for size
      trial.evaluation.time = 3*365
) {
      # Paul added the other params to this (risk here, others below):
      param <- param.dcm(epsilon = epsilon, rv144.log10lambda = rv144.log10lambda, hvtn702.log10lambda = hvtn702.log10lambda, rv144.high.risk.multiplier = rv144.high.risk.multiplier, hvtn702.high.risk.multiplier = hvtn702.high.risk.multiplier )

      # rv144 initial values
      rv144.Svh <- floor( rv144.highRiskProportion * vaccinatedProportion * trialSize );
      rv144.Sph <- floor( rv144.highRiskProportion * ( 1.0 - vaccinatedProportion ) * trialSize );
      rv144.Svl <- floor( ( 1.0 - rv144.highRiskProportion ) * rv144.lowRiskProportion * vaccinatedProportion * trialSize );
      rv144.Spl <- floor( ( 1.0 - rv144.highRiskProportion ) * rv144.lowRiskProportion * ( 1.0 - vaccinatedProportion ) * trialSize );
      rv144.Svm <- floor( ( 1.0 - rv144.highRiskProportion ) * ( 1.0 - rv144.lowRiskProportion ) * vaccinatedProportion * trialSize );
      rv144.Spm <- floor( ( 1.0 - rv144.highRiskProportion ) * ( 1.0 - rv144.lowRiskProportion ) * ( 1.0 - vaccinatedProportion ) * trialSize );

      rv144.Sp <- rv144.Spl + rv144.Spm + rv144.Sph;
      rv144.Sv <- rv144.Svl + rv144.Svm + rv144.Svh;
      if( vaccinatedProportion == 0.5 ) {
          stopifnot( rv144.Sp == rv144.Sv );
      }


      # hvtn702 initial values
      hvtn702.Svh <- floor( hvtn702.highRiskProportion * vaccinatedProportion * trialSize );
      hvtn702.Sph <- floor( hvtn702.highRiskProportion * ( 1.0 - vaccinatedProportion ) * trialSize );
      hvtn702.Svl <- floor( ( 1.0 - hvtn702.highRiskProportion ) * hvtn702.lowRiskProportion * vaccinatedProportion * trialSize );
      hvtn702.Spl <- floor( ( 1.0 - hvtn702.highRiskProportion ) * hvtn702.lowRiskProportion * ( 1.0 - vaccinatedProportion ) * trialSize );
      hvtn702.Svm <- floor( ( 1.0 - hvtn702.highRiskProportion ) * ( 1.0 - hvtn702.lowRiskProportion ) * vaccinatedProportion * trialSize );
      hvtn702.Spm <- floor( ( 1.0 - hvtn702.highRiskProportion ) * ( 1.0 - hvtn702.lowRiskProportion ) * ( 1.0 - vaccinatedProportion ) * trialSize );

      hvtn702.Sp <- hvtn702.Spl + hvtn702.Spm + hvtn702.Sph;
      hvtn702.Sv <- hvtn702.Svl + hvtn702.Svm + hvtn702.Svh;
      if( vaccinatedProportion == 0.5 ) {
          stopifnot( hvtn702.Sp == hvtn702.Sv );
      }

      init.rv144.hvtn702 <- init.dcm(rv144.Sph = rv144.Sph, rv144.Iph = 0,    #placebo, high risk
                       rv144.Spm = rv144.Spm, rv144.Ipm = 0,   #placebo, medium risk
                       rv144.Spl = rv144.Spl, rv144.Ipl = 0,    #placebo, low risk
                       rv144.Svh = rv144.Svh, rv144.Ivh = 0,    #vaccine
                       rv144.Svm = rv144.Svm, rv144.Ivm = 0,   #vaccine
                       rv144.Svl = rv144.Svl, rv144.Ivl = 0,    #vaccine
                       rv144.SIph.flow = 0, rv144.SIpm.flow = 0, rv144.SIpl.flow = 0,
                       rv144.SIvh.flow = 0, rv144.SIvm.flow = 0, rv144.SIvl.flow = 0,

                       hvtn702.Sph = hvtn702.Sph, hvtn702.Iph = 0,    #placebo, high risk
                       hvtn702.Spm = hvtn702.Spm, hvtn702.Ipm = 0,   #placebo, medium risk
                       hvtn702.Spl = hvtn702.Spl, hvtn702.Ipl = 0,    #placebo, low risk
                       hvtn702.Svh = hvtn702.Svh, hvtn702.Ivh = 0,    #vaccine
                       hvtn702.Svm = hvtn702.Svm, hvtn702.Ivm = 0,   #vaccine
                       hvtn702.Svl = hvtn702.Svl, hvtn702.Ivl = 0,    #vaccine
                       hvtn702.SIph.flow = 0, hvtn702.SIpm.flow = 0, hvtn702.SIpl.flow = 0,
                       hvtn702.SIvh.flow = 0, hvtn702.SIvm.flow = 0, hvtn702.SIvl.flow = 0
                       )
      
      control <- control.dcm( nsteps = trial.evaluation.time, new.mod = si.ode.rv144.hvtn702.fn );
      mod <- dcm( param, init.rv144.hvtn702, control );
      #print( mod )
      
      mod.with.stats <- mod.manipulate.rv144.hvtn702( mod );
      #print( mod.with.stats )
      mod.with.stats.df <- as.data.frame( mod.with.stats );
      
      # heterogeneous risk using cumulative VE:
      rv144.VE <- mod.with.stats.df$rv144.VE.cumul[ trial.evaluation.time ];
      hvtn702.VE <- mod.with.stats.df$hvtn702.VE.cumul[ trial.evaluation.time ];
      
      ## The placebo incidence out stat vector is the _cumulative_ incidence at time trial.evaluation.time.
      rv144.placeboIncidence <- mod.with.stats.df$cumul.rv144.rate.Placebo.het[ trial.evaluation.time ];
      hvtn702.placeboIncidence <- mod.with.stats.df$cumul.hvtn702.rate.Placebo.het[ trial.evaluation.time ];

      c( rv144.VE = rv144.VE, rv144.placeboIncidence = rv144.placeboIncidence, hvtn702.VE = hvtn702.VE, hvtn702.placeboIncidence = hvtn702.placeboIncidence );
} # run.and.compute.run.stats.rv144.hvtn702.old (..)


## This computes the distance as calculated internally in the abc function - but not 
# where the distances are normalized using "normalise", which we think centralizes also 
# (so makes each scaled stat have mean 0, sd 1) - the standard deviations are saved and 
# included in the abc output and need to be passed into here, because after keeping only the 
# closest X% of the samples, the remaining samples will have a different STDEV. So to get the 
# right distances it's important to use the correct stdev. These values are not centralized 
# before the distances are computed but this should not matter.

# example: calculate.abc.dist( fit.rej$stats, as.numeric( target.stats ), fit.rej$stats_normalization )
calculate.abc.dist <- function ( sampled.stats.matrix, target.stats, target.stat.stdevs ) {
    nss <- length( target.stats );
    stopifnot( ncol( sampled.stats.matrix ) == nss );

    scaled.sumstat <- sampled.stats.matrix;
    for (j in 1:nss) {
        scaled.sumstat[, j] <- ( sampled.stats.matrix[, j] / target.stat.stdevs[ j ] );
    }
    scaled.target <- target.stats;
    for (j in 1:nss) {
        scaled.target[j] <- ( target.stats[j] / target.stat.stdevs[ j ] );
    }
    sum1 <- 0
    for (j in 1:nss) {
        sum1 <- sum1 + (scaled.sumstat[, j] - scaled.target[j])^2
    }
    dist <- sqrt(sum1)
 
   return( dist );
} # calculate.abc.dist ( .. )


    # Ok, so there's two sets of model-specific parameters (rv144,
    # hvtn702), and then one parameter that is shared (epsilon, the
    # per-contact vaccine efficacy parameter), so the best way to
    # optimize is using conditional optimization. That is, we hold
    # epsilon fixed, and then do each model optimization separately,
    # then hold those params fixed and optimize epsilon, and iterate
    # until convergence.

    make.epsilon.abc.fn <- function ( other.parameters ) {
        function( x ) {
            run.and.compute.run.stats.rv144.hvtn702( epsilon = x, rv144.log10lambda = other.parameters[[ "rv144.log10lambda" ]], rv144.high.risk.multiplier = other.parameters[[ "rv144.high.risk.multiplier" ]], rv144.highRiskProportion = other.parameters[[ "rv144.highRiskProportion" ]], rv144.lowRiskProportion = other.parameters[[ "rv144.lowRiskProportion" ]], hvtn702.log10lambda = other.parameters[[ "hvtn702.log10lambda" ]], hvtn702.high.risk.multiplier = other.parameters[[ "hvtn702.high.risk.multiplier" ]], hvtn702.highRiskProportion = other.parameters[[ "hvtn702.highRiskProportion" ]], hvtn702.lowRiskProportion = other.parameters[[ "hvtn702.lowRiskProportion" ]] )
        }
    }
    make.rv144.abc.fn <- function ( other.parameters ) {
        function( x ) {
            run.and.compute.run.stats.rv144.hvtn702( epsilon = other.parameters[[ "epsilon" ]], rv144.log10lambda = x[ 1 ], rv144.high.risk.multiplier = x[ 2 ], rv144.highRiskProportion = x[ 3 ], rv144.lowRiskProportion = x[ 4 ], hvtn702.log10lambda = other.parameters[[ "hvtn702.log10lambda" ]], hvtn702.high.risk.multiplier = other.parameters[[ "hvtn702.high.risk.multiplier" ]], hvtn702.highRiskProportion = other.parameters[[ "hvtn702.highRiskProportion" ]], hvtn702.lowRiskProportion = other.parameters[[ "hvtn702.lowRiskProportion" ]] )
        }
    }
    make.hvtn702.abc.fn <- function ( other.parameters ) {
        function( x ) {
            run.and.compute.run.stats.rv144.hvtn702( epsilon = other.parameters[[ "epsilon" ]], rv144.log10lambda = other.parameters[[ "rv144.log10lambda" ]], rv144.high.risk.multiplier = other.parameters[[ "rv144.high.risk.multiplier" ]], rv144.highRiskProportion = other.parameters[[ "rv144.highRiskProportion" ]], rv144.lowRiskProportion = other.parameters[[ "rv144.lowRiskProportion" ]], hvtn702.log10lambda = x[ 1 ], hvtn702.high.risk.multiplier = x[ 2 ], hvtn702.highRiskProportion = x[ 3 ], hvtn702.lowRiskProportion = x[ 4 ] )
        }
    }

    make.epsilon.optim.fn <- function ( .f.epsilon.abc, .target.stats, .target.stat.stdevs = rep( 1, length( .target.stats ) ) ) {
        function( x ) {
            .stats <- .f.epsilon.abc( x );
            .stats.matrix <- matrix( .stats, nrow = 1 );
            c( dist = calculate.abc.dist( .stats.matrix, .target.stats, ifelse( is.na( .target.stat.stdevs ), 1, .target.stat.stdevs ) ) )
        }
    } # make.epsilon.optim.fn (..)

    make.rv144.optim.fn <- function ( .f.rv144.abc, .target.stats, .target.stat.stdevs = rep( 1, length( .target.stats ) ) ) {
        function( x ) {
            .stats <- .f.rv144.abc( x );
            .stats.matrix <- matrix( .stats, nrow = 1 );
            c( dist = calculate.abc.dist( .stats.matrix, .target.stats, ifelse( is.na( .target.stat.stdevs ), 1, .target.stat.stdevs ) ) )
        }
    } # make.rv144.optim.fn (..)

    make.hvtn702.optim.fn <- function ( .f.hvtn702.abc, .target.stats, .target.stat.stdevs = rep( 1, length( .target.stats ) ) ) {
        function( x ) {
            .stats <- .f.hvtn702.abc( x );
            .stats.matrix <- matrix( .stats, nrow = 1 );
            c( dist = calculate.abc.dist( .stats.matrix, .target.stats, ifelse( is.na( .target.stat.stdevs ), 1, .target.stat.stdevs ) ) )
        }
    } # make.hvtn702.optim.fn (..)

    # Filter the samples (perhaps subsetted samples from abc) to keep points up to max.dist.scaled units away, where one unit is eg the closest 5% (units.quantile) of the data.
    compute.trial.specific.candidate.modes.from.subset.of.sampled.points <-
        function (
                  the.fit.bin,
                  trial.parameters.names,
                  trial.target.stats.names,
                  target.stats,
                  target.stat.scale.units,
                  units.quantile,
                  max.dist.scaled,
                  pdfCluster.hmult,
                  Tukey.IQR.multiplier
                  ) {
        the.fit.trial.dist <-
            calculate.abc.dist( the.fit.bin$stats[ , trial.target.stats.names ], as.numeric( target.stats[ trial.target.stats.names ] ), ifelse( is.na( the.fit.bin$stats_normalization[ trial.target.stats.names ] ), 1, the.fit.bin$stats_normalization[ trial.target.stats.names ] ) * target.stat.scale.units[ trial.target.stats.names ] );

        # Scale distance by units defined by the max distance of the closest units.quantile fraction of the points.
        trial.dist.units <- quantile( the.fit.trial.dist, probs = units.quantile )
        the.fit.trial.dist.scaled <- the.fit.trial.dist / trial.dist.units;
    
        abc.trial.keep.sim <- the.fit.trial.dist.scaled < max.dist.scaled;
        cat( paste( "Keeping ", sum( abc.trial.keep.sim ), " trial-specific parameter sets because they are within ", max.dist.scaled, " scaled units on the trial-specific distance measure, where one unit has ", sprintf( "%0.2f", 100*units.quantile ), "% of the original ", nrow( the.fit.bin$param ), " draws in this epsilon bin.", sep = "" ), fill = TRUE );
    
        the.fit.trial <- the.fit.bin;
        the.fit.trial$param <- the.fit.trial$param[ abc.trial.keep.sim, c( "epsilon", trial.parameters.names ), drop = FALSE ];
        the.fit.trial$stats <- the.fit.trial$stats[ abc.trial.keep.sim, trial.target.stats.names, drop = FALSE ];
        the.fit.trial$weights <- the.fit.trial$weights[ abc.trial.keep.sim ];
        the.fit.trial.dist <- the.fit.trial.dist[ abc.trial.keep.sim ];

        cl.trial <- suppressWarnings( pdfCluster( the.fit.trial$param[ , trial.parameters.names, drop = FALSE ], bwtype="adaptive", hmult=pdfCluster.hmult, n.grid=nrow( the.fit.trial$param ) ) );
        trial.cluster.numbers <- groups( cl.trial );
        # Helpful when debugging - run through the next line.. it'll print when done
        cat( "trial cluster sizes:", fill = TRUE );
        print( table( trial.cluster.numbers ) );

        medians.by.trial.cluster <- sapply( 1:max( trial.cluster.numbers, na.rm = TRUE ), function( .cluster ) { apply( the.fit.trial$param[ !is.na( trial.cluster.numbers ) & ( trial.cluster.numbers == .cluster ), , drop = FALSE ], 2, median ) } );
        mins.by.trial.cluster <- sapply( 1:max( trial.cluster.numbers, na.rm = TRUE ), function( .cluster ) { apply( the.fit.trial$param[ !is.na( trial.cluster.numbers ) & ( trial.cluster.numbers == .cluster ), , drop = FALSE ], 2, min ) } );
        maxs.by.trial.cluster <- sapply( 1:max( trial.cluster.numbers, na.rm = TRUE ), function( .cluster ) { apply( the.fit.trial$param[ !is.na( trial.cluster.numbers ) & ( trial.cluster.numbers == .cluster ), , drop = FALSE ], 2, max ) } );
        IQRs.by.trial.cluster <- sapply( 1:max( trial.cluster.numbers, na.rm = TRUE ), function( .cluster ) { apply( the.fit.trial$param[ !is.na( trial.cluster.numbers ) & ( trial.cluster.numbers == .cluster ), , drop = FALSE ], 2, function( .col ) { diff( quantile( .col, c( 0.25, 0.75 ) ) ) } ) } );

        # Cluster-specific best sample is at this index, for each cluster:
        trial.cluster.member.minimizing.dist <- sapply( 1:max( trial.cluster.numbers, na.rm = TRUE ), function( .cluster.number ) { .cluster.indices <- which( !is.na( trial.cluster.numbers ) & ( trial.cluster.numbers == .cluster.number ) ); return( .cluster.indices[ which.min( the.fit.trial.dist[ .cluster.indices ] ) ] ); } );
    
        # instead of around the median or mode, center it around the cluster minimizer.
        low.Tukey.whisker.bound.by.trial.cluster <- sapply( 1:ncol( medians.by.trial.cluster ), function( .cluster ) { .tukey.low.whisker.candidate.values <- the.fit.trial$param[ trial.cluster.member.minimizing.dist[[ .cluster ]],  ] - ( Tukey.IQR.multiplier * IQRs.by.trial.cluster[ , .cluster ] ); ifelse( .tukey.low.whisker.candidate.values < mins.by.trial.cluster[ , .cluster ], mins.by.trial.cluster[ , .cluster ], .tukey.low.whisker.candidate.values ) } );
        high.Tukey.whisker.bound.by.trial.cluster <- sapply( 1:ncol( medians.by.trial.cluster ), function( .cluster ) { .tukey.low.whisker.candidate.values <- the.fit.trial$param[ trial.cluster.member.minimizing.dist[[ .cluster ]],  ] + ( Tukey.IQR.multiplier * IQRs.by.trial.cluster[ , .cluster ] ); ifelse( .tukey.low.whisker.candidate.values < mins.by.trial.cluster[ , .cluster ], mins.by.trial.cluster[ , .cluster ], .tukey.low.whisker.candidate.values ) } );
        return( list( low = low.Tukey.whisker.bound.by.trial.cluster, high = high.Tukey.whisker.bound.by.trial.cluster ) );
    } # compute.trial.specific.candidate.modes.from.subset.of.sampled.points (..)



runSim_rv144.hvtn702 <- function( reac = c( "numExecution" = 10 ) ) { # Use numExecution >>1000 for best results.
    stopifnot( all( c( "numExecution" ) %in% names( reac ) ) );

    ## Number of parameters to optimize (3, 4, or 5).
    num.params <- 9; # This is now fixed at 9 by definition, since we are exploraing a splace with a total of 9 = 1, 4, 4 params.

    num.sims <- unname( reac[ "numExecution" ] );

    ######################################################################
    ##### PARAMETERS / MAGIC #s -- MANY COULD BE EXPOSED FOR TUNING.
    ## TODO : REMOVE?
    #abc.keep.num.sims <- 1000; 
    if( num.sims <= 1000 ) {
        # Number of samples to run through the clustering and optimizing steps.
        abc.keep.num.sims <- floor( num.sims / 4 ); # MAGIC #
    } else {
        # Number of samples to run through the clustering and optimizing steps.
        abc.keep.num.sims <- 2500; # MAGIC #
    }
    stopifnot( num.sims > abc.keep.num.sims ); # It won't work to cluster uniformly drawn points. You first have to filter them by keeping those nearest the target.

    # MAGIC #s for the targets.
    rv144.VE <- 31.0;
    rv144.placeboIncidence <- 0.14;
    hvtn702.VE <- 0.0;
    hvtn702.placeboIncidence <- 3.3;

    rv144.placebo.incidence.target <- rv144.placeboIncidence; # incidence per 100 person years
    rv144.VE.target = rv144.VE; # cumulative VE observed by the end of the trial as percentage eg 31 for rv144
    hvtn702.placebo.incidence.target <- hvtn702.placeboIncidence; # incidence per 100 person years
    hvtn702.VE.target = hvtn702.VE; # cumulative VE observed by the end of the trial

    # For the abc and optimization we need a way to compute distances,
    # for which we use the abc default function which is Euclidean
    # distance after scaling each dimension; you can scale it by the
    # standard deviation as abc does, but for optimization we need a
    # better way to balance these, so instead of the observed SD we
    # explicitly weight the scales by setting the scale units such
    # that after scaling the dimensions by these corresponding
    # scale.unit values, the Euclidean distance of those scaled values
    # is approximately correctly weighting the multiple dimensions'
    # contributions. For example if you are finding the optimized
    # modes are very close in one dimension and not in another
    # dimension, you can modify these scales to help balance that
    # distance cost across dimensions better.
    placebo.incidence.target.scale.units <- 1; # MAGIC # (can tune the ratio placebo.incidence.target.scale.units:VE.target.scale.units to alter the relative weight in the distance function).
    VE.target.scale.units <- 0.1; # MAGIC #

    # MAGIC #s tune to modify bounds for the parameters; this could eliminate some degenerate parameter combos (degenerate in the sense that they make the model effectively a two-component model).
    high.risk.multiplier.max <- 50;
    log10lambda.min <- -7;
    log10lambda.max <- -3;

    # MAGIC # but this doesn't really seem to be a tunable parameter, but I basically just set it arbitrarily to a value that I think is enough to get some change; it might need tuning in future.
    smallest.discernable.amount <- 5E-4; # determined by trial and error this is the smallest amount you can change the parameters from 0 or 1 for it to register a difference from those extremes, eg to avoid NaN and Inf

    # We keep at least a fraction of ( abc.keep.num.sims / num.sims ) sims within each cluster, but we actually go a little further to help ensure we are not getting just isolated peaks but also some of the points nearby them to help estimate the IQRs of the clusters for the subsequent search. This is the scale factor by which we multiple the distance; so for example if we are keeping 2.5% of the samples, and the 0.025 quantile of trial-specific distances in a bin for rv144 is 0.67, then we would actually possibly keep more than 2.5% of the samples because we would keep all that samples fall within max.dist.scaled (which should be >= 1) times that distance (eg 1.05 * 0.67). Note that we separately keep samples in this manner for each trial for trial-specific clustering of the bin-specific samples.
    max.dist.scaled <- 1.05; # MAGIC #

    # Once we identify modes we determine their compatibility across studies by whether the epsilons are near enough to each other, based on the overlap or not of their search windows, which are defined by using the logic of identifying outliers in a box plot. That is we walk out from the mode some number of IQR-defined units and define the search window as anything within that zone. For Tukey-style boxplots the number is 1.5 IQRs from the median defines an outlier. Here we can tune this multiplier (Tukey.IQR.multiplier) to change the window width that we use to define the search space for optimization and also for this overlap detection.
    Tukey.IQR.multiplier <- 0.5; # MAGIC # governing the window size we use when optimizing and when merging trial-specific candidate optima (this is applied to epsilon only for merging).

    # If there are candidates from both studies but no bins overlap, do we force overlap of closest pair?
    ensure.overlap.in.epsilon.bins <- TRUE; # MAGIC # (keep TRUE to ensure candidates are explored in farther-distance epsilon bins.)

    # What to do with overlapping bounds ? Use the intersection of the tukey whisker bounds? Otherwise, use the union.
    merge.bounds.strategy.intersect <- FALSE; # MAGIC # that could be explored in conjuction with eg Tukey.IQR.multiplier -- when merging epsilon windows do we take the outer or inner bounds? We keep this FALSE for now indicating we use the UNION (so, the outer bounds, expanding the epsilon search window to be inclusive of both trials' windows, rather than shrinking it to only where the intervals overlap).

    # This helps when debugging...
    if( abc.keep.num.sims < 1000 ) {
        num.epsilon.bins <- 2; # MAGIC #
    } else {
        num.epsilon.bins <- 10; # MAGIC # (tune this to control how many total candidates we explore, since we ensure we keep at least one candidate per epsilon bin (if ensure.overlap.in.epsilon.bins is TRUE).
    }

    min.points.for.clustering <- 6; # MAGIC # >= 6, from "QH6214 qhull input error: not enough points(2) to construct initial simplex (need 6)"

    pdfCluster.hmult <- 1.05; # MAGIC #, tweaked it to get pdfCluster to run without crashing with the error message suggesting increasing n.grid -- even with max n.grid, hmult has to be just above 1, it seems. If it is too high, the clusters merge into one. Tune this (> 1 to prevent crashing) higher to get the clusters to merge more.

    be.verbose <- TRUE; # MAGIC # (just governs output printed to screen)
    ######################################################################


    ######################################################################
    # Code below here has no additional tunable parameters...
    ######################################################################
    target.stats <-
        data.frame( rv144.VE.target, rv144.placebo.incidence.target, hvtn702.VE.target, hvtn702.placebo.incidence.target );        
    rv144.target.stats.names <- c( "rv144.VE.target", "rv144.placebo.incidence.target" );
    hvtn702.target.stats.names <- c( "hvtn702.VE.target", "hvtn702.placebo.incidence.target" );

    ## Note that we target the cumulative incidence at the end of the trial.
    target.stat.scale.units <- c( "rv144.VE.target" = VE.target.scale.units, "rv144.placebo.incidence.target" = placebo.incidence.target.scale.units, "hvtn702.VE.target" = VE.target.scale.units, "hvtn702.placebo.incidence.target" = placebo.incidence.target.scale.units );

    # Specify bounds for the initial conditions; these are used as priors.
    # In order of "x" in the above function.
    bounds <- list(
                   epsilon = c(smallest.discernable.amount, 1-smallest.discernable.amount), # This is just the full range 0 to 1
                   rv144.log10lambda = c(log10lambda.min, log10lambda.max),
                   rv144.high.risk.multiplier = c(1, high.risk.multiplier.max), # risk multiplier for high risk group
                   rv144.highRiskProportion = c(smallest.discernable.amount, 0.5-smallest.discernable.amount), # This is the half range 0 to 0.5
                   rv144.lowRiskProportion = c(smallest.discernable.amount, 1-smallest.discernable.amount), # This is modified to avoid NaN/Inf
                   hvtn702.log10lambda = c(log10lambda.min, log10lambda.max),
                   hvtn702.high.risk.multiplier = c(1, high.risk.multiplier.max), # risk multiplier for high risk group
                   hvtn702.highRiskProportion = c(smallest.discernable.amount, 0.5-smallest.discernable.amount), # This is the half range 0 to 0.5
                   hvtn702.lowRiskProportion = c(smallest.discernable.amount, 1-smallest.discernable.amount) # This is modified to avoid NaN/Inf
                   );            
    stopifnot( all( names( bounds ) == all.parameters ) );

    # In order of "x" in the above function.
    priors  <- lapply( bounds, function( .bounds ) { c( "unif", unlist( .bounds ) ) } );

    ## PHASE 1: Find candidate complete 9-parameter starting places constructed from merging epsion-bin-specific, trial-specific local optima that share common ranges of epsilon across the two trials. Later (in phase 2) we will find 9-parameter local optima near each of these candidate starting places.
    ## First we draw num.sims draws, then for each candidate epsilon bin we compute the study-specific distances for points falling in that bin (computed using just that study's two target stats), and cluster all points within 1.05-fold (see max.dist.scaled) of the furthest of the set of the top (abc.keep.num.sims/num.sims) (eg 2.5%) of samples within that bin, to ensure a minimum number of points and the extra points up to max.dist.scaled ensures that we keep additional points to help flesh out the contours of the space just below these peaks. This might possibly help with the clustering but it's not entirely clear yet how sensitive things are to max.dist.scaled.

    # First we just draw num.sims draws from the priors, independently. Keep everything drawn. Compute the 4-parameter target stats from runs with these 9-parameter random starting places, and the standard deviations of these 4-parameter stats (which we use to scale the distance function on the target stats for balancing the optimization evenly across the parameters, see below).
    .f.abc <- function( x ) {
        run.and.compute.run.stats.rv144.hvtn702( epsilon = x[ 1 ], rv144.log10lambda = x[ 2 ], rv144.high.risk.multiplier = x[ 3 ], rv144.highRiskProportion = x[ 4 ], rv144.lowRiskProportion = x[ 5 ], hvtn702.log10lambda = x[ 6 ], hvtn702.high.risk.multiplier = x[ 7 ], hvtn702.highRiskProportion = x[ 8 ], hvtn702.lowRiskProportion = x[ 9 ] )
    }
    fit.rej <- ABC_rejection(model = .f.abc,
                         prior = priors,
                         nb_simul = num.sims,
                         summary_stat_target = as.numeric( target.stats ),
                         tol = 1.0, # Keep all of them for the moment; see below.
                         progress_bar = TRUE)
    colnames( fit.rej$param ) <- all.parameters;
    colnames( fit.rej$stats ) <- names( target.stats );
    names( fit.rej$stats_normalization ) <- names( target.stats );

    ## Great, we will now consider bins of epsilon (the only parameter that is shared across the two studies) and within each bin we will cluster the non-epsilon parameters for each study and construct a set of 9-parameter candidate starting places based on study-specific modes that share approximately common epsilon parameters across the studies. For example if there is a mode at around epsilon = 0.5 for both studies, we want to construct a 9-parameter starting place with epsilon = 0.5 and the study-specific maximizing parameters for the non-epsilon parameters when epsilon is 0.5.

   # Note that we use a sliding window approach so there's actually 2*num.epsilon.bins-1 total bins, and we may find the same modes multiple times. This is to balance focusing on relevant values of the other parameters while clustering the non-epsilon study-specific parameters by conditioning approximately on epsilon, and might need be tuned to consider different bin sizes (bin sizes are 1/num.epsilon.bins and overlap halfway through).

    ## Strategy within each bin: first compute the candidate parameter windows by clustering the points separately for each trial (unbounded, centered at the optima), and merging them if they are compatible, then compute the midpoints (prior to bounding) as the intial points for the optization, then finally apply global bounds to all of them (the ends and the midpoints jic). Note the order of when we apply the global bounds - it is just a hack to ensure that the midpoints represent the sampled optima even if the computed endpoints are out of bounds, lazily not keeping track of the sampled optima themselves, because of the symmetry of the way we compute the windows centered on those optima and can recover them (but not if we were to apply the bounds to the windows first, see?).

    # Walk up the epsilon range in overlapping windows (2*num.epsilon.bins - 1 of them). Within each window, find new candidate parameter sets by taking candidate sets from each trial separately and making complete parameter sets whenever the epsilon ranges of any two candidate parameter sets overlap. The width of the window is tunable using the Tukey.IQR.multiplier defined above (units are IQRs; center is the minimal sampled point within the cluster note this procedure does not guarantee anything and must be considered very hacky - one example is that the minimal point might be one of the clusters' outliers, this still will center the window at that point, so it is possibble technically for that point to be the only sampled point in the window - we ignore the sampled points henceforth anyway but I just want to be clear how hacky this is; sometimes candidate clusters are ignored because the epsilon windows do not overlap and this might turn out to be in need of tuning. Tuning parameters are highlighted as MAGIC #s above, but particularly note that the clustering is more likely to lump things together if you increase pdfMult, and epsilon windows are more likely to overlap across studies (and therefore complete a potential two-study parameter set) if you increase Tukey.IQR.multiplier).
    candidate.parameter.sets.low <- matrix( NA, nrow = 0, ncol = length( all.parameters ) );
    candidate.parameter.sets.high <- matrix( NA, nrow = 0, ncol = length( all.parameters ) );
    for( epsilon.bin in 1:(2*num.epsilon.bins - 1) ) {
        bin.min.epsilon <- max( 0, ( epsilon.bin - 1 )/(2*num.epsilon.bins) );
        bin.max.epsilon <- min( 1, ( epsilon.bin + 1 )/(2*num.epsilon.bins) );
        if( be.verbose ) {
            cat( paste( "epsilon.bin =", epsilon.bin ), fill = TRUE );
            cat( paste( "\tEpsilon range: from", bin.min.epsilon, "to", bin.max.epsilon ), fill = TRUE );
        }

        ### Determine the subset of sampled points that fall in this epsilon bin:
        sample.is.in.bin <-
            ( fit.rej$param[ , "epsilon" ] >= bin.min.epsilon ) & ( fit.rej$param[ , "epsilon" ] < bin.max.epsilon );
        if( be.verbose ) {
            cat( paste( "sum( sample.is.in.bin ) =", sum( sample.is.in.bin ) ), fill = TRUE );
        }
        if( sum( sample.is.in.bin ) < min.points.for.clustering ) {
            next; # Move on.
        }

        fit.rej.bin <- fit.rej;
        fit.rej.bin$param <- fit.rej$param[ sample.is.in.bin, , drop = FALSE ];
        fit.rej.bin$stats <- fit.rej$stats[ sample.is.in.bin, , drop = FALSE ];
        fit.rej.bin$weights <- fit.rej$weights[ sample.is.in.bin ];
    
        ### rv144
        rv144.Tukey.whisker.bounds.by.trial.cluster <-
            compute.trial.specific.candidate.modes.from.subset.of.sampled.points( fit.rej.bin, rv144.parameters, rv144.target.stats.names, target.stats, target.stat.scale.units, units.quantile = ( abc.keep.num.sims / num.sims ), max.dist.scaled = max.dist.scaled, pdfCluster.hmult = pdfCluster.hmult, Tukey.IQR.multiplier = Tukey.IQR.multiplier );
        low.Tukey.whisker.bound.by.rv144.cluster <- rv144.Tukey.whisker.bounds.by.trial.cluster[[ "low" ]];
        high.Tukey.whisker.bound.by.rv144.cluster <- rv144.Tukey.whisker.bounds.by.trial.cluster[[ "high" ]];

        ### hvtn702
        hvtn702.Tukey.whisker.bounds.by.trial.cluster <-
            compute.trial.specific.candidate.modes.from.subset.of.sampled.points( fit.rej.bin, hvtn702.parameters, hvtn702.target.stats.names, target.stats, target.stat.scale.units, units.quantile = ( abc.keep.num.sims / num.sims ), max.dist.scaled = max.dist.scaled, pdfCluster.hmult = pdfCluster.hmult, Tukey.IQR.multiplier = Tukey.IQR.multiplier );
        low.Tukey.whisker.bound.by.hvtn702.cluster <- hvtn702.Tukey.whisker.bounds.by.trial.cluster[[ "low" ]];
        high.Tukey.whisker.bound.by.hvtn702.cluster <- hvtn702.Tukey.whisker.bounds.by.trial.cluster[[ "high" ]];

        # Once we find a pair with overlapping epsilon windows (one from each trial) we use this to add a new set of bounds to the candidate.parameter.sets.high and candidate.parameter.sets.low values, which are altered by this (so it is not a pure function! it has side effects! We define this here for code reuse - and it is defined here intentionally so it can access the values low.Tukey.whisker.bound.by.hvtn702.cluster etc.
        add.candidate.complete.parameter.set <- function ( rv144.cluster.i, hvtn702.cluster.j ) {
            if( merge.bounds.strategy.intersect ) {
                low.epsilon.bound <-
                    max(
                        low.Tukey.whisker.bound.by.rv144.cluster[ "epsilon", rv144.cluster.i ],
                        low.Tukey.whisker.bound.by.hvtn702.cluster[ "epsilon", hvtn702.cluster.j ]
                        );
                high.epsilon.bound <-
                    min(
                        high.Tukey.whisker.bound.by.rv144.cluster[ "epsilon", rv144.cluster.i ],
                        high.Tukey.whisker.bound.by.hvtn702.cluster[ "epsilon", hvtn702.cluster.j ]
                        );
            } else { # if merge.bounds.strategy.intersect .. else ..
                # use the strategy "union" instead of "intersection"
                low.epsilon.bound <-
                    min(
                        low.Tukey.whisker.bound.by.rv144.cluster[ "epsilon", rv144.cluster.i ],
                        low.Tukey.whisker.bound.by.hvtn702.cluster[ "epsilon", hvtn702.cluster.j ]
                        );
                high.epsilon.bound <-
                    max(
                        high.Tukey.whisker.bound.by.rv144.cluster[ "epsilon", rv144.cluster.i ],
                        high.Tukey.whisker.bound.by.hvtn702.cluster[ "epsilon", hvtn702.cluster.j ]
                        );
            } # End if merge.bounds.strategy.intersect .. else ..
            stopifnot( high.epsilon.bound > low.epsilon.bound );
            candidate.parameter.sets.low <<-
                rbind( candidate.parameter.sets.low,
                      c( "epsilon" = low.epsilon.bound,
                        low.Tukey.whisker.bound.by.rv144.cluster[ rv144.parameters, rv144.cluster.i ],
                        low.Tukey.whisker.bound.by.hvtn702.cluster[ hvtn702.parameters, hvtn702.cluster.j ]
                        ) );
            candidate.parameter.sets.high <<-
                rbind( candidate.parameter.sets.high,
                      c( "epsilon" = high.epsilon.bound,
                        high.Tukey.whisker.bound.by.rv144.cluster[ rv144.parameters, rv144.cluster.i ],
                        high.Tukey.whisker.bound.by.hvtn702.cluster[ hvtn702.parameters, hvtn702.cluster.j ]
                        ) );
 
            return( NULL );
        } # add.candidate.complete.parameter.set (..)

        # Returns TRUE iff the two clusters have compatible epsilon values. Determined by overlap of epsilon windows.
        parameters.are.compatible <- function ( rv144.cluster.i, hvtn702.cluster.j ) {
            ( high.Tukey.whisker.bound.by.rv144.cluster[ "epsilon", rv144.cluster.i ] < low.Tukey.whisker.bound.by.hvtn702.cluster[ "epsilon", hvtn702.cluster.j ] ) || ( high.Tukey.whisker.bound.by.hvtn702.cluster[ "epsilon", hvtn702.cluster.j ] < low.Tukey.whisker.bound.by.rv144.cluster[ "epsilon", rv144.cluster.i ] )
        } # parameters.are.compatible (..)

        # Returns the absolute minimum distance between epsilon window bounds
        compute.distance.of.non.overlapping.pairs <- function ( rv144.cluster.i, hvtn702.cluster.j ) {
            min( abs( low.Tukey.whisker.bound.by.hvtn702.cluster[ "epsilon", hvtn702.cluster.j ] - high.Tukey.whisker.bound.by.rv144.cluster[ "epsilon", rv144.cluster.i ] ), abs( low.Tukey.whisker.bound.by.rv144.cluster[ "epsilon", rv144.cluster.i ] - high.Tukey.whisker.bound.by.hvtn702.cluster[ "epsilon", hvtn702.cluster.j ] ) )
        } # compute.distance.of.non.overlapping.pairs (..)

        closest.pair.rv144.cluster.i <- NA;
        closest.pair.hvtn702.cluster.j <- NA;
        closest.pair.dist <- NA;
        compatible.pair.found <- FALSE;
        for( rv144.cluster.i in 1:ncol( low.Tukey.whisker.bound.by.rv144.cluster ) ) {
            if( be.verbose ) {
                cat( paste( "rv144.cluster.i = ", rv144.cluster.i ), fill = TRUE );
            }
            for( hvtn702.cluster.j in 1:ncol( low.Tukey.whisker.bound.by.hvtn702.cluster ) ) {
                if( be.verbose ) {
                    cat( paste( "hvtn702.cluster.j = ", hvtn702.cluster.j ), fill = TRUE );
                }
                # Is it incompatible?
                if( parameters.are.compatible( rv144.cluster.i, hvtn702.cluster.j ) ) {
                    # The epsilon windows do not overlap. Move on.
                    if( be.verbose ) {
                        cat( paste( "Epsilon windows for rv144 cluster ", rv144.cluster.i, " and hvtn702 cluster ", hvtn702.cluster.j, " in epsilon bin ", epsilon.bin, " do not overlap.", sep = "" ), fill = TRUE );
                    }
                    # But if they are the closest non-overlapping pair, keep note.
                    .dist <- compute.distance.of.non.overlapping.pairs( rv144.cluster.i, hvtn702.cluster.j );
                    if( is.na( closest.pair.dist ) || ( .dist < closest.pair.dist ) ) {
                        if( be.verbose ) {
                            cat( paste( "Found new closest non-overlapping pair: distance is", .dist ), fill = TRUE  );
                        }
                        closest.pair.dist <- .dist;
                        closest.pair.rv144.cluster.i <- rv144.cluster.i;
                        closest.pair.hvtn702.cluster.j <- hvtn702.cluster.j;
                    }
                    next;
                }
                if( be.verbose ) {
                    cat( paste( "Epsilon windows for rv144 cluster ", rv144.cluster.i, " and hvtn702 cluster ", hvtn702.cluster.j, " in epsilon bin ", epsilon.bin, " ARE COMPATIBLE.", sep = "" ), fill = TRUE );
                }
                compatible.pair.found <- TRUE;
                add.candidate.complete.parameter.set( rv144.cluster.i, hvtn702.cluster.j );
            } # End foreach hvtn702.cluster.j
        } # End foreach rv144.cluster.i
        if( !compatible.pair.found && ensure.overlap.in.epsilon.bins ) {
            # Force a pairing even though no pairs overlapped.
            rv144.cluster.i <- closest.pair.rv144.cluster.i;
            hvtn702.cluster.j <- closest.pair.hvtn702.cluster.j;
            ## TODO: Make this a function instead of copied code. It is verbatim copied from above.
            if( be.verbose ) {
                cat( paste( "Epsilon windows for rv144 cluster ", rv144.cluster.i, " and hvtn702 cluster ", hvtn702.cluster.j, " in epsilon bin ", epsilon.bin, " ARE NOT COMPATIBLE BUT WE WILL PAIR THEM AS CLOSEST NON-COMPATIBLE PAIR.", sep = "" ), fill = TRUE );
            }
            add.candidate.complete.parameter.set( rv144.cluster.i, hvtn702.cluster.j );
        } # End if we need to force a pairing.
    } # End foreach epsilon.bin
    cat( paste( "There are", nrow( candidate.parameter.sets.high ), "candidate parameter sets." ), fill = TRUE );

    ## TODO: Filter/merge the overlapping candidates so they are not redundant.

    ### PHASE 1.5: Apply bounds to the computed parameter set search windows and compute initial values.

    ## Compute midpoints of candidates as medians of bounds; these
    ## will be the starting places for the optimizations below, but
    ## first we need to ensure everything is within the specified
    ## bounds.
    candidate.parameter.sets.midpoint <-
        candidate.parameter.sets.low + ( candidate.parameter.sets.high - candidate.parameter.sets.low ) / 2;

    ## Ensure global bounds
    bounds.low <- sapply( bounds, function ( .lower.and.higher ) { .lower.and.higher[ 1 ] } );
    bounds.high <- sapply( bounds, function ( .lower.and.higher ) { .lower.and.higher[ 2 ] } );
    candidate.parameter.sets.low.bounded <- t( apply( candidate.parameter.sets.low, 1, function( .candidate.parameter.set.low ) { ifelse( .candidate.parameter.set.low < bounds.low, bounds.low, .candidate.parameter.set.low ) } ) );
    candidate.parameter.sets.high.bounded <- t( apply( candidate.parameter.sets.high, 1, function( .candidate.parameter.set.high ) { ifelse( .candidate.parameter.set.high > bounds.high, bounds.high, .candidate.parameter.set.high ) } ) );
    candidate.parameter.sets.midpoint.bounded <- t( apply( candidate.parameter.sets.midpoint, 1, function( .candidate.parameter.set.midpoint ) { ifelse( .candidate.parameter.set.midpoint > bounds.high, bounds.high, ifelse( .candidate.parameter.set.midpoint < bounds.low, bounds.low, .candidate.parameter.set.midpoint ) ) } ) );

    ### PHASE 2:  from these candidate starting places constructed from epsion-bin-specific, trial-specific local optima, find modes in the 9-parameter space by conditional optimization.

    optima.by.candidate <- sapply( 1:nrow( candidate.parameter.sets.midpoint.bounded ), function( .candidate ) {
        print( .candidate );
        .lower.bounds <- candidate.parameter.sets.low.bounded[ .candidate, ];
        .upper.bounds <- candidate.parameter.sets.high.bounded[ .candidate, ];
        .midpoint <- candidate.parameter.sets.midpoint.bounded[ .candidate, ];
        print( .midpoint );
        .iterative.result <- optimize.iteratively( .midpoint, lower = .lower.bounds, upper = .upper.bounds, current.value = NULL, be.verbose = TRUE );
        # current.parameters <- .iterative.result[ all.parameters ];
        # current.stats <- .iterative.result[[ setdiff( names( .iterative.result ), all.parameters, "dist" ) ]];
        # current.value <- .iterative.result[[ "dist" ]];
        return( .iterative.result );
    } );
    optima.by.candidate.sorted <- optima.by.candidate[ , order( as.numeric( optima.by.candidate[ "dist", ] ) ), drop = FALSE ];
        
    return( list( fit = fit.rej, priors = priors, bounds = bounds, target.stats = target.stats, fn = .f.abc, sampled.modes = optima.by.candidate.sorted ) );
} # runSim_rv144.hvtn702 (..)

### ERE I AM testing...
the.seed <- 98103;
# To test replicability of the identified modes, uncomment this:
# set.seed( the.seed ); the.seed <- floor( runif( 1, max = 1E5 ) );
# num.sims <- 1000; # Fast for debugging.
num.sims <- 10000; # For reals.

set.seed( the.seed );

## Ok, I think that we need to first determine the reasonable epsilon values from the RV144-like setting, then ask could you get zero (or very low) VE with the same epsilon by increasing only lambda and risk values, with placebo incidence going from an rv144-like value (~0.14) to a 702-like value (~3.3).
## rv144 placebo incidence in the prespecified analysis (MITT) cohort was 0.1397 per 100 person years, 100*(74/52985) from N Engl J Med 2009; 361:2209-2220 DOI: 10.1056/NEJMoa0908492 December 3, 2009 (https://www.nejm.org/doi/full/10.1056/nejmoa0908492):
# "HIV-1 infection was diagnosed in 132 subjects (56 in the vaccine group and 76 in the placebo group) during 52,985 person-years of follow-up in the intention-to-treat analysis, in 86 subjects (36 in the vaccine group and 50 in the placebo group) during 36,720 person-years of follow-up in the per-protocol analysis, and in 125 subjects (51 in the vaccine group and 74 in the placebo group) during 52,985 person-years of follow-up in the modified intention-to-treat analysis. One subject in the placebo group who was identified by hospital record as being seropositive for HIV after dying from Pneumocystis jirovecii pneumonia was included in the analysis before the unblinding of the study. This diagnosis of HIV-1 infection was the only one that occurred outside planned procedures."

## hvtn 702 placebo incidence was 3.3 per 100 person-years (95% CI, 2.8 to 3.9), from n engl j med 384;12 nejm.org March 25, 2021 (https://www.nejm.org/doi/pdf/10.1056/NEJMoa2031499), page 1092:
# "During the first 24 months of follow-up, 138 HIV-1 infections occurred in the vaccine group and 133 in the placebo group, for an estimated incidence rate of 3.4 per 100 person-years (95% confidence interval [CI], 2.8 to 4.0) and 3.3 per 100 person-years (95% CI, 2.8 to 3.9), respectively (hazard ratio, 1.02; 95% CI, 0.81 to 1.30; P=0.84) (Fig. 1A and Table 2). The incidence of HIV-1 infection was similar in the vaccine group and the placebo group in secondary analyses during 36 months of follow-up (hazard ratio, 1.05; 95% CI, 0.83 to 1.31), in the month 6.5 at-risk cohort between 6.5 months and 24 months (hazard ratio, 1.15; 95% CI, 0.84 to 1.58), and in the perprotocol cohort, as well as in other secondary analyses (Figs. S3 through S9)."

# .sim <- runSim_rv144.hvtn702( reac = c( "numExecution" = num.sims ) );

######
## Some plotting. Run manually. See above.
if( FALSE ) {
    the.sim <- .sim;

    fit <- the.sim$fit
    bounds <- the.sim$bounds;
    
    # Contour plotting
    fit.dist <-
        calculate.abc.dist( fit$stats, as.numeric( target.stats ), ifelse( is.na( fit$stats_normalization ), 1, fit$stats_normalization ) * target.stat.scale.units );

    dist.units <- quantile( fit.dist, probs = 0.025 ) # top 2.5%
    dist.units.for.contour.plot <- 1.0; # Show all the data within 1 unit (so, the top 2.5%).
    fit.dist.scaled <- fit.dist / dist.units;
    
    abc.keep.sim <- fit.dist.scaled < dist.units.for.contour.plot;
    fit.filtered <- fit;
    fit.filtered$param <- fit$param[ abc.keep.sim, , drop = FALSE ];
    fit.filtered$stats <- fit$stats[ abc.keep.sim, , drop = FALSE ];
    fit.filtered$weights <- fit$weights[ abc.keep.sim ];
    fit.filtered.dist <- fit.dist[ abc.keep.sim ];

    # These are the original bounds, separated for use in optim:
    lower.bounds <- sapply( bounds, function ( .bounds.for.x ) { .bounds.for.x[ 1 ] } );
    upper.bounds <- sapply( bounds, function ( .bounds.for.x ) { .bounds.for.x[ 2 ] } );

    plot(density(fit.filtered$param[, 1], from = lower.bounds[1],  to = upper.bounds[1]),
         main = "epsilon", 
         xlim = c(lower.bounds[1], upper.bounds[1]),
         #ylim = c(0, 10),
         col=2)
    lines(density(fit.filtered$param[, 1], from = lower.bounds[1],  to = upper.bounds[1]), col = 2)
    # There is not really a target for this parameter. This might be misleading.
    #abline(v = the.sim$target.stats[[ "rv144.VE.target" ]], lty = 2, col = 1)
    legend("topright", legend = c("per-contact VE", "Posterior"),
           lty = c(1, 2), col = 1:2, lwd = 2)
    
    plot(density(fit.filtered$param[, 2], from = lower.bounds[2],  to = upper.bounds[2]),
         main = "RV144 log10lambda", 
         xlim = c(lower.bounds[2], upper.bounds[2]),
         col=2)
    lines(density(fit.filtered$param[, 2], from = lower.bounds[2],  to = upper.bounds[2]), col = 2)
    #abline(v = VE.target, lty = 2, col = 1) # This is a bug, it plots VE target, not lambda
    legend("topright", legend = c("Truth", "Posterior"),
           lty = c(1, 2), col = 1:2, lwd = 2)
    
    plot(density(fit.filtered$param[, 3], from = lower.bounds[3],  to = upper.bounds[3]),
         main = "RV144 high risk multiplier", 
         xlim = c(lower.bounds[3], upper.bounds[3]),
         col=2)
    lines(density(fit.filtered$param[, 3], from = lower.bounds[3],  to = upper.bounds[3]), col = 2)
    #abline(v = VE.target, lty = 2, col = 1) # Another bug
    legend("topright", legend = c("Truth", "Posterior"),
           lty = c(1, 2), col = 1:2, lwd = 2)
    
    ## This is just printing the rv144 contours, see below for 702 contours.
    .rv144.df <- as.data.frame( fit.filtered$param[ , c( "epsilon", "rv144.log10lambda", "rv144.high.risk.multiplier", "rv144.highRiskProportion", "rv144.lowRiskProportion" ) ] );
    names( .rv144.df ) <- c( "epsilon", "log10lambda", "risk", "highProp", "lowPropOfNonHigh" );
    .rv144.df <- cbind( .rv144.df, data.frame( "lowProp" = .rv144.df$lowPropOfNonHigh * ( 1 - .rv144.df$highProp ) ) );
    pdf( "rv144log10lambda_rv144risk.pdf" ); ggplot( .rv144.df, aes( x=log10lambda,y=risk ) ) + geom_point() + stat_density2d_filled() + ggtitle( "RV144 log risk (log10lambda) vs RV144 high-risk FC over normal (risk)" ); dev.off()
    pdf( "epsilon_rv144risk.pdf" ); ggplot( .rv144.df, aes( x=epsilon,y=risk ) ) + geom_point() + stat_density2d_filled() + ggtitle( "Per-contact VE (epsilon) vs RV144 high-risk group FC over normal (risk)" ); dev.off()
    pdf( "epsilon_rv144log10lambda.pdf" ); ggplot( .rv144.df, aes( x=epsilon,y=log10lambda ) ) + geom_point() + stat_density2d_filled() + ggtitle( "Per-contact VE (epsilon) vs RV144 baseline log risk (log10lambda)" ); dev.off()
    pdf( "epsilon_rv144highProp.pdf" ); ggplot( .rv144.df, aes( x=epsilon,y=highProp ) ) + geom_point() + stat_density2d_filled() + ggtitle( "Per-contact VE (epsilon) vs RV144 % of pop at high risk (highProp)" ); dev.off()
    pdf( "epsilon_rv144lowProp.pdf" ); ggplot( .rv144.df, aes( x=epsilon,y=lowProp ) ) + geom_point() + stat_density2d_filled() + ggtitle( "Per-contact VE (epsilon) vs RV144 % of pop at low risk (lowProp)" ); dev.off()
    pdf( "rv144log10lambda_rv144highProp.pdf" ); ggplot( .rv144.df, aes( x=log10lambda,y=highProp ) ) + geom_point() + stat_density2d_filled() + ggtitle( "RV144 baseline log risk (log10lambda) vs % of pop at high risk (highProp)" ); dev.off()
    pdf( "rv144log10lambda_rv144lowProp.pdf" ); ggplot( .rv144.df, aes( x=log10lambda,y=lowProp ) ) + geom_point() + stat_density2d_filled() + ggtitle( "RV144 baseline log risk (log10lambda) vs % of pop at low risk (lowProp)" ); dev.off()
    pdf( "rv144lowProp_rv144highProp.pdf" ); ggplot( .rv144.df, aes( x=lowProp,y=highProp ) ) + geom_point() + stat_density2d_filled() + ggtitle( "RV144 % of pop at low risk (lowProp) vs % of pop at high risk (highProp)" ); dev.off()
    pdf( "rv144risk_rv144highProp.pdf" ); ggplot( .rv144.df, aes( x=risk,y=highProp ) ) + geom_point() + stat_density2d_filled() + ggtitle( "RV144 high-risk FC over normal (risk) vs % of pop at high risk (highProp)" ); dev.off()
    pdf( "rv144risk_rv144lowProp.pdf" ); ggplot( .rv144.df, aes( x=risk,y=lowProp ) ) + geom_point() + stat_density2d_filled() + ggtitle( "RV144 high-risk group FC over normal (risk) vs % of pop at zero risk (lowProp)" ); dev.off()

    ## This is just printing the hvtn702 contours, see above for rv144 contours.
    .hvtn702.df <- as.data.frame( fit.filtered$param[ , c( "epsilon", "hvtn702.log10lambda", "hvtn702.high.risk.multiplier", "hvtn702.highRiskProportion", "hvtn702.lowRiskProportion" ) ] );
    names( .hvtn702.df ) <- c( "epsilon", "log10lambda", "risk", "highProp", "lowPropOfNonHigh" );
    .hvtn702.df <- cbind( .hvtn702.df, data.frame( "lowProp" = .hvtn702.df$lowPropOfNonHigh * ( 1 - .hvtn702.df$highProp ) ) );
    pdf( "hvtn702log10lambda_hvtn702risk.pdf" ); ggplot( .hvtn702.df, aes( x=log10lambda,y=risk ) ) + geom_point() + stat_density2d_filled() + ggtitle( "HVTN702 log risk (log10lambda) vs HVTN702 high-risk FC over normal (risk)" ); dev.off()
    pdf( "epsilon_hvtn702risk.pdf" ); ggplot( .hvtn702.df, aes( x=epsilon,y=risk ) ) + geom_point() + stat_density2d_filled() + ggtitle( "Per-contact VE (epsilon) vs HVTN702 high-risk group FC over normal (risk)" ); dev.off()
    pdf( "epsilon_hvtn702log10lambda.pdf" ); ggplot( .hvtn702.df, aes( x=epsilon,y=log10lambda ) ) + geom_point() + stat_density2d_filled() + ggtitle( "Per-contact VE (epsilon) vs HVTN702 baseline log risk (log10lambda)" ); dev.off()
    pdf( "epsilon_hvtn702highProp.pdf" ); ggplot( .hvtn702.df, aes( x=epsilon,y=highProp ) ) + geom_point() + stat_density2d_filled() + ggtitle( "Per-contact VE (epsilon) vs HVTN702 % of pop at high risk (highProp)" ); dev.off()
    pdf( "epsilon_hvtn702lowProp.pdf" ); ggplot( .hvtn702.df, aes( x=epsilon,y=lowProp ) ) + geom_point() + stat_density2d_filled() + ggtitle( "Per-contact VE (epsilon) vs HVTN702 % of pop at low risk (lowProp)" ); dev.off()
    pdf( "hvtn702log10lambda_hvtn702highProp.pdf" ); ggplot( .hvtn702.df, aes( x=log10lambda,y=highProp ) ) + geom_point() + stat_density2d_filled() + ggtitle( "HVTN702 baseline log risk (log10lambda) vs % of pop at high risk (highProp)" ); dev.off()
    pdf( "hvtn702log10lambda_hvtn702lowProp.pdf" ); ggplot( .hvtn702.df, aes( x=log10lambda,y=lowProp ) ) + geom_point() + stat_density2d_filled() + ggtitle( "HVTN702 baseline log risk (log10lambda) vs % of pop at low risk (lowProp)" ); dev.off()
    pdf( "hvtn702lowProp_hvtn702highProp.pdf" ); ggplot( .hvtn702.df, aes( x=lowProp,y=highProp ) ) + geom_point() + stat_density2d_filled() + ggtitle( "HVTN702 % of pop at low risk (lowProp) vs % of pop at high risk (highProp)" ); dev.off()
    pdf( "hvtn702risk_hvtn702highProp.pdf" ); ggplot( .hvtn702.df, aes( x=risk,y=highProp ) ) + geom_point() + stat_density2d_filled() + ggtitle( "HVTN702 high-risk FC over normal (risk) vs % of pop at high risk (highProp)" ); dev.off()
    pdf( "hvtn702risk_hvtn702lowProp.pdf" ); ggplot( .hvtn702.df, aes( x=risk,y=lowProp ) ) + geom_point() + stat_density2d_filled() + ggtitle( "HVTN702 high-risk group FC over normal (risk) vs % of pop at zero risk (lowProp)" ); dev.off()

    ## Finally, print some cross-study contours. Some thought will have to be given to this. We could maybe color the modes somehow by the epsilon, which I think differs for the two rv144 lambda modes we see here, and it's the only thing that creates non-indpendence so it must also correspond to the multiple two-dimensional modes for the other parameters, too.
    .log10lambda.df <- as.data.frame( fit.filtered$param[ , c( "rv144.log10lambda", "hvtn702.log10lambda" ) ] );
    pdf( "rv144log10lambda_hvtn702log10lambda.pdf" ); ggplot( .log10lambda.df, aes( x=rv144.log10lambda,y=hvtn702.log10lambda ) ) + geom_point() + stat_density2d_filled() + ggtitle( "RV144 vs HVTN702 baseline log risk (log10lambda)" ); dev.off()

    .risk.df <- as.data.frame( fit.filtered$param[ , c( "rv144.high.risk.multiplier", "hvtn702.high.risk.multiplier" ) ] );
    names( .risk.df ) <- c( "rv144.risk", "hvtn702.risk" );
    pdf( "rv144risk_hvtn702risk.pdf" ); ggplot( .risk.df, aes( x=rv144.risk,y=hvtn702.risk ) ) + geom_point() + stat_density2d_filled() + ggtitle( "RV144 vs HVTN702 high-risk FC over normal (risk)" ); dev.off()
    
    .highProp.df <- as.data.frame( fit.filtered$param[ , c( "rv144.highRiskProportion", "hvtn702.highRiskProportion" ) ] );
    names( .highProp.df ) <- c( "rv144.highProp", "hvtn702.highProp" );
    pdf( "rv144highProp_hvtn702highProp.pdf" ); ggplot( .highProp.df, aes( x=rv144.highProp,y=hvtn702.highProp ) ) + geom_point() + stat_density2d_filled() + ggtitle( "RV144 vs HVTN702 % of prop at high risk (highProp)" ); dev.off()
    
    .lowProp.df <- data.frame( "rv144.lowProp" = .rv144.df$lowProp, "hvtn702.lowProp" = .hvtn702.df$lowProp );
    pdf( "rv144lowProp_hvtn702lowProp.pdf" ); ggplot( .lowProp.df, aes( x=rv144.lowProp,y=hvtn702.lowProp ) ) + geom_point() + stat_density2d_filled() + ggtitle( "RV144 vs HVTN702 % of prop at low risk (lowProp)" ); dev.off()
    
} # END IF FALSE

