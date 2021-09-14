library(deSolve)
library(tidyverse)
library(EpiModel)
library(survival)
library(EasyABC)
library(ggplot2)
library(pdfCluster)

si.ode.twogroup.fn <- function ( times, init, param ) {
  with(as.list(c(init, param)), {

    # Flows
    # the number of people moving from the S to I compartment at each time step

    # The one trial

    #PLACEBO arm
    #Susceptible, Infected, placebo, high, medium, low
    #SIph.flow <- risk*lambda*Sph # line from original model, FYI
    SIph.flow <- unname( high.risk.multiplier*(10^log10lambda)*Sph ) # Because this can't go negative, high.risk.multiplier must be <= 1/(10^log10lambda). See bounds defs.
    SIpl.flow <- unname( (10^log10lambda)*Spl )

    #VACCINE arm
    #Susceptible, Infected, vaccine, high, medium, low
    SIvh.flow <- unname( high.risk.multiplier*(10^log10lambda)*(1-epsilon)*Svh )
    SIvl.flow <- unname( (10^log10lambda)*(1-epsilon)*Svl )

    # ODEs
    # placebo; heterogeneous high.risk.multiplier
    # original ODE:  dSph <- -SIph.flow
    dSph <- -SIph.flow
    dIph <- SIph.flow  #high.risk.multiplier*lambda*Sph
    dSpl <- -SIpl.flow
    dIpl <- SIpl.flow  #lambda*Spl

    # vaccine; heterogeneous high.risk.multiplier
    dSvh <- -SIvh.flow
    dIvh <- SIvh.flow  #high.risk.multiplier*lambda*(1-epsilon)*Svh
    dSvl <- -SIvl.flow
    dIvl <- SIvl.flow  #lambda*Svl


    #Output
    list(c(
           dSph,dIph,
           dSpl,dIpl,
           dSvh,dIvh,
           dSvl,dIvl,
           SIph.flow,SIpl.flow,
           SIvh.flow,SIvl.flow
           ))
  })
} # si.ode.twogroup.fn (..)
mod.manipulate.twogroup <- function( mod ) {
  #browser()

  # TWOGROUP
  mod <- mutate_epi(mod, total.Svh.Svl = Svh + Svl) #all susceptible in heterogeneous risk vaccine pop
  mod <- mutate_epi(mod, total.Sph.Spl = Sph + Spl) #all susceptible in heterogeneous risk placebo pop
  mod <- mutate_epi(mod, total.Ivh.Ivl = Ivh + Ivl) #all infected in heterogeneous risk vaccine pop
  mod <- mutate_epi(mod, total.Iph.Ipl = Iph + Ipl) #all infected in heterogeneous risk placebo pop
  mod <- mutate_epi(mod, total.SIvh.SIvl.flow = SIvh.flow + SIvl.flow) #all infections per day in heterogeneous risk vaccine pop
  mod <- mutate_epi(mod, total.SIph.SIpl.flow = SIph.flow + SIpl.flow) #all infections in heterogeneous risk placebo pop

  #Incidence estimates, per 100 person years


  #Instantaneous incidence / hazard
  mod <- mutate_epi(mod, rate.Vaccine.het = (total.SIvh.SIvl.flow/total.Svh.Svl)*365*100)
  mod <- mutate_epi(mod, rate.Placebo.het = (total.SIph.SIpl.flow/total.Sph.Spl)*365*100)

  #Cumulative incidence
  mod <- mutate_epi(mod, cumul.Svh.Svl = cumsum(total.Svh.Svl))
  mod <- mutate_epi(mod, cumul.Sph.Spl = cumsum(total.Sph.Spl))
  mod <- mutate_epi(mod, cumul.rate.Vaccine.het = (total.Ivh.Ivl/cumul.Svh.Svl)*365*100)
  mod <- mutate_epi(mod, cumul.rate.Placebo.het = (total.Iph.Ipl/cumul.Sph.Spl)*365*100)

  #Vaccine efficacy (VE) estimates
  #VE <- 1 - Relative Risk; this is VE for instantaneous incidence / hazard
  mod <- mutate_epi(mod, VE.inst = 100 * ( 1 - rate.Vaccine.het/rate.Placebo.het ) )

  #VE <- 1 - Relative Risk; this is VE from cumulative incidence
  mod <- mutate_epi(mod, VE.cumul = 100 * ( 1 - cumul.rate.Vaccine.het/cumul.rate.Placebo.het ) )

  return( mod );
} # mod.manipulate.twogroup (..)

run.and.compute.run.stats.twogroup <- function (
      epsilon,   #per contact vaccine efficacy
      log10lambda,     #log10( beta*c*prev ),
      high.risk.multiplier,          # Risk multiplier for high risk group
      highRiskProportion,
      vaccinatedProportion = 0.5,  # In lieu of naming vaccine and placebo arms separately (and their N)
      trialSize = 10000,  # Now we just have to add this magic number for size
      trial.evaluation.time = 3*365
) {
      param <- param.dcm(epsilon = epsilon, log10lambda = log10lambda, high.risk.multiplier = high.risk.multiplier );
 
      # initial values
      Svh <- floor( highRiskProportion * vaccinatedProportion * trialSize );
      Sph <- floor( highRiskProportion * ( 1.0 - vaccinatedProportion ) * trialSize );
      Svl <- floor( ( 1.0 - highRiskProportion ) * vaccinatedProportion * trialSize );
      Spl <- floor( ( 1.0 - highRiskProportion ) * ( 1.0 - vaccinatedProportion ) * trialSize );

      Sp <- Spl + Sph;
      Sv <- Svl + Svh;
      if( vaccinatedProportion == 0.5 ) {
          stopifnot( Sp == Sv );
      }

      init <- init.dcm(Sph = Sph, Iph = 0,    #placebo, high risk
                       Spl = Spl, Ipl = 0,    #placebo, low risk
                       Svh = Svh, Ivh = 0,    #vaccine
                       Svl = Svl, Ivl = 0,    #vaccine
                       SIph.flow = 0, SIpl.flow = 0,
                       SIvh.flow = 0, SIvl.flow = 0
                       );

      control <- control.dcm( nsteps = trial.evaluation.time, new.mod = si.ode.twogroup.fn );
      mod <- dcm( param, init, control );
      #print( mod )
      mod.with.stats <- mod.manipulate.twogroup( mod );
      #print( mod.with.stats )
      mod.with.stats.df <- as.data.frame( mod.with.stats );

      # heterogeneous risk using cumulative VE:
      VE <- mod.with.stats.df$VE.cumul[ trial.evaluation.time ];

      ## The placebo incidence out stat vector is the _cumulative_ incidence at time trial.evaluation.time.
      placeboIncidence <- mod.with.stats.df$cumul.rate.Placebo.het[ trial.evaluation.time ];

      c( VE = VE, placeboIncidence = placeboIncidence );
} # run.and.compute.run.stats.twogroup (..)

common.parameters <- c( "epsilon" );
trial.parameters <- c( "log10lambda", "high.risk.multiplier", "highRiskProportion" );
