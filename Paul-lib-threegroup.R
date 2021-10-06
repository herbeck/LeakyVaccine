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
    SIph.flow <- unname( (10^log10riskmultiplier)*(10^log10lambda)*Sph )  # Because this can't go negative, (10^log10riskmultiplier) must be <= 1/(10^log10lambda). See bounds defs.
    SIpm.flow <- unname( (10^log10lambda)*Spm )
    SIpl.flow <- unname( 0*(10^log10lambda)*Spl )  #0 to give this group zero exposures
       # Could also use "1/(10^log10riskmultiplier)" if we don't want ZERO exposures
    
    #VACCINE arm
    #Susceptible, Infected, vaccine, high, medium, low
    SIvh.flow <- unname( (10^log10riskmultiplier)*(10^log10lambda)*(1-epsilon)*Svh )
    SIvm.flow <- unname( (10^log10lambda)*(1-epsilon)*Svm )
    SIvl.flow <- unname( 0*(10^log10lambda)*(1-epsilon)*Svl )  #0 to give this group zero exposures
       # Could also use "1/(10^log10riskmultiplier)" if we don't want ZERO exposures
    
    # ODEs
    # placebo; heterogeneous (10^log10riskmultiplier)
    # original ODE:  dSph <- -SIph.flow
    dSph <- -SIph.flow
    dIph <- SIph.flow  #(10^log10riskmultiplier)*lambda*Sph
    dSpm <- -SIpm.flow
    dIpm <- SIpm.flow  #lambda*Spm
    dSpl <- -SIpl.flow
    dIpl <- SIpl.flow  #0*lambda*Spl
    
    # vaccine; heterogeneous (10^log10riskmultiplier)
    dSvh <- -SIvh.flow
    dIvh <- SIvh.flow  #(10^log10riskmultiplier)*lambda*(1-epsilon)*Svh
    dSvm <- -SIvm.flow
    dIvm <- SIvm.flow  #lambda*Svm
    dSvl <- -SIvl.flow
    dIvl <- SIvl.flow  #0*lambda*(1-epsilon)*Svl


    #Output
    list(c(0,0, # hack/test
           dSph,dIph,
           dSpm,dIpm,
           dSpl,dIpl,
           dSvh,dIvh,
           dSvm,dIvm,
           dSvl,dIvl,
           0, # hack/test
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
      log10riskmultiplier,          # log10( risk multiplier for high risk group )
      highRiskProportion,
      lowRiskProportion,             # This is a proportion among those not high risk
      vaccinatedProportion = 0.5,  # In lieu of naming vaccine and placebo arms separately (and their N)
      VE.start = 0, # must be scalar
      VE.end = 0, # can be a vector
      trialSize = 10000,  # Now we just have to add this magic number for size
      sim.ndays = base::max( VE.end )
) {
    stopifnot( length( VE.start ) == 1 );

    param <- param.dcm( epsilon = epsilon, log10lambda = log10lambda, log10riskmultiplier = log10riskmultiplier );

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
          #stopifnot( Sp == Sv );
          if( Sp != Sv ) {
              cat( paste( "NOTE: Sp=",Sp, " is not equal to Sv=", Sv, sep = "" ), fill = TRUE );
          }
      }

      init <- init.dcm(Sxx = 0, Ixx = 0, #hack/test
                       Sph = Sph, Iph = 0,    #placebo, high risk
                       Spm = Spm, Ipm = 0,   #placebo, medium risk
                       Spl = Spl, Ipl = 0,    #placebo, low risk
                       Svh = Svh, Ivh = 0,    #vaccine
                       Svm = Svm, Ivm = 0,   #vaccine
                       Svl = Svl, Ivl = 0,    #vaccine
                       SIxx.flow = 0, # hack/test
                       SIph.flow = 0, SIpm.flow = 0, SIpl.flow = 0,
                       SIvh.flow = 0, SIvm.flow = 0, SIvl.flow = 0
                       );
      
      control <- control.dcm( nsteps = sim.ndays, new.mod = si.ode.threegroup.fn );
      mod <- dcm( param, init, control );
      #print( mod )
      
      mod.with.stats <- mod.manipulate.threegroup( mod );
      #print( mod.with.stats )
      mod.with.stats.df <- as.data.frame( mod.with.stats );
      
      # heterogeneous risk using cumulative VE:
      if( VE.start <= 1 ) { # 0 or 1 because the start and end are inclusive.
          VE <- mod.with.stats.df$VE.cumul[ VE.end ];
      } else {
          cum.since.start.rate.Vaccine.het <-
              ( # VE.start - 1 because the start and end are inclusive.
               ( mod.with.stats.df$total.Ivh.Ivm.Ivl[ VE.end ] - mod.with.stats.df$total.Ivh.Ivm.Ivl[ max( 1, VE.start - 1 ) ] ) /
               ( mod.with.stats.df$cumul.Svh.Svm.Svl[ VE.end ] - mod.with.stats.df$cumul.Svh.Svm.Svl[ max( 1, VE.start - 1 ) ] )
              );
          cum.since.start.rate.Placebo.het <-
              (
               ( mod.with.stats.df$total.Iph.Ipm.Ipl[ VE.end ] - mod.with.stats.df$total.Iph.Ipm.Ipl[ max( 1, VE.start - 1 ) ] ) /
               ( mod.with.stats.df$cumul.Sph.Spm.Spl[ VE.end ] - mod.with.stats.df$cumul.Sph.Spm.Spl[ max( 1, VE.start - 1 ) ] )
              );
          VE <- 100 *
              ( 1 -
               ( cum.since.start.rate.Vaccine.het / cum.since.start.rate.Placebo.het )
              );
      }
      ## The placebo incidence out stat vector is the _cumulative_ incidence at time sim.ndays.
      placeboIncidence <- mod.with.stats.df$cumul.rate.Placebo.het[ VE.end ];

    .rv <- data.frame( VE = VE, placeboIncidence = placeboIncidence );
    rownames( .rv ) <- VE.end;
    return( .rv );
} # run.and.compute.run.stats.threegroup (..)

common.parameters <- c( "epsilon" );
trial.parameters <- c( "log10lambda", "log10riskmultiplier", "highRiskProportion", "lowRiskProportion" );
