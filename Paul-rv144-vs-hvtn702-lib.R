library(deSolve)
library(tidyverse)
library(EpiModel)
library(survival)
library(EasyABC)
library(ggplot2)
library(pdfCluster)

si.ode.rv144.hvtn702.fn <- function ( times, init, param ) {
  with(as.list(c(init, param)), {
    
    # Flows

    # rv144
    # the number of people moving from S to I at each time step
    #Susceptible, Infected, placebo, high, medium, low
    rv144.SIph.flow <- rv144.high.risk.multiplier*rv144.lambda*rv144.Sph
    rv144.SIpm.flow <- rv144.lambda*rv144.Spl
    rv144.SIpl.flow <- 0*rv144.lambda*rv144.Spl  #0 to give this group zero exposures
    
    #Susceptible, Infected, vaccine, high, medium, low
    rv144.SIvh.flow <- rv144.high.risk.multiplier*rv144.lambda*(1-epsilon)*rv144.Svh
    ### Paul found this line has a bug:
    ### rv144.SIvm.flow <- rv144.lambda*(1-epsilon)*rv144.Spl
    rv144.SIvm.flow <- rv144.lambda*(1-epsilon)*rv144.Svl
    rv144.SIvl.flow <- 0*rv144.lambda*(1-epsilon)*rv144.Svl  #0 to give this group zero exposures
    
    # ODEs
    # placebo; heterogeneous rv144.high.risk.multiplier
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
    


    # hvtn702
    # the number of people moving from S to I at each time step
    #Susceptible, Infected, placebo, high, medium, low
    hvtn702.SIph.flow <- hvtn702.high.risk.multiplier*hvtn702.lambda*hvtn702.Sph
    hvtn702.SIpm.flow <- hvtn702.lambda*hvtn702.Spl
    hvtn702.SIpl.flow <- 0*hvtn702.lambda*hvtn702.Spl  #0 to give this group zero exposures
    
    #Susceptible, Infected, vaccine, high, medium, low
    hvtn702.SIvh.flow <- hvtn702.high.risk.multiplier*hvtn702.lambda*(1-epsilon)*hvtn702.Svh
    ### Paul found this line has a bug:
    ### hvtn702.SIvm.flow <- hvtn702.lambda*(1-epsilon)*hvtn702.Spl
    hvtn702.SIvm.flow <- hvtn702.lambda*(1-epsilon)*hvtn702.Svl
    hvtn702.SIvl.flow <- 0*hvtn702.lambda*(1-epsilon)*hvtn702.Svl  #0 to give this group zero exposures
    
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
} # si.ode.rv144.hvtn702.fn (..)


mod.manipulate <- function( mod ) {
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
  mod <- mutate_epi(mod, rv144.VE.inst = 1 - rv144.rate.Vaccine.het/rv144.rate.Placebo.het)
  
  #VE <- 1 - Relative Risk; this is VE from cumulative incidence
  mod <- mutate_epi(mod, rv144.VE.cumul = 1 - cumul.rv144.rate.Vaccine.het/cumul.rv144.rate.Placebo.het)


  
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
  mod <- mutate_epi(mod, hvtn702.VE.inst = 1 - hvtn702.rate.Vaccine.het/hvtn702.rate.Placebo.het)
  
  #VE <- 1 - Relative Risk; this is VE from cumulative incidence
  mod <- mutate_epi(mod, hvtn702.VE.cumul = 1 - cumul.hvtn702.rate.Vaccine.het/cumul.hvtn702.rate.Placebo.het)

  #return(mod)
} # mod.manipulate (..)

## This computes the distance as calculated internally in the abc function - but not there the distances are normalized using "normalise" which we think centralizes too (so makes each scaled stat have mean 0, sd 1) - the standard deviations are saved and included in the abc output and need to be passed into here, because after keeping only the closest X% of the samples, the remaining samples will have a different STDEV. So to get the right distances it's important to use the correct stdev. These values are not centralized before the distances are computed but this should not matter.
# example: calculate.abc.dist( fit.rej$stats, as.numeric( target.stats[-1] ), fit.rej$stats_normalization )
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

#------------------------------------------------------------------------------
# sim execution
#------------------------------------------------------------------------------
runSim_rv144.hvtn702 <- function(reac = c( "numExecution" = 1000 ) ) {
    stopifnot( all( c( "numExecution" ) %in% names( reac ) ) );

    ## Number of parameters to optimize (3, 4, or 5).
    num.params <- 9; # This is now fixed at 9.

    num.sims <- unname( reac[ "numExecution" ] );
    abc.keep.num.sims <- 250; # Number of samples to run through the clustering and optimizing steps.
    stopifnot( num.sims > abc.keep.num.sims ); # It won't work to cluster uniformly drawn points. You first have to filter them by keeping those nearest the target.

    rv144.VE <- 0.31;
    rv144.placeboIncidence <- 0.14;
    hvtn702.VE <- 0;
    hvtn702.placeboIncidence <- 3.3;

    rv144.placebo.incidence.target <- rv144.placeboIncidence; # incidence per 100 person years, eg 3.5 for 702 or 0.3 for RV144 [todo: double-check these values!]
    rv144.VE.target = rv144.VE; # instantaneous VE observed by the end of the trial
    hvtn702.placebo.incidence.target <- hvtn702.placeboIncidence; # incidence per 100 person years, eg 3.5 for 702 or 0.3 for HVTN702 [todo: double-check these values!]
    hvtn702.VE.target = hvtn702.VE; # instantaneous VE observed by the end of the trial

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
    placebo.incidence.target.scale.units <- 0.1; # These values mean that a 0.1 difference in placebo incidence from the target placebo incidence should cost about the same as a 0.01 difference in VE from its target. VE is meansured on a scale of -1 to 1, with units typically refered to as percentages, eg rv144 had a 31% efficacy is 0.31, so this means a "1% VE difference" has about the same cost in our distance metrix as does a 0.1 per 100-person year difference (so that's the difference between 3.3 as seen in 702 vs 3.2 or 3.3 per 100 person-years -- note that for rv144 the estimated placebo incidence was about 0.14 per 100 person-years, see below for references).
    VE.target.scale.units <- 0.01;

    target.stat.scale.units <- c( "rv144.VE.target" = VE.target.scale.units, "rv144.placebo.incidence.target" = placebo.incidence.target.scale.units, "hvtn702.VE.target" = VE.target.scale.units, "hvtn702.placebo.incidence.target" = placebo.incidence.target.scale.units );

    #browser()

    ## Note that we target the average incidence over time, rather than the incidence at the end of the trial, since with multiple risk groups there is waning expected.

    # MAGIC NUMBERS
    trial.size <- 10000; # Number of participants (vaccine & placebo recipients)
    trial.frac.vaccinated <- 0.5; # 1:1 vaccine:placebo treatment assigment ratio
    trial.evaluation.time <- 3*365; # End of the trial.
    default.high.risk.proportion <- 0.1; # For 3 parameter version.
    default.low.risk.proportion.among.those.not.high.risk <- 0.1; # For 3 and 4 parameter venrsions.
    high.risk.multiplier.max <- 50;
    lambda.min <- 1E-7;
    lambda.max <- 1E-2;
    smallest.discernable.amount <- 5E-4; # determined by trial and error this is the smallest amount you can change the parameters from 0 or 1 for it to register a difference from those extremes, eg to avoid NaN and Inf

    target.stats <- data.frame( trial.evaluation.time, rv144.VE.target, rv144.placebo.incidence.target, hvtn702.VE.target, hvtn702.placebo.incidence.target );
    nsteps <- trial.evaluation.time;

    run.and.compute.run.stats <- function (
      epsilon,   #per contact vaccine efficacy
      rv144.lambda,     #beta*c*prev,
      hvtn702.lambda,     #beta*c*prev,
      rv144.high.risk.multiplier,          #risk multiplier for high risk group
      hvtn702.high.risk.multiplier,
      rv144.highRiskProportion = default.high.risk.proportion, # these have defaults, so you can run this as 3 or 4 or 5 parameters.
      hvtn702.highRiskProportion = default.high.risk.proportion,
      rv144.lowRiskProportion = ( default.low.risk.proportion.among.those.not.high.risk/(1-rv144.highRiskProportion) ), # Among those not high risk.
      hvtn702.lowRiskProportion = ( default.low.risk.proportion.among.those.not.high.risk/(1-hvtn702.highRiskProportion) ),
      vaccinatedProportion = trial.frac.vaccinated,
      trialSize = trial.size
    ) {
      
      # Paul added the other params to this (risk here, others below):
      param <- param.dcm(epsilon = epsilon, rv144.lambda = rv144.lambda, hvtn702.lambda = hvtn702.lambda, rv144.high.risk.multiplier = rv144.high.risk.multiplier, hvtn702.high.risk.multiplier = hvtn702.high.risk.multiplier )

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

      init <- init.dcm(rv144.Sph = rv144.Sph, rv144.Iph = 0,    #placebo, high risk
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
      
      control <- control.dcm( nsteps = nsteps, new.mod = si.ode.rv144.hvtn702.fn );
      mod <- dcm( param, init, control );
      #print( mod )
      
      mod.with.stats <- mod.manipulate( mod );
      #print( mod.with.stats )
      mod.with.stats.df <- as.data.frame( mod.with.stats );
      
      # heterogeneous risk using cumulative VE:
      rv144.VE <- mod.with.stats.df$rv144.VE.cumul[target.stats$trial.evaluation.time];
      hvtn702.VE <- mod.with.stats.df$hvtn702.VE.cumul[target.stats$trial.evaluation.time];
      
      ## The placebo incidence out stat vector is the mean over time up to each time in target.stats$trial.evaluation.time.
      rv144.meanPlaceboIncidence <-
        sapply( target.stats$trial.evaluation.time, function( .time ) { mean( mod.with.stats.df$rv144.rate.Placebo.het[ 1:.time ] ) } );
      hvtn702.meanPlaceboIncidence <-
        sapply( target.stats$trial.evaluation.time, function( .time ) { mean( mod.with.stats.df$hvtn702.rate.Placebo.het[ 1:.time ] ) } );

      c( rv144.VE = rv144.VE, rv144.meanPlaceboIncidence = rv144.meanPlaceboIncidence, hvtn702.VE = hvtn702.VE, hvtn702.meanPlaceboIncidence = hvtn702.meanPlaceboIncidence );
    } # run.and.compute.run.stats (..)

    # Specify bounds for the initial conditions; these are used as priors.
    # In order of "x" in the above function.
    bounds <- list(
                   epsilon = c(smallest.discernable.amount, 1-smallest.discernable.amount), # This is just the full range 0 to 1
                   rv144.lambda = c(lambda.min, lambda.max),
                   rv144.high.risk.multiplier = c(1, high.risk.multiplier.max), # risk multiplier for high risk group
                   rv144.highRiskProportion = c(smallest.discernable.amount, 1-smallest.discernable.amount), # This is just the full range 0 to 1
                   rv144.lowRiskProportion = c(smallest.discernable.amount, 1-smallest.discernable.amount), # This is modified to avoid NaN/Inf
                   hvtn702.lambda = c(lambda.min, lambda.max),
                   hvtn702.high.risk.multiplier = c(1, high.risk.multiplier.max), # risk multiplier for high risk group
                   hvtn702.highRiskProportion = c(smallest.discernable.amount, 1-smallest.discernable.amount), # This is just the full range 0 to 1
                   hvtn702.lowRiskProportion = c(smallest.discernable.amount, 1-smallest.discernable.amount) # This is modified to avoid NaN/Inf
                   );            

    # In order of "x" in the above function.
    priors  <- lapply( bounds, function( .bounds ) { c( "unif", unlist( .bounds ) ) } );

    make.abc.fn <- function ( target ) {
        function( x ) {
            run.and.compute.run.stats( epsilon = x[ 1 ], rv144.lambda = x[ 2 ], rv144.high.risk.multiplier = x[ 3 ], rv144.highRiskProportion = x[ 4 ], rv144.lowRiskProportion = x[ 5 ], hvtn702.lambda = x[ 6 ], hvtn702.high.risk.multiplier = x[ 7 ], hvtn702.highRiskProportion = x[ 8 ], hvtn702.lowRiskProportion = x[ 9 ] )
        }
    }

    .f.abc <- make.abc.fn( target.stats[-1] )
    fit.rej <- ABC_rejection(model = .f.abc,
                         prior = priors,
                         nb_simul = num.sims,
                         summary_stat_target = as.numeric( target.stats[-1] ),
                         tol = 1.0, # Keep all of them for the moment; see below.
                         progress_bar = TRUE)

    # Compute the distances for each retained sample
    fit.rej.dist <-
        calculate.abc.dist( fit.rej$stats, as.numeric( target.stats[-1] ), ifelse( is.na( fit.rej$stats_normalization ), 1, fit.rej$stats_normalization ) );

    # Filter to keep fewer points.
    # Divide distance into units the width of the top abc.keep.num.sims values.
    dist.units <- quantile( fit.rej.dist, probs = ( abc.keep.num.sims / num.sims ) )
    fit.rej.dist.scaled <- fit.rej.dist / dist.units;
    fit.rej.dist.scaled.int <- floor( fit.rej.dist.scaled );

    ## Make a contour-like plot
    ## ERE I AM
    # .df <- as.data.frame( fit.rej$param ) #df of just the parameter combinations ABC sampled
    # names( .df ) <- c( "epsilon", "lambda", "risk" );
    # .df <- cbind( .df, fit.rej.dist, fit.rej.dist.scaled );
    # .df <- .df[ fit.rej.dist.scaled.int <= 2, , drop = FALSE ];
    # ggplot( .df, aes( x=lambda,y=risk, alpha = 1-(fit.rej.dist.scaled*(1/2)) ) ) + geom_point()

    # Filter to keep points up to max.fit.rej.dist.scaled units away.
    max.fit.rej.dist.scaled <- 1.5; # MAGIC #
    fit.orig <- fit.rej;
    abc.keep.sim <- fit.rej.dist.scaled < max.fit.rej.dist.scaled;
    fit.rej$param <- fit.rej$param[ abc.keep.sim, , drop = FALSE ];
    fit.rej$stats <- fit.rej$stats[ abc.keep.sim, , drop = FALSE ];
    fit.rej$weights <- fit.rej$weights[ abc.keep.sim ];
    fit.rej.dist <- fit.rej.dist[ abc.keep.sim ];
    # Now find optimal points near the sampled modes.
    
    # Globally best sample is at this index:
    # sample.index.minimizing.dist <- which.min( fit.rej.dist );
    # This might fail with the error that there's too course of a grid.
    # from help( "pdfCluster" ):
    #     # Warning:
    # 
    #      It may happen that the variability of the estimated density is so
    #      high that not all jumps in the mode function can be detected by
    #      the selected grid scanning the density function. In that case, no
    #      output is produced and a message is displayed. As this may be
    #      associated to the occurrence of some spurious connected
    #      components, which appear and disappear within the range between
    #      two subsequent values of the grid, a natural solution is to
    #      increase the value of ‘n.grid’.  Alternatively either ‘lambda’ or
    #      ‘hmult’ may be increased to alleviate the chance of detecting
    #      spurious connected components.
    #
    # I TRIED first and it worked for 3 params just adding: bwtype="adaptive" 
    # For 4 params that was not sufficient so I TRIED also adding (and it worked, despite actually just setting it to n): n.grid=1E5
    # For 4 params that was too slow so I TRIED instead adding hmult = 1.25 and that worked but had too few modes, since I knew there were others. So I played around and hmult = 1 which I thought was the default gives 5 modes which feels ok.
    # For 5 params that was also sufficient.
    # if( ncol( fit.rej$param ) < 5 ) {
    #     cl <- suppressWarnings( pdfCluster( fit.rej$param, bwtype="adaptive", hmult=1 ) );
    # } else {
    #     cl <- pdfCluster( fit.rej$param, bwtype="adaptive", hmult=1, n.grid=nrow(fit.rej$param) );
    # }
    # Update: with current params (keeping 250 draws only, and with wide enough parameter ranges), see abc.keep.num.sims, this seems to work fine and better, even, in that the modes are good:
    # cl <- suppressWarnings( pdfCluster( fit.rej$param ) );
    # Further update: with new 9-param model, we need this:
    cl <- pdfCluster( fit.rej$param, bwtype="adaptive", hmult=1, n.grid=nrow(fit.rej$param) );

    cluster.numbers <- groups( cl );
    # boxplot( fit.rej.dist ~ cluster.numbers )
    
    # Cluster-specific best sample is at this index, for each cluster:
    cluster.member.minimizing.dist <- sapply( 1:max( cluster.numbers ), function( .cluster.number ) { cluster.indices <- which( cluster.numbers == .cluster.number ); return( cluster.indices[ which.min( fit.rej.dist[ cluster.indices ] ) ] ); } );
    
    make.optim.fn <- function ( .target.stats, .target.stat.stdevs = rep( 1, length( .target.stats ) ) ) {
        function( x ) {
            .stats <- .f.abc( x );
            .stats.matrix <- matrix( .stats, nrow = 1 );
            c( dist = calculate.abc.dist( .stats.matrix, .target.stats, ifelse( is.na( .target.stat.stdevs ), 1, .target.stat.stdevs ) ) )
        }
    } # make.optim.fn (..)

    # OLD: .f <- make.optim.fn( target.stats[-1], fit.rej$stats_normalization );
    # NEW:
    .f <- make.optim.fn( target.stats[-1], target.stat.scale.units ); # For better balancing of the costs; see above.

    # These are the original bounds, separated for use in optim:
    lower.bounds <- sapply( bounds, function ( .bounds.for.x ) { .bounds.for.x[ 1 ] } );
    upper.bounds <- sapply( bounds, function ( .bounds.for.x ) { .bounds.for.x[ 2 ] } );

    optima.by.cluster <- sapply( cluster.member.minimizing.dist, function( .minimizer.index ) {
        print( .minimizer.index );
        print( fit.rej$param[ .minimizer.index, ] );
        print( .f( fit.rej$param[ .minimizer.index, ] ) );
        .optim.result <- optim( fit.rej$param[ .minimizer.index, ], .f, lower = lower.bounds, upper = upper.bounds, method = "L-BFGS-B" );
        .par <- .optim.result$par;
        names( .par ) <- names( bounds );
        .stats <- .f.abc( .par );
        print( c( .par, .stats, dist = .optim.result$value ) );
        return( c( .par, .stats, dist = .optim.result$value ) );
    } );
    optima.by.cluster.sorted <- optima.by.cluster[ , order( optima.by.cluster[ "dist", ] ), drop = FALSE ];

    return( list( fit = fit.rej, priors = priors, bounds = bounds, target.stats = target.stats, fn = .f.abc, sampled.modes = optima.by.cluster.sorted ) );
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
    
    # These are the original bounds, separated for use in optim:
    lower.bounds <- sapply( bounds, function ( .bounds.for.x ) { .bounds.for.x[ 1 ] } );
    upper.bounds <- sapply( bounds, function ( .bounds.for.x ) { .bounds.for.x[ 2 ] } );

    plot(density(fit$param[, 1], from = lower.bounds[1],  to = upper.bounds[1]),
         main = "epsilon", 
         xlim = c(lower.bounds[1], upper.bounds[1]),
         #ylim = c(0, 10),
         col=2)
    lines(density(fit$param[, 1], from = lower.bounds[1],  to = upper.bounds[1]), col = 2)
    abline(v = VE.target, lty = 2, col = 1)
    legend("topright", legend = c("per-contact VE", "Posterior"),
           lty = c(1, 2), col = 1:2, lwd = 2)
    
    plot(density(fit$param[, 2], from = lower.bounds[2],  to = upper.bounds[2]),
         main = "lambda", 
         xlim = c(lower.bounds[2], upper.bounds[2]),
         col=2)
    lines(density(fit$param[, 2], from = lower.bounds[2],  to = upper.bounds[2]), col = 2)
    #abline(v = VE.target, lty = 2, col = 1) # This is a bug, it plots VE target, not lambda
    legend("topright", legend = c("Truth", "Posterior"),
           lty = c(1, 2), col = 1:2, lwd = 2)
    
    plot(density(fit$param[, 3], from = lower.bounds[3],  to = upper.bounds[3]),
         main = "risk", 
         xlim = c(lower.bounds[3], upper.bounds[3]),
         col=2)
    lines(density(fit$param[, 3], from = lower.bounds[3],  to = upper.bounds[3]), col = 2)
    #abline(v = VE.target, lty = 2, col = 1) # Another bug
    legend("topright", legend = c("Truth", "Posterior"),
           lty = c(1, 2), col = 1:2, lwd = 2)
    
    .df <- as.data.frame( fit$param ) #df of just the parameter combinations ABC sampled
    names( .df ) <- c( "epsilon", "lambda", "risk" )
    pdf( "lambda_risk.pdf" ); ggplot( .df, aes( x=lambda,y=risk ) ) + geom_point() + stat_density2d_filled(); dev.off()
    pdf( "epsilon_risk.pdf" ); ggplot( .df, aes( x=epsilon,y=risk ) ) + geom_point() + stat_density2d_filled(); dev.off()
    pdf( "epsilon_lambda.pdf" ); ggplot( .df, aes( x=epsilon,y=lambda ) ) + geom_point() + stat_density2d_filled(); dev.off()
} # END IF FALSE

