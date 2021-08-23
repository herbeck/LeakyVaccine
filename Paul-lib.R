library(deSolve)
library(tidyverse)
library(EpiModel)
library(survival)
library(EasyABC)
library(ggplot2)

si_odePaul <- function(times, init, param){
  with(as.list(c(init, param)), {
    
    # Flows
    # the number of people moving from S to I at each time step
    #Susceptible, Infected, placebo
    SIp.flow <- lambda*Sp
    SIv.flow <- lambda*(1-epsilon)*Sv
    
    #Susceptible, Infected, placebo, high, medium, low
    SIph.flow <- risk*lambda*Sph
    SIpm.flow <- lambda*Spl
    SIpl.flow <- 0*lambda*Spl  #0 to give this group zero exposures
    
    #Susceptible, Infected, vaccine, high, medium, low
    SIvh.flow <- risk*lambda*(1-epsilon)*Svh
    ### Paul found this line has a bug:
    ### SIvm.flow <- lambda*(1-epsilon)*Spl
    SIvm.flow <- lambda*(1-epsilon)*Svl
    SIvl.flow <- 0*lambda*(1-epsilon)*Svl  #0 to give this group zero exposures
    
    # ODEs
    # placebo; homogeneous risk
    dSp <- -SIp.flow
    dIp <- SIp.flow  #lambda*Sp
    
    # vaccine; homogeneous risk
    dSv <- -SIv.flow
    dIv <- SIv.flow  #lambda*epsilon*Sv
    
    # placebo; heterogeneous risk
    dSph <- -SIph.flow
    dIph <- SIph.flow  #risk*lambda*Sph
    dSpm <- -SIpm.flow
    dIpm <- SIpm.flow  #lambda*Spm
    dSpl <- -SIpl.flow
    dIpl <- SIpl.flow  #0*lambda*Spl
    
    # vaccine; heterogeneous risk
    dSvh <- -SIvh.flow
    dIvh <- SIvh.flow  #risk*lambda*(1-epsilon)*Svh
    dSvm <- -SIvm.flow
    dIvm <- SIvm.flow  #lambda*Svm
    dSvl <- -SIvl.flow
    dIvl <- SIvl.flow  #0*lambda*(1-epsilon)*Svl
    
    #Output
    list(c(dSp,dIp,
           dSv,dIv,
           dSph,dIph,
           dSpm,dIpm,
           dSpl,dIpl,
           dSvh,dIvh,
           dSvm,dIvm,
           dSvl,dIvl,
           SIp.flow,SIv.flow,
           SIph.flow,SIpm.flow,SIpl.flow,
           SIvh.flow,SIvm.flow,SIvl.flow))
  })
}


mod.manipulate <- function(mod){
  #browser()
  mod <- mutate_epi(mod, total.Svh.Svm.Svl = Svh + Svm + Svl) #all susceptible in heterogeneous risk vaccine pop
  mod <- mutate_epi(mod, total.Sph.Spm.Spl = Sph + Spm + Spl) #all susceptible in heterogeneous risk placebo pop
  mod <- mutate_epi(mod, total.Ivh.Ivm.Ivl = Ivh + Ivm + Ivl) #all infected in heterogeneous risk vaccine pop
  mod <- mutate_epi(mod, total.Iph.Ipm.Ipl = Iph + Ipm + Ipl) #all infected in heterogeneous risk placebo pop
  mod <- mutate_epi(mod, total.SIvh.SIvm.SIvl.flow = SIvh.flow + SIvm.flow + SIvl.flow) #all infections per day in heterogeneous risk vaccine pop
  mod <- mutate_epi(mod, total.SIph.SIpm.SIpl.flow = SIph.flow + SIpm.flow + SIpl.flow) #all infections in heterogeneous risk placebo pop
  
  #Incidence estimates, per 100 person years
  #Instantaneous incidence / hazard
  mod <- mutate_epi(mod, rate.Vaccine = (SIv.flow/Sv)*365*100)
  mod <- mutate_epi(mod, rate.Placebo = (SIp.flow/Sp)*365*100)
  mod <- mutate_epi(mod, rate.Vaccine.het = (total.SIvh.SIvm.SIvl.flow/total.Svh.Svm.Svl)*365*100)
  mod <- mutate_epi(mod, rate.Placebo.het = (total.SIph.SIpm.SIpl.flow/total.Sph.Spm.Spl)*365*100)
  
  #Cumulative incidence
  mod <- mutate_epi(mod, cumul.Sv = cumsum(Sv))
  mod <- mutate_epi(mod, cumul.Sp = cumsum(Sp))
  mod <- mutate_epi(mod, cumul.Svh.Svm.Svl = cumsum(total.Svh.Svm.Svl))
  mod <- mutate_epi(mod, cumul.Sph.Spm.Spl = cumsum(total.Sph.Spm.Spl))
  mod <- mutate_epi(mod, cumul.rate.Vaccine = (Iv/cumul.Sv)*365*100)
  mod <- mutate_epi(mod, cumul.rate.Placebo = (Ip/cumul.Sp)*365*100)
  mod <- mutate_epi(mod, cumul.rate.Vaccine.het = (total.Ivh.Ivm.Ivl/cumul.Svh.Svm.Svl)*365*100)
  mod <- mutate_epi(mod, cumul.rate.Placebo.het = (total.Iph.Ipm.Ipl/cumul.Sph.Spm.Spl)*365*100)
  
  #Vaccine efficacy (VE) estimates
  #VE <- 1 - Relative Risk; this is VE for instantaneous incidence / hazard
  mod <- mutate_epi(mod, VE1.inst = 1 - rate.Vaccine/rate.Placebo)
  mod <- mutate_epi(mod, VE2.inst = 1 - rate.Vaccine.het/rate.Placebo.het)
  
  #VE <- 1 - Relative Risk; this is VE from cumulative incidence
  mod <- mutate_epi(mod, VE1.cumul = 1 - cumul.rate.Vaccine/cumul.rate.Placebo)
  mod <- mutate_epi(mod, VE2.cumul = 1 - cumul.rate.Vaccine.het/cumul.rate.Placebo.het)
  
  #return(mod)
}


#------------------------------------------------------------------------------
# sim execution
#------------------------------------------------------------------------------
runSim_Paul <- function(reac = c( "numExecution" = 1000 )) {
  
    #browser()
    #time <- c(180,360,540,720,900,1080)  # every 6 months for 3 years
    time <- 3*365; # End of the trial.

    ## Note here actually want to target the average incidence over time, rather than the incidence at the end of the trial, since with multiple risk groups there is waning expected. I have changed this in the ABC function f, below.
    placebo.incidence.target <- rep( 3.5, length( time ) )    # flat incidence of 3.5% per 100 person years

    ### Paul notes this is not going to work, if efficacy is waning you can't match a constant efficacy over time. We should match just one time, probably duration of 702 trial.
    VE.target <- rep(0.3, length( time ))
    target.stats <- data.frame(time, VE.target, placebo.incidence.target )

    run.and.compute.run.stats <- function (
      epsilon,   #per contact vaccine efficacy
      lambda,     #beta*c*prev,
      risk          #risk multiplier
    ) {
      
      # Paul added the other params to this (risk):
      param <- param.dcm(lambda = lambda, epsilon = epsilon, risk = risk )
      init <- init.dcm(Sp = 5000, Ip = 0,
                       Sv = 5000, Iv = 0,
                       Sph = 500, Iph = 0,    #placebo, high risk
                       Spm = 4000, Ipm = 0,   #placebo, medium risk
                       Spl = 500, Ipl = 0,    #placebo, low risk
                       Svh = 500, Ivh = 0,    #vaccine
                       Svm = 4000, Ivm = 0,   #vaccine
                       Svl = 500, Ivl = 0,    #vaccine
                       SIp.flow = 0, SIv.flow = 0, 
                       SIph.flow = 0, SIpm.flow = 0, SIpl.flow = 0,
                       SIvh.flow = 0, SIvm.flow = 0, SIvl.flow = 0)
      
      control <- control.dcm(nsteps = 365*3, new.mod = si_odePaul)
      mod <- dcm(param, init, control)
      #print( mod )
      
      mod.with.stats <- mod.manipulate( mod )
      #print( mod.with.stats )
      mod.with.stats.df <- as.data.frame( mod.with.stats )
      
      # homogeneous risk:
      #hom.VE <- mod.with.stats.df$VE1.inst[target.stats$time]
      # heterogeneous risk:
      het.VE <- mod.with.stats.df$VE2.inst[target.stats$time]
      
      #out <- mod.with.stats.df$rate.Placebo.het[target.stats$time]
      ## Paul changed the placebo incidence out stat as the mean over time up to each time in target.stats$time.
      .x <- mod.with.stats.df$rate.Placebo.het;
      het.meanPlaceboIncidence <-
        sapply( target.stats$time, function( .time ) { mean( .x[ 1:.time ] ) } );
      #out <- .x[ target.stats$time ];
      c( het.VE = het.VE, het.meanPlaceboIncidence = het.meanPlaceboIncidence );
    } # run.and.compute.run.stats (..)

    # Specify bounds for the initial conditions; these are used as priors.
    # In order of "x" in the above function.
    bounds <- list(#c(0.003, 0.008),      # beta
                    #c(60/365, 120/365))    # c
                   epsilon = c(1E-10, 1-(1E-10)),
    #               lambda = c(0.000005, 0.0001),  # lambda
                   lambda = c(5E-6, 1E-4),#1E-5),  # lambda
                   risk = c(1, 20))            # risk multiplier

    ## Trying to use optim. Best for 2 dimensions.
    # make.optim.fn <- function ( target ) {
    #     function( x ) {
    #         abs( mean( run.and.compute.run.stats( epsilon = x[ 1 ], lambda = x[ 2 ], risk = x[ 3 ] ) - as.numeric( target ) ) );
    #     }
    # }
    # 
    # .f <- make.optim.fn( target.stats[-1] )
    # 
    # lower.bounds <- sapply( bounds, function ( .bounds.for.x ) { .bounds.for.x[ 1 ] } );
    # upper.bounds <- sapply( bounds, function ( .bounds.for.x ) { .bounds.for.x[ 2 ] } );
    # optim( c( epsilon = 0.30, lambda = beta*c*prev, risk = 10 ), .f, lower = lower.bounds, upper = upper.bounds, method = "L-BFGS-B" )
    
    # With three params, optim fails here, presumably because it is nonconvex (multiple optima, saddles, as I have suspected; a potential option is reparameterization, I will think on that, but here we are going with ABC as suggested by Sam.

    ### ABC
    make.abc.fn <- function ( target ) {
        function( x ) {
            run.and.compute.run.stats( epsilon = x[ 1 ], lambda = x[ 2 ], risk = x[ 3 ] )
        }
    }

    # In order of "x" in the above function.
    priors  <- lapply( bounds, function( .bounds ) { c( "unif", unlist( .bounds ) ) } );

    .f.abc <- make.abc.fn( target.stats[-1] )
    fit.rej <- ABC_rejection(model = .f.abc,
                         prior = priors,
                         nb_simul = reac[ "numExecution" ],
                         summary_stat_target = as.numeric( target.stats[-1] ),
                         tol = 0.25,
                         progress_bar = TRUE)

    fit <- fit.rej

    return( fit );
} # runSim_Paul (..)

### Here is where Paul and Josh are working on August 23, 2021.
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
