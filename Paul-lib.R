library(deSolve)
library(tidyverse)
library(EpiModel)
library(survival)
library(EasyABC)
library(ggplot2)
library(pdfCluster)

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
runSim_Paul <- function(reac = c( "numExecution" = 1000, "numParams" = 3 )) {
  
    ## Number of parameters to optimize (3, 4, or 5).
    num.params <- reac[ "numParams" ];

    #browser()
    #time <- c(180,360,540,720,900,1080)  # every 6 months for 3 years
    time <- 3*365; # End of the trial.
    nsteps <- time;
    ## Note here actually want to target the average incidence over time, rather than the incidence at the end of the trial, since with multiple risk groups there is waning expected. I have changed this in the ABC function f, below.
    placebo.incidence.target <- rep( 3.5, length( time ) )    # flat incidence of 3.5% per 100 person years

    ### Paul notes this is not going to work, if efficacy is waning you can't match a constant efficacy over time. We should match just one time, probably duration of 702 trial.
    VE.target <- rep(0.3, length( time ))
    target.stats <- data.frame(time, VE.target, placebo.incidence.target )

    run.and.compute.run.stats <- function (
      epsilon,   #per contact vaccine efficacy
      lambda,     #beta*c*prev,
      risk,          #risk multiplier for high risk group
      highRiskProportion = 0.1,
      immuneProportion = 0.1,
      vaccinatedProportion = 0.5,
      trialSize = 10000
    ) {
      
      # Paul added the other params to this (risk here, others below):
      param <- param.dcm(lambda = lambda, epsilon = epsilon, risk = risk )

      Svl <- floor( immuneProportion * vaccinatedProportion * trialSize );
      Spl <- floor( immuneProportion * ( 1.0 - vaccinatedProportion ) * trialSize );
      Svh <- floor( highRiskProportion * vaccinatedProportion * trialSize );
      Sph <- floor( highRiskProportion * ( 1.0 - vaccinatedProportion ) * trialSize );
      Svm <- floor( ( 1.0 - ( immuneProportion + highRiskProportion ) ) * vaccinatedProportion * trialSize );
      Spm <- floor( ( 1.0 - ( immuneProportion + highRiskProportion ) ) * ( 1.0 - vaccinatedProportion ) * trialSize );

      Sp <- Spl + Spm + Sph;
      Sv <- Svl + Svm + Svh;
      if( vaccinatedProportion == 0.5 ) {
          stopifnot( Sp == Sv );
      }
      init <- init.dcm(Sp = Sp, Ip = 0,
                       Sv = Sv, Iv = 0,
                       Sph = Sph, Iph = 0,    #placebo, high risk
                       Spm = Spm, Ipm = 0,   #placebo, medium risk
                       Spl = Spl, Ipl = 0,    #placebo, low risk
                       Svh = Svh, Ivh = 0,    #vaccine
                       Svm = Svm, Ivm = 0,   #vaccine
                       Svl = Svl, Ivl = 0,    #vaccine
                       SIp.flow = 0, SIv.flow = 0, 
                       SIph.flow = 0, SIpm.flow = 0, SIpl.flow = 0,
                       SIvh.flow = 0, SIvm.flow = 0, SIvl.flow = 0)
      
      control <- control.dcm(nsteps = nsteps, new.mod = si_odePaul)
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
    bounds <- list(
                   epsilon = c(1E-10, 1-(1E-10)),
                   lambda = c(1E-6, 1E-3),#1E-5),  # lambda
                   risk = c(1, 50), # risk multiplier for high risk group
                   highRiskProportion = c(1E-10, 1-(1E-10)),
                   immuneProportion = c(1E-10, 1-(1E-10))
                   );            
    # In order of "x" in the above function.
    priors  <- lapply( bounds, function( .bounds ) { c( "unif", unlist( .bounds ) ) } );

    ### ABC
    # For the three-parameter version use this one, and change above too.
    make.abc.fn.3params <- function ( target ) {
        function( x ) {
            run.and.compute.run.stats( epsilon = x[ 1 ], lambda = x[ 2 ], risk = x[ 3 ] )
        }
    }
    # For the four-parameter version use this one, and change above too.
    make.abc.fn.4params <- function ( target ) {
        function( x ) {
            run.and.compute.run.stats( epsilon = x[ 1 ], lambda = x[ 2 ], risk = x[ 3 ], highRiskProportion = x[ 4 ] )
        }
    }
    # For the five-parameter version use this one, and change above too.
    make.abc.fn.5params <- function ( target ) {
        function( x ) {
            run.and.compute.run.stats( epsilon = x[ 1 ], lambda = x[ 2 ], risk = x[ 3 ], highRiskProportion = x[ 4 ], immuneProportion = x[ 5 ] )
        }
    }

    bounds <- bounds[ 1:num.params ];
    priors <- priors[ 1:num.params ];
    if( num.params == 3 ) {
        .f.abc <- make.abc.fn.3params( target.stats[-1] )
    } else if( num.params == 4 ) {
        .f.abc <- make.abc.fn.4params( target.stats[-1] )
    } else {
        .f.abc <- make.abc.fn.5params( target.stats[-1] )
    }
    fit.rej <- ABC_rejection(model = .f.abc,
                         prior = priors,
                         nb_simul = reac[ "numExecution" ],
                         summary_stat_target = as.numeric( target.stats[-1] ),
                         tol = 0.25,
                         progress_bar = TRUE)

    # Compute the distances for each retained sample
    fit.rej.dist <-
        calculate.abc.dist( fit.rej$stats, as.numeric( target.stats[-1] ), fit.rej$stats_normalization );

    # Now find optimal points near the sampled modes.
    
    # Globally best sample is at this index:
    # sample.index.minimizing.dist <- which.min( fit.rej.dist );
    
    cl <- pdfCluster( fit.rej$param );
    cluster.numbers <- groups( cl );
    # boxplot( fit.rej.dist ~ cluster.numbers )
    
    # Cluster-specific best sample is at this index, for each cluster:
    cluster.member.minimizing.dist <- sapply( 1:max( cluster.numbers ), function( .cluster.number ) { cluster.indices <- which( cluster.numbers == .cluster.number ); return( cluster.indices[ which.min( fit.rej.dist[ cluster.indices ] ) ] ); } );
    
    make.optim.fn <- function ( target.stats, target.stat.stdevs = rep( 1, length( target.stats ) ) ) {
        function( x ) {
            .stats <- .f.abc( x );
            .stats.matrix <- matrix( .stats, nrow = 1 );
            c( dist = calculate.abc.dist( .stats.matrix, target.stats, target.stat.stdevs ) )
        }
    } # make.optim.fn (..)

    .f <- make.optim.fn( target.stats[-1], fit.rej$stats_normalization );

    # These are the original bounds, separated for use in optim:
    lower.bounds <- sapply( bounds, function ( .bounds.for.x ) { .bounds.for.x[ 1 ] } );
    upper.bounds <- sapply( bounds, function ( .bounds.for.x ) { .bounds.for.x[ 2 ] } );

    optima.by.cluster <- sapply( cluster.member.minimizing.dist, function( .minimizer.index ) {
        print( .minimizer.index );
        print( fit.rej$param[ .minimizer.index, ] );
        print( .f( fit.rej$param[ .minimizer.index, ] ) );
        .optim.result <- optim( fit.rej$param[ .minimizer.index, ], .f, lower = lower.bounds, upper = upper.bounds, method = "L-BFGS-B" );
        print( c( param = .optim.result$par, dist = .optim.result$value ) );
        return( c( param = .optim.result$par, dist = .optim.result$value ) );
    } );
    rownames( optima.by.cluster )[ 1:length( bounds ) ] <- names( bounds );
    optima.by.cluster.sorted <- optima.by.cluster[ , order( optima.by.cluster[ "dist", ] ), drop = FALSE ];

    return( list( fit = fit.rej, priors = priors, target.stats = target.stats, fn = .f.abc, sampled.modes = optima.by.cluster.sorted ) );
} # runSim_Paul (..)

### ERE I AM testing...
# set.seed( 98103 );
# .sim3 <- runSim_Paul( reac = c( "numExecution" = 10000, "numParams" = 3 ));

# set.seed( 98103 );
# .sim4 <- runSim_Paul( reac = c( "numExecution" = 10000, "numParams" = 4 ));

# set.seed( 98103 );
# .sim5 <- runSim_Paul( reac = c( "numExecution" = 10000, "numParams" = 5 ));


