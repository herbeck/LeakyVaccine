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
    list(c(
           dSph,dIph,
           dSpm,dIpm,
           dSpl,dIpl,
           dSvh,dIvh,
           dSvm,dIvm,
           dSvl,dIvl,
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
  mod <- mutate_epi(mod, rate.Vaccine.het = (total.SIvh.SIvm.SIvl.flow/total.Svh.Svm.Svl)*365*100)
  mod <- mutate_epi(mod, rate.Placebo.het = (total.SIph.SIpm.SIpl.flow/total.Sph.Spm.Spl)*365*100)
  
  #Cumulative incidence
  mod <- mutate_epi(mod, cumul.Svh.Svm.Svl = cumsum(total.Svh.Svm.Svl))
  mod <- mutate_epi(mod, cumul.Sph.Spm.Spl = cumsum(total.Sph.Spm.Spl))
  mod <- mutate_epi(mod, cumul.rate.Vaccine.het = (total.Ivh.Ivm.Ivl/cumul.Svh.Svm.Svl)*365*100)
  mod <- mutate_epi(mod, cumul.rate.Placebo.het = (total.Iph.Ipm.Ipl/cumul.Sph.Spm.Spl)*365*100)
  
  #Vaccine efficacy (VE) estimates
  #VE <- 1 - Relative Risk; this is VE for instantaneous incidence / hazard
  mod <- mutate_epi(mod, VE2.inst = 1 - rate.Vaccine.het/rate.Placebo.het)
  
  #VE <- 1 - Relative Risk; this is VE from cumulative incidence
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
runSim_Paul <- function(reac = c( "numExecution" = 10000, "numParams" = 3, "VE" = 0.1, "placeboIncidence" = 3.5 )) {
    stopifnot( all( c( "numExecution", "numParams", "VE", "placeboIncidence" ) %in% names( reac ) ) );

    ## Number of parameters to optimize (3, 4, or 5).
    num.params <- unname( reac[ "numParams" ] );

    num.sims <- unname( reac[ "numExecution" ] );
    placebo.incidence.target <- unname( reac[ "placeboIncidence" ] ); # incidence per 100 person years, eg 3.5 for 702 or 0.3 for RV144 [todo: double-check these values!]
    VE.target = unname( reac[ "VE" ] ); # instantaneous VE observed by the end of the trial

    #browser()

    ## Note that we target the average incidence over time, rather than the incidence at the end of the trial, since with multiple risk groups there is waning expected.

    # MAGIC NUMBERS
    trial.evaluation.time <- 3*365; # End of the trial.
    risk.max <- 50;
    lambda.min <- 1E-6;
    lambda.max <- 1E-3;
    smallest.discernable.amount <- 5E-4; # determined by trial and error this is the smallest amount you can change the parameters from 0 or 1 for it to register a difference from those extremes, eg to avoid NaN and Inf

    target.stats <- data.frame( trial.evaluation.time, VE.target, placebo.incidence.target );
    nsteps <- trial.evaluation.time;

    run.and.compute.run.stats <- function (
      epsilon,   #per contact vaccine efficacy
      lambda,     #beta*c*prev,
      risk,          #risk multiplier for high risk group
      highRiskProportion = 0.1, # these have defaults, so you can run this as 3 or 4 or 5 parameters.
      lowRiskProportion = ( 0.1/(1-highRiskProportion) ), # Among those not high risk.
      vaccinatedProportion = 0.5,
      trialSize = 10000
    ) {
      
      # Paul added the other params to this (risk here, others below):
      param <- param.dcm(lambda = lambda, epsilon = epsilon, risk = risk )

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
                       SIvh.flow = 0, SIvm.flow = 0, SIvl.flow = 0)
      
      control <- control.dcm(nsteps = nsteps, new.mod = si_odePaul)
      mod <- dcm(param, init, control)
      #print( mod )
      
      mod.with.stats <- mod.manipulate( mod );
      #print( mod.with.stats )
      mod.with.stats.df <- as.data.frame( mod.with.stats );
      
      # OLD: heterogeneous risk using instantaneous VE:
      #VE <- mod.with.stats.df$VE2.inst[target.stats$trial.evaluation.time];
      # NEW: heterogeneous risk using cumulative VE:
      VE <- mod.with.stats.df$VE2.cumul[target.stats$trial.evaluation.time];
      
      ## The placebo incidence out stat vector is the mean over time up to each time in target.stats$trial.evaluation.time.
      .x <- mod.with.stats.df$rate.Placebo.het;
      meanPlaceboIncidence <-
        sapply( target.stats$trial.evaluation.time, function( .time ) { mean( .x[ 1:.time ] ) } );
      #out <- .x[ target.stats$trial.evaluation.time ];
      c( VE = VE, meanPlaceboIncidence = meanPlaceboIncidence );
    } # run.and.compute.run.stats (..)

    # Specify bounds for the initial conditions; these are used as priors.
    # In order of "x" in the above function.
    bounds <- list(
                   epsilon = c(smallest.discernable.amount, 1-smallest.discernable.amount), # This is just the full range 0 to 1
                   lambda = c(lambda.min, lambda.max),
                   risk = c(1, risk.max), # risk multiplier for high risk group
                   highRiskProportion = c(smallest.discernable.amount, 1-smallest.discernable.amount), # This is just the full range 0 to 1
                   lowRiskProportion = c(smallest.discernable.amount, 1-smallest.discernable.amount) # This is modified to avoid NaN/Inf
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
            run.and.compute.run.stats( epsilon = x[ 1 ], lambda = x[ 2 ], risk = x[ 3 ], highRiskProportion = x[ 4 ], lowRiskProportion = x[ 5 ] )
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
                         nb_simul = num.sims,
                         summary_stat_target = as.numeric( target.stats[-1] ),
                         tol = 0.10, # Just keep top 10%
                         progress_bar = TRUE)

    # Compute the distances for each retained sample
    fit.rej.dist <-
        calculate.abc.dist( fit.rej$stats, as.numeric( target.stats[-1] ), ifelse( is.na( fit.rej$stats_normalization ), 1, fit.rej$stats_normalization ) );

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
    if( ncol( fit.rej$param ) < 5 ) {
        cl <- suppressWarnings( pdfCluster( fit.rej$param, bwtype="adaptive", hmult=1 ) );
    } else {
        cl <- pdfCluster( fit.rej$param, bwtype="adaptive", hmult=1, n.grid=nrow(fit.rej$param) );
    }

     cluster.numbers <- groups( cl );
     boxplot( fit.rej.dist ~ cluster.numbers )
    
    # Cluster-specific best sample is at this index, for each cluster:
    cluster.member.minimizing.dist <- sapply( 1:max( cluster.numbers ), function( .cluster.number ) { cluster.indices <- which( cluster.numbers == .cluster.number ); return( cluster.indices[ which.min( fit.rej.dist[ cluster.indices ] ) ] ); } );
    
    make.optim.fn <- function ( .target.stats, .target.stat.stdevs = rep( 1, length( .target.stats ) ) ) {
        function( x ) {
            .stats <- .f.abc( x );
            .stats.matrix <- matrix( .stats, nrow = 1 );
            c( dist = calculate.abc.dist( .stats.matrix, .target.stats, ifelse( is.na( .target.stat.stdevs ), 1, .target.stat.stdevs ) ) )
        }
    } # make.optim.fn (..)

    .f <- make.optim.fn( target.stats[-1], fit.rej$stats_normalization );

    # These are the original bounds, separated for use in optim:
    lower.bounds <- sapply( bounds, function ( .bounds.for.x ) { .bounds.for.x[ 1 ] } );
    upper.bounds <- sapply( bounds, function ( .bounds.for.x ) { .bounds.for.x[ 2 ] } );

    optima.by.cluster <- sapply( cluster.member.minimizing.dist[1:2], function( .minimizer.index ) {
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
} # runSim_Paul (..)

### ERE I AM testing...
the.seed <- 98103;
# To test replicability of the identified modes, uncomment this:
# set.seed( the.seed ); the.seed <- floor( runif( 1, max = 1E5 ) );
num.sims <- 1000; # Fast for debugging.
# num.sims <- 10000; # For reals.

set.seed( the.seed );

## Ok, I think that we need to first determine the reasonable epsilon values from the RV144-like setting, then ask could you get zero (or very low) VE with the same epsilon by increasing only lambda and risk values, with placebo incidence going from an rv144-like value (~0.14) to a 702-like value (~3.3).
## rv144 placebo incidence in the prespecified analysis (MITT) cohort was 0.1397 per 100 person years, 100*(74/52985) from N Engl J Med 2009; 361:2209-2220 DOI: 10.1056/NEJMoa0908492 December 3, 2009 (https://www.nejm.org/doi/full/10.1056/nejmoa0908492):
# "HIV-1 infection was diagnosed in 132 subjects (56 in the vaccine group and 76 in the placebo group) during 52,985 person-years of follow-up in the intention-to-treat analysis, in 86 subjects (36 in the vaccine group and 50 in the placebo group) during 36,720 person-years of follow-up in the per-protocol analysis, and in 125 subjects (51 in the vaccine group and 74 in the placebo group) during 52,985 person-years of follow-up in the modified intention-to-treat analysis. One subject in the placebo group who was identified by hospital record as being seropositive for HIV after dying from Pneumocystis jirovecii pneumonia was included in the analysis before the unblinding of the study. This diagnosis of HIV-1 infection was the only one that occurred outside planned procedures."

## hvtn 702 placebo incidence was 3.3 per 100 person-years (95% CI, 2.8 to 3.9), from n engl j med 384;12 nejm.org March 25, 2021 (https://www.nejm.org/doi/pdf/10.1056/NEJMoa2031499), page 1092:
# "During the first 24 months of follow-up, 138 HIV-1 infections occurred in the vaccine group and 133 in the placebo group, for an estimated incidence rate of 3.4 per 100 person-years (95% confidence interval [CI], 2.8 to 4.0) and 3.3 per 100 person-years (95% CI, 2.8 to 3.9), respectively (hazard ratio, 1.02; 95% CI, 0.81 to 1.30; P=0.84) (Fig. 1A and Table 2). The incidence of HIV-1 infection was similar in the vaccine group and the placebo group in secondary analyses during 36 months of follow-up (hazard ratio, 1.05; 95% CI, 0.83 to 1.31), in the month 6.5 at-risk cohort between 6.5 months and 24 months (hazard ratio, 1.15; 95% CI, 0.84 to 1.58), and in the perprotocol cohort, as well as in other secondary analyses (Figs. S3 through S9)."

## FOR RV144:
.sim3 <- runSim_Paul( reac = c( "numExecution" = num.sims, "numParams" = 3, "VE" = 0.31, "placeboIncidence" = 0.14 ));
# .sim4 <- runSim_Paul( reac = c( "numExecution" = num.sims, "numParams" = 4, "VE" = 0.31, "placeboIncidence" = 0.14 ));
# .sim5 <- runSim_Paul( reac = c( "numExecution" = num.sims, "numParams" = 5, "VE" = 0.31, "placeboIncidence" = 0.14 ));


## FOR 702:
# .sim3 <- runSim_Paul( reac = c( "numExecution" = num.sims, "numParams" = 3, "VE" = 0, "placeboIncidence" = 3.3 ));
# .sim4 <- runSim_Paul( reac = c( "numExecution" = num.sims, "numParams" = 4, "VE" = 0, "placeboIncidence" = 3.3 ));
# .sim5 <- runSim_Paul( reac = c( "numExecution" = num.sims, "numParams" = 5, "VE" = 0, "placeboIncidence" = 3.3 ));

######
## Some plotting. Run manually. See above.
if( FALSE ) {
    the.sim <- .sim3;

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
    legend("topright", legend = c("VE target", "Posterior"),
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

