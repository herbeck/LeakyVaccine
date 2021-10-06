source( "Paul-lib-threegroup.R" )

hvtn705.parameters <- paste( "hvtn705", trial.parameters, sep = "." );
all.parameters <- c( common.parameters, hvtn705.parameters );

# Note we can now have targets with multiple end times for VE, with the placebo incidence target associated with the same end time - but NOTE that the VE is evaluated from VE.start whereas the cumulative placebo incidence is evaluated over the whole trial (from time 0).
hvtn705.VE <- 25.2;
hvtn705.placeboIncidence <- 4.32;

#------------------------------------------------------------------------------
# sim execution
#------------------------------------------------------------------------------
run.and.compute.run.stats.hvtn705 <- function (
      epsilon,   #per contact vaccine efficacy
      hvtn705.log10lambda,
      hvtn705.log10riskmultiplier,
      hvtn705.highRiskProportion,
      hvtn705.lowRiskProportion,
      hvtn705.VE.start = ceiling( 365 / 2 ), # start of 7th month
      hvtn705.VE.end = ( 2 * 365 ), # end of 24th month
      vaccinatedProportion = 0.5,  # In lieu of naming vaccine and placebo arms separately (and their N)
      trialSize = 2200,  # Size of HVTN 705, 1:1 allocation
      sim.ndays = hvtn705.VE.end
) {
    hvtn705.results <- run.and.compute.run.stats.threegroup(
      epsilon = epsilon,   # common
      log10lambda = hvtn705.log10lambda,
      log10riskmultiplier = hvtn705.log10riskmultiplier,
      highRiskProportion = hvtn705.highRiskProportion,
      lowRiskProportion = hvtn705.lowRiskProportion,
      vaccinatedProportion = vaccinatedProportion,
      VE.start = hvtn705.VE.start,
      VE.end = hvtn705.VE.end,
      trialSize = trialSize,
      sim.ndays = sim.ndays
    );
    return( list( hvtn705 = hvtn705.results ) );
} # run.and.compute.run.stats.hvtn705 (..)

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

# Ok, so there's one set of model-specific parameters (
# hvtn705), and then one parameter that is shared (epsilon, the
# per-contact vaccine efficacy parameter), so we
# optimize this using conditional optimization. That is, we hold
# epsilon fixed, and then do model-specific params optimization separately,
# then hold those params fixed and optimize epsilon, and iterate
# until convergence.
make.epsilon.abc.fn <- function ( other.parameters ) {
    function( x ) {
        run.and.compute.run.stats.hvtn705( epsilon = x, hvtn705.log10lambda = other.parameters[[ "hvtn705.log10lambda" ]], hvtn705.log10riskmultiplier = other.parameters[[ "hvtn705.log10riskmultiplier" ]], hvtn705.highRiskProportion = other.parameters[[ "hvtn705.highRiskProportion" ]], hvtn705.lowRiskProportion = other.parameters[[ "hvtn705.lowRiskProportion" ]] )
    }
}
make.hvtn705.abc.fn <- function ( other.parameters ) {
    function( x ) {
        run.and.compute.run.stats.hvtn705( epsilon = other.parameters[[ "epsilon" ]], hvtn705.log10lambda = x[ 1 ], hvtn705.log10riskmultiplier = x[ 2 ], hvtn705.highRiskProportion = x[ 3 ], hvtn705.lowRiskProportion = x[ 4 ] )
    }
}
 
make.epsilon.optim.fn <- function ( .f.epsilon.abc, .target.stats, .target.stat.stdevs = rep( 1, length( .target.stats ) ) ) {
    function( x ) {
        .stats <- .f.epsilon.abc( x );
        .stats.matrix <- matrix( unlist( .stats ), nrow = 1 );
        c( dist = unname( calculate.abc.dist( .stats.matrix, .target.stats, ifelse( is.na( .target.stat.stdevs ), 1, .target.stat.stdevs ) ) ) )
    }
} # make.epsilon.optim.fn (..)
 
make.hvtn705.optim.fn <- function ( .f.hvtn705.abc, .target.stats, .target.stat.stdevs = rep( 1, length( .target.stats ) ) ) {
    function( x ) {
        .stats <- .f.hvtn705.abc( x );
        .stats.matrix <- matrix( unlist( .stats ), nrow = 1 );
        c( dist = calculate.abc.dist( .stats.matrix, .target.stats, ifelse( is.na( .target.stat.stdevs ), 1, .target.stat.stdevs ) ) )
    }
} # make.hvtn705.optim.fn (..)
 
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
              Tukey.IQR.multiplier,
              be.verbose = FALSE
              ) {
    trial.target.stat.stdevs <-
        ifelse( is.na( the.fit.bin$stats_normalization[ trial.target.stats.names ] ), 1, the.fit.bin$stats_normalization[ trial.target.stats.names ] ) * target.stat.scale.units[ trial.target.stats.names ];
    the.fit.trial.dist <-
        calculate.abc.dist( the.fit.bin$stats[ , trial.target.stats.names ], as.numeric( target.stats[ trial.target.stats.names ] ), trial.target.stat.stdevs );

    # Scale distance by units defined by the max distance of the closest units.quantile fraction of the points.
    trial.dist.units <- quantile( the.fit.trial.dist, probs = units.quantile )
    the.fit.trial.dist.scaled <- the.fit.trial.dist / trial.dist.units;

    abc.trial.keep.sim <- the.fit.trial.dist.scaled < max.dist.scaled;
    if( be.verbose ) {
        cat( paste( "Keeping ", sum( abc.trial.keep.sim ), " trial-specific parameter sets because they are within ", max.dist.scaled, " scaled units on the trial-specific distance measure, where one unit has ", sprintf( "%0.2f", 100*units.quantile ), "% of the original ", nrow( the.fit.bin$param ), " draws in this epsilon bin.", sep = "" ), fill = TRUE );
    }

    the.fit.trial <- the.fit.bin;
    the.fit.trial$param <- the.fit.trial$param[ abc.trial.keep.sim, c( "epsilon", trial.parameters.names ), drop = FALSE ];
    the.fit.trial$stats <- the.fit.trial$stats[ abc.trial.keep.sim, trial.target.stats.names, drop = FALSE ];
    the.fit.trial$weights <- the.fit.trial$weights[ abc.trial.keep.sim ];
    the.fit.trial.dist <- the.fit.trial.dist[ abc.trial.keep.sim ];

    cl.trial <- suppressWarnings( pdfCluster( the.fit.trial$param[ , trial.parameters.names, drop = FALSE ], bwtype="adaptive", hmult=pdfCluster.hmult, n.grid=nrow( the.fit.trial$param ) ) );
    trial.cluster.numbers <- groups( cl.trial );
    if( be.verbose ) {
        # Helpful when debugging - run through the next line.. it'll print when done
        print( table( trial.cluster.numbers ) );
    }

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

# For better balancing of the costs we use target.stat.scale.units; see above.
optimize.step <- function ( current.parameters, target.stats, target.stat.stdevs, lower, upper, current.value = NULL, be.verbose = FALSE ) {
    is.anything.changed <- FALSE;
      
    ## First, update epsilon
    .f.epsilon.abc <- make.epsilon.abc.fn( current.parameters );
    .f.epsilon <- make.epsilon.optim.fn( .f.epsilon.abc, target.stats, target.stat.stdevs );

    if( is.null( current.value ) ) {
        current.value <- unlist( .f.epsilon( current.parameters[[ "epsilon" ]] ) );
    }
    # Save the starting parameters and value, jic.
    starting.parameters <- current.parameters;
    starting.value <- current.value;

    if( be.verbose ) {
        cat( paste( apply( cbind( names( starting.parameters ), starting.parameters ), 1, paste, collapse = "=" ), collapse = ", " ), fill = TRUE )
        cat( starting.value, fill = TRUE );
    }

    .optim.result <- optim( current.parameters[[ "epsilon" ]], .f.epsilon, lower = lower[[ "epsilon" ]], upper = upper[[ "epsilon" ]], method = "L-BFGS-B" );
    .par <- .optim.result$par;
    names( .par ) <- "epsilon";

    .new.value <- .optim.result$value;
    if( be.verbose ) {
        cat( paste( "epsilon =", .par ), fill = TRUE );
        cat( .new.value, fill = TRUE );
    }

    # If they are printing the same to 6 digits, they are the same afaict
    # MAGIC # (6 digits)
    if( ( .new.value < current.value ) && ( sprintf( "%0.6f", .new.value ) != sprintf( "%0.6f", current.value ) ) ) {
        ## Update the current.parameters with this new epsilon value.
        if( be.verbose ) {
            cat( paste( "ACCEPT epsilon change from ", current.parameters[[ "epsilon" ]], " to ", .par, sep = "" ), fill = TRUE );
        }
        current.parameters[[ "epsilon" ]] <- .par;
        current.value <- .new.value;
        is.anything.changed <- TRUE;
    } else {
        if( be.verbose ) {
            cat( "REJECT epsilon change", fill = TRUE );
            ## TODO: REMOVE
            print( .optim.result );
        }
    }


    ## Next, update hvtn705 parameters
    .f.hvtn705.abc <- make.hvtn705.abc.fn( current.parameters );
    .f.hvtn705 <- make.hvtn705.optim.fn( .f.hvtn705.abc, target.stats, target.stat.stdevs );

    .optim.result <- optim( current.parameters[ hvtn705.parameters ], .f.hvtn705, lower = lower[ hvtn705.parameters ], upper = upper[ hvtn705.parameters ], method = "L-BFGS-B" );
    .par <- .optim.result$par;
    names( .par ) <- hvtn705.parameters;

    .new.value <- .optim.result$value;
    if( be.verbose ) {
        cat( paste( apply( cbind( hvtn705.parameters, .par ), 1, paste, collapse = "=" ), collapse = ", " ), fill = TRUE )
        cat( .new.value, fill = TRUE );
    }

    # If they are printing the same to 6 digits, they are the same afaict
    # MAGIC # (6 digits)
    if( ( .new.value < current.value ) && ( sprintf( "%0.6f", .new.value ) != sprintf( "%0.6f", current.value ) ) ) {
        if( be.verbose ) {
            cat( paste( "ACCEPT hvtn705 parameters change from ", paste( apply( cbind( hvtn705.parameters, current.parameters[ hvtn705.parameters ] ), 1, paste, collapse = "=" ), collapse = ", " ), " to ", paste( apply( cbind( hvtn705.parameters, .par ), 1, paste, collapse = "=" ), collapse = ", " ), sep = "" ), fill = TRUE );
        }
        ## Update the current.parameters with this new epsilon value.
        current.parameters[ hvtn705.parameters ] <- .par;
        current.value <- .new.value;
        is.anything.changed <- TRUE;
    } else {
        if( be.verbose ) {
            cat( "REJECT hvtn705 parameters change", fill = TRUE );
            ## TODO: REMOVE
            print( .optim.result );
        }
    }

    .f.epsilon.abc <- make.epsilon.abc.fn( current.parameters );
    .new.stats <- .f.epsilon.abc( current.parameters[[ "epsilon" ]] );

    if( !is.anything.changed ) {
        return( NULL );
    }
    if( be.verbose ) {
        print( c( current.parameters, .new.stats, dist = current.value ) );
    }
    return( c( current.parameters, .new.stats, dist = current.value ) );
} # optimize.step (..)
 
# First we optimize epsilon, then hvtn705, then back to epsilon again
#      ‘reltol’ Relative convergence tolerance.  The algorithm stops if
# it is unable to reduce the value by a factor of ‘reltol *
# (abs(val) + reltol)’ at a step.
optimize.iteratively <- function ( current.parameters, target.stats, target.stat.stdevs, lower, upper, current.value = NULL, reltol = 1E-2, step.i = 1, max.steps = 50, be.verbose = FALSE ) {
    .converged <- FALSE;
    last.dist <- current.value;
    while( !.converged && ( step.i < max.steps ) ) {
        if( be.verbose ) {
            cat( paste( "optimize.iteratively( step.i = ", step.i, " )", sep = "" ), fill = TRUE );
        }
        .step.i.result <-
            optimize.step( current.parameters = current.parameters, target.stats = target.stats, target.stat.stdevs = target.stat.stdevs, lower = lower, upper = upper, current.value = current.value, be.verbose = be.verbose );
        if( is.null( .step.i.result ) ) {
            .converged <- TRUE;
        } else {
            current.parameters <- .step.i.result[ all.parameters ];
            current.value <- .step.i.result[[ "dist" ]];
            current.stats <-
                .step.i.result[ setdiff( names( .step.i.result ), c( all.parameters, "dist" ) ) ];
            if( !is.null( last.dist ) ) {
                if( last.dist == current.value ) {
                    .converged <- TRUE;
                } else {
                    .converged <- abs( last.dist - current.value ) < ( reltol * ( abs( last.dist ) + reltol ) )
                }
            } else {
                last.dist <- current.value;
            }
        }
        step.i <- step.i + 1;
    } # End while !.converged && step.i < max.steps
    if( be.verbose ) {
        if( .converged ) {
            cat( "CONVERGED.", fill = TRUE );
            cat( paste( apply( cbind( names( current.parameters ), current.parameters ), 1, paste, collapse = "=" ), collapse = ", " ), fill = TRUE )
            cat( current.value, fill = TRUE );
        } else {
            cat( "DID NOT CONVERGE (max steps reached).", fill = TRUE );
        }
    }
    return( c( current.parameters, unlist( current.stats ), dist = unname( current.value ) ) );
} # optimize.iteratively (..)

# turn the input params (eg from enclosing Shiny app), for some reason called "reac" (maybe something to do with Shiny?), into the full parameter set that we need in order to run the analysis.
runSim_hvtn705.create.params.list <- function ( reac ) {
    stopifnot( all( c( "numExecution" ) %in% names( reac ) ) );

    num.sims <- unname( reac[ "numExecution" ] );

    ######################################################################
    ##### PARAMETERS / MAGIC #s -- MANY COULD BE EXPOSED FOR TUNING.
    if( num.sims <= 1000 ) {
        # Number of samples to run through the clustering and optimizing steps.
        abc.keep.num.sims <- floor( num.sims / 4 ); # MAGIC #
    } else {
        # Number of samples to run through the clustering and optimizing steps.
        abc.keep.num.sims <- 2500; # MAGIC #
    }
    stopifnot( num.sims > abc.keep.num.sims ); # It won't work to cluster uniformly drawn points. You first have to filter them by keeping those nearest the target.

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
    log10riskmultiplier.max <- log10( 50 ); # Note that log10riskmultiplier.max must be < 1/(10^log10lambda.max) to ensure that the flow does not go negative within the ode - that just means log10riskmultiplier+log10lambda must be < 0; see ode def in Paul-lib-threegroup.R.
    log10lambda.min <- -7;
    log10lambda.max <- -3;
    # Note restriction - see note above on log10riskmultiplier.max
    stopifnot( ( log10riskmultiplier.max + log10lambda.max ) < 0 );

    # MAGIC # but this doesn't really seem to be a tunable parameter, but I basically just set it arbitrarily to a value that I think is enough to get some change; it might need tuning in future.
    smallest.discernable.amount <- 5E-4; # determined by trial and error this is the smallest amount you can change the parameters from 0 or 1 for it to register a difference from those extremes, eg to avoid NaN and Inf

    # We keep at least a fraction of ( abc.keep.num.sims / num.sims ) sims within each cluster, but we actually go a little further to help ensure we are not getting just isolated peaks but also some of the points nearby them to help estimate the IQRs of the clusters for the subsequent search. This is the scale factor by which we multiple the distance; so for example if we are keeping 2.5% of the samples, and the 0.025 quantile of trial-specific distances in a bin for rv144 is 0.67, then we would actually possibly keep more than 2.5% of the samples because we would keep all that samples fall within max.dist.scaled (which should be >= 1) times that distance (eg 1.05 * 0.67). Note that we separately keep samples in this manner for each trial for trial-specific clustering of the bin-specific samples.
    max.dist.scaled <- 1.05; # MAGIC #

    # Once we identify modes we determine their compatibility across studies by whether the epsilons are near enough to each other, based on the overlap or not of their search windows, which are defined by using the logic of identifying outliers in a box plot. That is we walk out from the mode some number of IQR-defined units and define the search window as anything within that zone. For Tukey-style boxplots the number is 1.5 IQRs from the median defines an outlier. Here we can tune this multiplier (Tukey.IQR.multiplier) to change the window width that we use to define the search space for optimization and also for this overlap detection.
    Tukey.IQR.multiplier <- 1.5; # MAGIC # governing the window size we use when optimizing and when merging trial-specific candidate optima (this is applied to epsilon only for merging).

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

    optimize.iteratively.max.steps <- 50; # MAGIC # governs limit on times optimize.iteratively(..) is called.
    optimize.iteratively.reltol <- 1E-2; # MAGIC # governs convergence criteria for optimize.iteratively(..) (when the dist target changes by less than this fraction in one iteration, the optimization is deemed to have converged).

    pdfCluster.hmult <- 1.05; # MAGIC #, tweaked it to get pdfCluster to run without crashing with the error message suggesting increasing n.grid -- even with max n.grid, hmult has to be just above 1, it seems. If it is too high, the clusters merge into one. Tune this (> 1 to prevent crashing) higher to get the clusters to merge more.

    be.verbose <- TRUE; # MAGIC # (just governs output printed to screen)
    ######################################################################


    ######################################################################
    # Code below here has no additional tunable parameters...
    ######################################################################
    hvtn705.placebo.incidence.target <- hvtn705.placeboIncidence; # incidence per 100 person years
    hvtn705.VE.target <- hvtn705.VE; # cumulative VE observed by the end of the trial

    target.stats <-
        c( hvtn705.VE.target = hvtn705.VE.target, hvtn705.placebo.incidence.target = hvtn705.placebo.incidence.target );        
    hvtn705.target.stats.names <- grep( "hvtn705", names( target.stats ), value = TRUE );

    ## Note that we target the cumulative incidence at the end of the trial.
    target.stat.scale.units <- sapply( names( target.stats ), function( .target.stat ) { ifelse( length( grep( "VE", .target.stat ) ) > 0, VE.target.scale.units, ifelse( length( grep( "placebo.incidence.target", .target.stat ) ) > 0, placebo.incidence.target.scale.units, NA ) ) } );
    stopifnot( all( !is.na( target.stat.scale.units ) ) );

    # Specify bounds for the initial conditions; these are used as priors.
    # In order of "x" in the above function.
    bounds <- list(
                   epsilon = c(smallest.discernable.amount, 1-smallest.discernable.amount), # This is just the full range 0 to 1
                   hvtn705.log10lambda = c(log10lambda.min, log10lambda.max),
                   hvtn705.log10riskmultiplier = c(log10(1), log10riskmultiplier.max), # log10( risk multiplier for high risk group )
                   hvtn705.highRiskProportion = c(smallest.discernable.amount, 0.5-smallest.discernable.amount), # This is the half range 0 to 0.5
                   hvtn705.lowRiskProportion = c(smallest.discernable.amount, 1-smallest.discernable.amount) # This is modified to avoid NaN/Inf
                   );            
    stopifnot( all( names( bounds ) == all.parameters ) );

    # In order of "x" in the above function.
    priors  <- lapply( bounds, function( .bounds ) { c( "unif", unlist( .bounds ) ) } );
    return( as.list( environment() ) );
} # runSim_hvtn705.create.params.list ()

# uses params, so do within "with( runSim_hvtn705.create.params.list( reac ), .."
draw.from.priors <- function ( .f.abc ) {
    ## PHASE 1: Find candidate complete 9-parameter starting places constructed from merging epsion-bin-specific, trial-specific local optima that share common ranges of epsilon across the two trials. Later (in phase 2) we will find 9-parameter local optima near each of these candidate starting places.
    ## First we draw num.sims draws, then for each candidate epsilon bin we compute the study-specific distances for points falling in that bin (computed using just that study's two target stats), and cluster all points within 1.05-fold (see max.dist.scaled) of the furthest of the set of the top (abc.keep.num.sims/num.sims) (eg 2.5%) of samples within that bin, to ensure a minimum number of points and the extra points up to max.dist.scaled ensures that we keep additional points to help flesh out the contours of the space just below these peaks. This might possibly help with the clustering but it's not entirely clear yet how sensitive things are to max.dist.scaled.

    # First we just draw num.sims draws from the priors, independently. Keep everything drawn. Compute the 4-parameter target stats from runs with these 9-parameter random starting places, and the standard deviations of these 4-parameter stats (which we use to scale the distance function on the target stats for balancing the optimization evenly across the parameters, see below).
    fit.rej <- ABC_rejection(
                             model = .f.abc,
                             prior = priors,
                             nb_simul = num.sims,
                             summary_stat_target = as.numeric( target.stats ),
                             tol = 1.0, # Keep all of them for the moment; see below.
                             progress_bar = be.verbose
               );
    colnames( fit.rej$param ) <- all.parameters;
    colnames( fit.rej$stats ) <- names( target.stats );
    names( fit.rej$stats_normalization ) <- names( target.stats );

    return( fit.rej );
} # draw.from.priors ( .. )

# Define functions that operate on hvtn705.Tukey.whisker.bounds.by.trial.cluster:
# Once we find a pair with overlapping epsilon windows (one from each trial) we use this to add a new set of bounds to the candidate.parameter.sets.high and candidate.parameter.sets.low values, which are returned by this.
add.candidate.parameter.set <-
    function ( candidate.parameter.sets, hvtn705.cluster.j, hvtn705.Tukey.whisker.bounds.by.trial.cluster ) {
    if( merge.bounds.strategy.intersect ) {
        low.epsilon.bound <-
            hvtn705.Tukey.whisker.bounds.by.trial.cluster[[ "low" ]][ "epsilon", hvtn705.cluster.j ];
        high.epsilon.bound <-
            hvtn705.Tukey.whisker.bounds.by.trial.cluster[[ "high" ]][ "epsilon", hvtn705.cluster.j ];
    } else { # if merge.bounds.strategy.intersect .. else ..
        # use the strategy "union" instead of "intersection"
        low.epsilon.bound <-
            hvtn705.Tukey.whisker.bounds.by.trial.cluster[[ "low" ]][ "epsilon", hvtn705.cluster.j ];
        high.epsilon.bound <-
            hvtn705.Tukey.whisker.bounds.by.trial.cluster[[ "high" ]][ "epsilon", hvtn705.cluster.j ];
    } # End if merge.bounds.strategy.intersect .. else ..
    stopifnot( high.epsilon.bound > low.epsilon.bound );
    candidate.parameter.sets.low <-
        rbind( candidate.parameter.sets[[ "low" ]],
              c( "epsilon" = unname( low.epsilon.bound ),
                hvtn705.Tukey.whisker.bounds.by.trial.cluster[[ "low" ]][ hvtn705.parameters, hvtn705.cluster.j ]
                ) );
    candidate.parameter.sets.high <-
        rbind( candidate.parameter.sets[[ "high" ]],
              c( "epsilon" = unname( high.epsilon.bound ),
                hvtn705.Tukey.whisker.bounds.by.trial.cluster[[ "high" ]][ hvtn705.parameters, hvtn705.cluster.j ]
                ) );
    
    return( list( "low" = candidate.parameter.sets.low, "high" = candidate.parameter.sets.high ) );
} # add.candidate.parameter.set (..)

# uses params, so do within "with( runSim_hvtn705.create.params.list( reac ), .."
get.candidate.parameter.sets <- function ( fit.rej ) {
    target.stat.stdevs <-
        ifelse( is.na( fit.rej$stats_normalization ), 1, fit.rej$stats_normalization ) * target.stat.scale.units;

    candidate.parameter.sets <-
        list(
             "low" = matrix( NA, nrow = 0, ncol = length( all.parameters ) ),
             "high" = matrix( NA, nrow = 0, ncol = length( all.parameters ) )
            );
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

        ### hvtn705
        hvtn705.Tukey.whisker.bounds.by.trial.cluster <-
            compute.trial.specific.candidate.modes.from.subset.of.sampled.points( fit.rej.bin, hvtn705.parameters, hvtn705.target.stats.names, target.stats, target.stat.scale.units, units.quantile = base::min( 0.25, ( abc.keep.num.sims / num.sims ) ), max.dist.scaled = max.dist.scaled, pdfCluster.hmult = pdfCluster.hmult, Tukey.IQR.multiplier = Tukey.IQR.multiplier );

        for( hvtn705.cluster.j in 1:ncol( hvtn705.Tukey.whisker.bounds.by.trial.cluster[[ "low" ]] ) ) {
            if( be.verbose ) {
                cat( paste( "hvtn705.cluster.j = ", hvtn705.cluster.j ), fill = TRUE );
            }
            candidate.parameter.sets <-
                add.candidate.parameter.set( candidate.parameter.sets, hvtn705.cluster.j, hvtn705.Tukey.whisker.bounds.by.trial.cluster );
        } # End foreach hvtn705.cluster.j
    } # End foreach epsilon.bin
    return( candidate.parameter.sets );
} # get.candidate.parameter.sets (..)

bound.candidate.parameter.sets <- function ( candidate.parameter.sets, bounds.low = sapply( bounds, function ( .lower.and.higher ) { .lower.and.higher[ 1 ] } ), bounds.high = sapply( bounds, function ( .lower.and.higher ) { .lower.and.higher[ 2 ] } ) ) {

    ## Compute midpoints of candidates as medians of bounds; these
    ## will be the starting places for the optimizations below, but
    ## first we need to ensure everything is within the specified
    ## bounds.
    candidate.parameter.sets.midpoint <-
        candidate.parameter.sets[[ "low" ]] + ( candidate.parameter.sets[[ "high" ]] - candidate.parameter.sets[[ "low" ]] ) / 2;

    candidate.parameter.sets.low.bounded <- t( apply( candidate.parameter.sets[[ "low" ]], 1, function( .candidate.parameter.set.low ) { ifelse( .candidate.parameter.set.low < bounds.low, bounds.low, .candidate.parameter.set.low ) } ) );
    candidate.parameter.sets.high.bounded <- t( apply( candidate.parameter.sets[[ "high" ]], 1, function( .candidate.parameter.set.high ) { ifelse( .candidate.parameter.set.high > bounds.high, bounds.high, .candidate.parameter.set.high ) } ) );
    candidate.parameter.sets.midpoint.bounded <- t( apply( candidate.parameter.sets.midpoint, 1, function( .candidate.parameter.set.midpoint ) { ifelse( .candidate.parameter.set.midpoint > bounds.high, bounds.high, ifelse( .candidate.parameter.set.midpoint < bounds.low, bounds.low, .candidate.parameter.set.midpoint ) ) } ) );

    return( list( "low" = candidate.parameter.sets.low.bounded, "high" = candidate.parameter.sets.high.bounded, "midpoint" = candidate.parameter.sets.midpoint.bounded ) );
} # bound.candidate.parameter.sets

# uses params, so do within "with( runSim_hvtn705.create.params.list( reac ), .."
find.local.optimum <- function ( init, lower, upper, target.stat.stdevs ) {
    if( be.verbose ) {
        cat( "Starting midpoint: ", fill = FALSE );
        print( init );
        cat( "Optimization range: ", fill = FALSE );
        print( rbind( lower, upper ) );
    }
    .iterative.result <- optimize.iteratively( init, target.stats, target.stat.stdevs, lower = lower, upper = upper, current.value = NULL, reltol = optimize.iteratively.reltol, step.i = 1, max.steps = optimize.iteratively.max.steps, be.verbose = be.verbose );
    # current.parameters <- .iterative.result[ all.parameters ];
    # current.stats <- .iterative.result[[ setdiff( names( .iterative.result ), all.parameters, "dist" ) ]];
    # current.value <- .iterative.result[[ "dist" ]];
    return( .iterative.result );
} # find.local.optimum (..)

## This is the main function to run. It takes one argument, "reac" which is a named vector that must contain "numExecution" which is the number of initial random samples to start the process with. The result is some identified modes that are as close as we could get to the target stats, and their distances to the target stats. This is contained in the result, as a matrix called "sampled.modes", with results in columns, sorted by distance to target, low to high.
# MAGIC #s for the targets.
runSim_hvtn705 <- function( reac = c( "numExecution" = 10 ) ) { # Use numExecution >>1000 for best results.

    params.list <- runSim_hvtn705.create.params.list( reac );
    attach( params.list, name = "params.list" );

    #### Phase 1: Create starting points and optimization ranges: candidate.parameter.sets.bounded

    ## Phase 1a, slow step: draw numExecution draws from the independent priors.
    ## To restart from a saved one do this instead:

    ## TODO: REMOVE
    #load( file = "fit.rej.hvtn705.10k.threegroup.Rda" )
    .f.abc <- function( x ) {
        unlist( run.and.compute.run.stats.hvtn705( epsilon = x[ 1 ], hvtn705.log10lambda = x[ 2 ], hvtn705.log10riskmultiplier = x[ 3 ], hvtn705.highRiskProportion = x[ 4 ], hvtn705.lowRiskProportion = x[ 5 ] ) )
    }
    fit.rej <- draw.from.priors( .f.abc ); 
    #save( fit.rej, file = "fit.rej.hvtn705.10k.threegroup.Rda" )
    
    # Note that this is all we need to retain from the original fit.rej for phase 2:
    target.stat.stdevs <-
        ifelse( is.na( fit.rej$stats_normalization ), 1, fit.rej$stats_normalization ) * target.stat.scale.units;

    ## Phase 1b, slow step: use these draws to get candidate parameter sets for optimization.
    candidate.parameter.sets <-
        get.candidate.parameter.sets( fit.rej );
    if( be.verbose ) {
        cat( paste( "There are", nrow( candidate.parameter.sets[[ "high" ]] ), "candidate parameter sets." ), fill = TRUE );
    }

    ### PHASE 1c: Filter/merge the overlapping candidates so they are not redundant.
    ## TODO: [something; for now we just skip this step]

    ### PHASE 1d: Apply bounds to the computed parameter set search windows and compute initial values.
    candidate.parameter.sets.bounded <-
        bound.candidate.parameter.sets( candidate.parameter.sets );
    num.candidates <- nrow( candidate.parameter.sets.bounded[[ "midpoint" ]] );
    stopifnot( nrow( candidate.parameter.sets[[ "high" ]] ) == num.candidates );
    ## These are the only PHASE 1 outputs that are used below in subsequent phases:
    # target.stat.stdevs
    # candidate.parameter.sets.bounded

    ### PHASE 2:  from these candidate starting places constructed from epsion-bin-specific, trial-specific local optima, find modes in the 9-parameter space by conditional optimization.
    optimize.candidate <- function ( .candidate ) {
        if( be.verbose ) {
            cat( paste( "Optimizing candidate", .candidate ), fill = TRUE );
        }
        .lower.bounds <- candidate.parameter.sets.bounded[[ "low" ]][ .candidate, ];
        .upper.bounds <- candidate.parameter.sets.bounded[[ "high" ]][ .candidate, ];
        .midpoint <- candidate.parameter.sets.bounded[[ "midpoint" ]][ .candidate, ];
        return( find.local.optimum( init = .midpoint, lower = .lower.bounds, upper = .upper.bounds, target.stat.stdevs ) );
    } # optimize.candidate (..)
    optima.by.candidate <- sapply( 1:num.candidates, optimize.candidate );
    optima.by.candidate.sorted <-
        optima.by.candidate[ , order( as.numeric( optima.by.candidate[ "dist", ] ) ), drop = FALSE ];
        
    sim.results <- list( fit = fit.rej, priors = priors, bounds = bounds, target.stats = target.stats, fn = .f.abc, sampled.modes = optima.by.candidate.sorted );

    detach( name = "params.list" );

    return( sim.result );
} # runSim_hvtn705 (..)

the.seed <- 98103;
# To test replicability of the identified modes, uncomment this:
# set.seed( the.seed ); the.seed <- floor( runif( 1, max = 1E5 ) );
# num.sims <- 1000; # Fast for debugging.
num.sims <- 10000; # For reals.

set.seed( the.seed );

# This runs it. I've commented it out so you can "source" this file safely.
# .sim <- runSim_hvtn705( reac = c( "numExecution" = num.sims ) );

###################################################################################################
######## Documentation on the provenance of the target statistics:
###################################################################################################
### HVTN 705 results: see HVTN705_cumulative_infections.pdf for the figure showing the placebo incidence (and other data), and below for the VE point estimate.
### From https://www.prnewswire.com/news-releases/johnson--johnson-and-global-partners-announce-results-from-phase-2b-imbokodo-hiv-vaccine-clinical-trial-in-young-women-in-sub-saharan-africa-301365918.html
## What the Imbokodo Data Tell Us
# The Imbokodo vaccine regimen was administered to participants through four vaccination visits over one year. The primary analysis was conducted 24 months after participants received their first vaccinations. The study's primary endpoint was based on the difference in number of new HIV infections between the placebo and vaccine groups from month seven (one month after the third vaccination timepoint) through month 24. These data found that through 24 months of follow up, 63 of 1,109 participants who received placebo compared to 51 of 1,079 participants who received active vaccine acquired HIV. This analysis demonstrated a vaccine efficacy point estimate of 25.2% (95% confidence interval of -10.5% to 49.3%). The vaccine regimen did not cause harm and was generally well-tolerated.
#####################
