#load("~/Dropbox/R01.MIDAS.modeling/LeakyVaccine/fit.rej.unif.10k.Rda")
#load("~/Dropbox/R01.MIDAS.modeling/LeakyVaccine/sim.results.unif.10k.Rda")

fit <- fit.rej
bounds <- sim.results$bounds
target.stats <- sim.results$target.stats
RV144.VE.target <- sim.results$target.stats$rv144.VE.target
HVTN702.VE.target <- sim.results$target.stats$hvtn702.VE.target

placebo.incidence.target.scale.units <- 1
VE.target.scale.units <- 0.1

target.stat.scale.units <- c( "rv144.VE.target" = VE.target.scale.units, 
                              "rv144.placebo.incidence.target" = placebo.incidence.target.scale.units, 
                              "hvtn702.VE.target" = VE.target.scale.units, 
                              "hvtn702.placebo.incidence.target" = placebo.incidence.target.scale.units )


########

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

########

# These are the original bounds, separated for use in optim:
lower.bounds <- sapply( bounds, function ( .bounds.for.x ) { .bounds.for.x[ 1 ] } );
upper.bounds <- sapply( bounds, function ( .bounds.for.x ) { .bounds.for.x[ 2 ] } );

# epsilon (per-exposure VE; true VE) (used for both RV144 and HVTN702 sims)

plot(density(fit.filtered$param[, 1], from = lower.bounds[1],  to = upper.bounds[1]),
     main = "epsilon", 
     xlim = c(lower.bounds[1], upper.bounds[1]),
     #ylim = c(0, 10),
     col=1)
abline(v = RV144.VE.target/100, lty = 2, col = 2)
abline(v = HVTN702.VE.target/100, lty = 2, col = 4)
lines(density(fit$param[, 1], from = lower.bounds[1],  to = upper.bounds[1]), col = 1, lty=3)
legend("bottomright", legend = c("RV144 target VE", "HVTN702 target VE", "Posterior"),
       lty = c(2, 2, 1), col = c(2, 4, 1), lwd = 1, cex = 0.6)

#### RV144 plots and 702 plots, overlaid

# lambda (force of infection; beta*prev*c)
plot(density(fit.filtered$param[, 6], from = lower.bounds[6],  to = upper.bounds[6]),
     main = "log10lambda", 
     xlim = c(lower.bounds[6], upper.bounds[6]),
     #ylim = 
     col=2)
lines(density(fit.filtered$param[, 2], from = lower.bounds[2],  to = upper.bounds[2]), col = 1)
lines(density(fit$param[, 2], from = lower.bounds[2],  to = upper.bounds[2]), col = 1, lty=3)
lines(density(fit$param[, 6], from = lower.bounds[6],  to = upper.bounds[6]), col = 1, lty=3)
legend("topleft", legend = c("RV144", "HVTN702"),
       lty = 1, col = 1:2, lwd = 1, cex = 0.75)

# high risk multiplier
plot(density(fit.filtered$param[, 7], from = lower.bounds[7],  to = upper.bounds[7]),
     main = "risk multiplier", 
     xlim = c(lower.bounds[7], upper.bounds[7]),
     col=2)
lines(density(fit.filtered$param[, 3], from = lower.bounds[3],  to = upper.bounds[3]), col = 1)
legend("bottomright", legend = c("RV144", "HVTN702"),
       lty = 1, col = 1:2, lwd = 1, cex = 0.75)

# high risk proportion
plot(density(fit.filtered$param[, 8], from = lower.bounds[8],  to = upper.bounds[8]),
     main = "high risk proportion", 
     xlim = c(lower.bounds[8], upper.bounds[8]),
     col=2)
lines(density(fit.filtered$param[, 4], from = lower.bounds[4],  to = upper.bounds[4]), col = 1)
legend("bottomright", legend = c("RV144", "HVTN702"),
       lty = 1, col = 1:2, lwd = 1, cex = 0.75)


########################
# contour plots        #
########################

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

## This is just printing the rv144 contours, note.
.df <- as.data.frame( fit.filtered$param[ , c( "epsilon", "rv144.log10lambda", "rv144.high.risk.multiplier" ) ] );
names( .df ) <- c( "epsilon", "log10lambda", "risk" )
pdf( "rv144log10lambda_rv144risk.pdf" ); ggplot( .df, aes( x=log10lambda,y=risk ) ) + geom_point() + stat_density2d_filled(); dev.off()
pdf( "epsilon_rv144risk.pdf" ); ggplot( .df, aes( x=epsilon,y=risk ) ) + geom_point() + stat_density2d_filled(); dev.off()
pdf( "epsilon_rv144log10lambda.pdf" ); ggplot( .df, aes( x=epsilon,y=log10lambda ) ) + geom_point() + stat_density2d_filled(); dev.off()

.df <- as.data.frame( fit.filtered$param[ , c( "epsilon", "hvtn702.log10lambda", "hvtn702.high.risk.multiplier" ) ] );
names( .df ) <- c( "epsilon", "log10lambda", "risk" )
pdf( "hvtn702log10lambda_rv144risk.pdf" ); ggplot( .df, aes( x=log10lambda,y=risk ) ) + geom_point() + stat_density2d_filled(); dev.off()
pdf( "epsilon_hvtn702risk.pdf" ); ggplot( .df, aes( x=epsilon,y=risk ) ) + geom_point() + stat_density2d_filled(); dev.off()
pdf( "epsilon_hvtn702log10lambda.pdf" ); ggplot( .df, aes( x=epsilon,y=log10lambda ) ) + geom_point() + stat_density2d_filled(); dev.off()


