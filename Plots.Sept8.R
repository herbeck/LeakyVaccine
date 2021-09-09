#load("~/Dropbox/R01.MIDAS.modeling/LeakyVaccine/fit.rej.unif.10k.Rda")
#load("~/Dropbox/R01.MIDAS.modeling/LeakyVaccine/sim.results.unif.10k.Rda")

fit <- fit.rej
bounds <- sim.results$bounds
RV144.VE.target <- sim.results$target.stats$rv144.VE.target
HVTN702.VE.target <- sim.results$target.stats$hvtn702.VE.target

# These are the original bounds, separated for use in optim:
lower.bounds <- sapply( bounds, function ( .bounds.for.x ) { .bounds.for.x[ 1 ] } );
upper.bounds <- sapply( bounds, function ( .bounds.for.x ) { .bounds.for.x[ 2 ] } );

plot(density(fit$param[, 1], from = lower.bounds[1],  to = upper.bounds[1]),
     main = "epsilon", 
     xlim = c(lower.bounds[1], upper.bounds[1]),
     #ylim = c(0, 10),
     col=1)
#lines(density(fit$param[, 1], from = lower.bounds[1],  to = upper.bounds[1]), col = 2)
abline(v = RV144.VE.target/100, lty = 2, col = 2)
abline(v = HVTN702.VE.target/100, lty = 2, col = 4)
legend("bottomright", legend = c("RV144 target VE", "HVTN702 target VE", "Posterior"),
       lty = c(2, 2, 1), col = c(2, 4, 1), lwd = 1, cex = 0.6)

plot(density(fit$param[, 2], from = lower.bounds[2],  to = upper.bounds[2]),
     main = "log10lambda", 
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
names( .df ) <- c( "epsilon", "log10lambda", "risk" )
pdf( "log10lambda_risk.pdf" ); ggplot( .df, aes( x=log10lambda,y=risk ) ) + geom_point() + stat_density2d_filled(); dev.off()
pdf( "epsilon_risk.pdf" ); ggplot( .df, aes( x=epsilon,y=risk ) ) + geom_point() + stat_density2d_filled(); dev.off()
pdf( "epsilon_log10lambda.pdf" ); ggplot( .df, aes( x=epsilon,y=log10lambda ) ) + geom_point() + stat_density2d_filled();

dev.off()

