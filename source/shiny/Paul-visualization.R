

source("model/Paul-lib.R") 

createVisualization <- function(output, reac) {
  
  fit_model <- reactive({
    fit <-runSim_Paul(reac);
  });
  
  bounds <- list(#c(0.003, 0.008),      # beta
    #c(60/365, 120/365))    # c
    epsilon = c(1E-10, 1-1E-10),
    #               lambda = c(0.000005, 0.0001),  # lambda
    lambda = c(0.000005, 0.00001),  # lambda
    risk = c(1, 20))            # risk multiplier
  
  
  time <- 3*365; # End of the trial.
  
  ## Note here actually want to target the average incidence over time, rather than the incidence at the end of the trial, since with multiple risk groups there is waning expected. I have changed this in the ABC function f, below.
  placebo.incidence.target <- rep( 3, length( time ) )    # flat incidence of 3.5% per 100 person years
  
  ### Paul notes this is not going to work, if efficacy is waning you can't match a constant efficacy over time. We should match just one time, probably duration of 702 trial.
  VE.target <- rep(0.3, length( time ))
  target.stats <- data.frame(time, VE.target, placebo.incidence.target )
  
  lower.bounds <- sapply( bounds, function ( .bounds.for.x ) { .bounds.for.x[ 1 ] } );
  upper.bounds <- sapply( bounds, function ( .bounds.for.x ) { .bounds.for.x[ 2 ] } );
  
  #------------------------------------------------------------------------------
  # for creating Lambda risk plot
  #------------------------------------------------------------------------------
  output$plotTestLambdaRisk <- renderPlot({
    
    fit <- fit_model()
    
    plot(density(fit$param[, 1], from = lower.bounds[1],  to = upper.bounds[1]),
         main = "epsilon",
         xlim = c(lower.bounds[1], upper.bounds[1]),
         #ylim = c(0, 10),
         col=2)
    lines(density(fit$param[, 1], from = lower.bounds[1],  to = upper.bounds[1]), col = 2)
    abline(v = VE.target, lty = 2, col = 1)
    legend("topright", legend = c("Truth", "Posterior"),
           lty = c(1, 2), col = 1:2, lwd = 2)
    
    plot(density(fit$param[, 2], from = lower.bounds[2],  to = upper.bounds[2]),
         main = "lambda",
         xlim = c(lower.bounds[2], upper.bounds[2]),
         col=2)
    lines(density(fit$param[, 2], from = lower.bounds[2],  to = upper.bounds[2]), col = 2)
    abline(v = VE.target, lty = 2, col = 1)
    legend("topright", legend = c("Truth", "Posterior"),
           lty = c(1, 2), col = 1:2, lwd = 2)
    
    plot(density(fit$param[, 3], from = lower.bounds[3],  to = upper.bounds[3]),
         main = "risk",
         xlim = c(lower.bounds[3], upper.bounds[3]),
         col=2)
    lines(density(fit$param[, 3], from = lower.bounds[3],  to = upper.bounds[3]), col = 2)
    abline(v = VE.target, lty = 2, col = 1)
    legend("topright", legend = c("Truth", "Posterior"),
           lty = c(1, 2), col = 1:2, lwd = 2)
    
    .df <- as.data.frame( fit$param )
    names( .df ) <- c( "epsilon", "lambda", "risk" )
    ggplot( .df, aes( x=lambda,y=risk ) ) + geom_point() + stat_density2d_filled(); #dev.off()
    
  });
  
  #------------------------------------------------------------------------------
  # for creating Epsilon risk plot
  #------------------------------------------------------------------------------
  output$plotTestEspilonRisk <- renderPlot({
    
    fit <- fit_model()
    
    plot(density(fit$param[, 1], from = lower.bounds[1],  to = upper.bounds[1]),
         main = "epsilon",
         xlim = c(lower.bounds[1], upper.bounds[1]),
         #ylim = c(0, 10),
         col=2)
    lines(density(fit$param[, 1], from = lower.bounds[1],  to = upper.bounds[1]), col = 2)
    abline(v = VE.target, lty = 2, col = 1)
    legend("topright", legend = c("Truth", "Posterior"),
           lty = c(1, 2), col = 1:2, lwd = 2)
    
    plot(density(fit$param[, 2], from = lower.bounds[2],  to = upper.bounds[2]),
         main = "lambda",
         xlim = c(lower.bounds[2], upper.bounds[2]),
         col=2)
    lines(density(fit$param[, 2], from = lower.bounds[2],  to = upper.bounds[2]), col = 2)
    abline(v = VE.target, lty = 2, col = 1)
    legend("topright", legend = c("Truth", "Posterior"),
           lty = c(1, 2), col = 1:2, lwd = 2)
    
    plot(density(fit$param[, 3], from = lower.bounds[3],  to = upper.bounds[3]),
         main = "risk",
         xlim = c(lower.bounds[3], upper.bounds[3]),
         col=2)
    lines(density(fit$param[, 3], from = lower.bounds[3],  to = upper.bounds[3]), col = 2)
    abline(v = VE.target, lty = 2, col = 1)
    legend("topright", legend = c("Truth", "Posterior"),
           lty = c(1, 2), col = 1:2, lwd = 2)
    
    .df <- as.data.frame( fit$param )
    names( .df ) <- c( "epsilon", "lambda", "risk" )
    ggplot( .df, aes( x=epsilon,y=risk ) ) + geom_point() + stat_density2d_filled(); #dev.off()
    
  });
  
  
  #------------------------------------------------------------------------------
  # for creating Epsilon risk plot
  #------------------------------------------------------------------------------
  output$plotTestEpsilonLambda <- renderPlot({
    
    fit <- fit_model()
    
    plot(density(fit$param[, 1], from = lower.bounds[1],  to = upper.bounds[1]),
         main = "epsilon",
         xlim = c(lower.bounds[1], upper.bounds[1]),
         #ylim = c(0, 10),
         col=2)
    lines(density(fit$param[, 1], from = lower.bounds[1],  to = upper.bounds[1]), col = 2)
    abline(v = VE.target, lty = 2, col = 1)
    legend("topright", legend = c("Truth", "Posterior"),
           lty = c(1, 2), col = 1:2, lwd = 2)
    
    plot(density(fit$param[, 2], from = lower.bounds[2],  to = upper.bounds[2]),
         main = "lambda",
         xlim = c(lower.bounds[2], upper.bounds[2]),
         col=2)
    lines(density(fit$param[, 2], from = lower.bounds[2],  to = upper.bounds[2]), col = 2)
    abline(v = VE.target, lty = 2, col = 1)
    legend("topright", legend = c("Truth", "Posterior"),
           lty = c(1, 2), col = 1:2, lwd = 2)
    
    plot(density(fit$param[, 3], from = lower.bounds[3],  to = upper.bounds[3]),
         main = "risk",
         xlim = c(lower.bounds[3], upper.bounds[3]),
         col=2)
    lines(density(fit$param[, 3], from = lower.bounds[3],  to = upper.bounds[3]), col = 2)
    abline(v = VE.target, lty = 2, col = 1)
    legend("topright", legend = c("Truth", "Posterior"),
           lty = c(1, 2), col = 1:2, lwd = 2)
    
    .df <- as.data.frame( fit$param )
    names( .df ) <- c( "epsilon", "lambda", "risk" )
    
    ggplot( .df, aes( x=epsilon,y=lambda ) ) + geom_point() + stat_density2d_filled(); #dev.off()
    
  });
  
  
  
  ## TODO: 3d plots
}