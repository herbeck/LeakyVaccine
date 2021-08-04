

source("Paul-lib.R") 

createVisualization <- function(output, reac) {
  
  fit_model <- reactive({
    fit <-runSim_Paul(reac);
  });
  
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