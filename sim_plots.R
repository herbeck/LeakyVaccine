#------------------------------------------------------------------------------
# for creating cumulative infection plot
#------------------------------------------------------------------------------
createCInfectionPlot <- function(output, reac) {
  output$plotCInfect <- renderPlot({
    
    mod <- runSim(reac)
    
    plot(mod, y = c("Ip", "total.Iph.Ipm.Ipl"), 
         alpha = 0.8, 
         main = "Cumulative infections",
         legend = FALSE,
         ylab = "infected",
         xlab = "days",
         col = c("blue", "red"))
    
    legend("bottomright", legend = c("homogeneous risk", "heterogeneous risk"), 
           col = c("blue", "red"), lwd = 2, cex = 0.9, bty = "n")
  })
}

#------------------------------------------------------------------------------
# for creating placebo risk plot
#------------------------------------------------------------------------------
createPlaceboRiskPlot <- function(output, reac) {
  output$plot2 <- renderPlot({
    
    mod <- runSim(reac)
    
    plot(mod, y=c("rate.Placebo", "rate.Placebo.het"),
         alpha = 0.8,
         ylim = c(0, 4.5),
         main = "Hazard",
         xlab = "days",
         ylab = "infections per 100 person yrs",
         legend = FALSE,
         col = c("blue", "red"))
    legend("bottomright", legend = c("homogen. risk", "heterogen. risk"), col = c("blue", "red"), lwd = 2, cex = 0.9, bty = "n")
    
    par(mfrow=c(1,1))
  })
}

#------------------------------------------------------------------------------
# for creating cumulative incidence plot
#------------------------------------------------------------------------------
createCIncidencePlot <- function(output, reac) {
  output$plot3 <- renderPlot({
    
    mod <- runSim(reac)
    
    #mod <- mod.manipulate(mod)
    
    plot(mod, y=c("cumul.rate.Placebo", "cumul.rate.Vaccine"),
         alpha = 0.8,
         ylim = c(0, 4.5),
         main = "Cumulative incidence",
         xlab = "days",
         ylab = "infections per 100 person yrs",
         legend = FALSE,
         col = c("blue", "green"))
    legend("bottomright", legend = c("placebo", "vaccine"), col = c("blue", "green"), lwd = 2, cex = 0.9, bty = "n")
  })
}

#------------------------------------------------------------------------------
# for creating placebo vs vaccine plot
#------------------------------------------------------------------------------
createPlaceboVaccinePlot <- function(output, reac) {
  output$plot4 <- renderPlot({
    mod <- runSim(reac)
    
    plot(mod, y=c("rate.Placebo", "rate.Vaccine"),
         alpha = 0.8,
         ylim = c(0, 4.5),
         main = "Hazard",
         xlab = "days",
         ylab = "infections per 100 person yrs",
         legend = FALSE,
         col = c("blue", "green"))
    legend("bottomright", legend = c("placebo", "vaccine"), col = c("blue", "green"), lwd = 2, cex = 0.9, bty = "n")
  })
}

#------------------------------------------------------------------------------
# for creating placebo vaccine risk plot
#------------------------------------------------------------------------------
createPlaceboVaccineRiskPlot <- function(output, reac) {
  output$plot5 <- renderPlot({
    mod <- runSim(reac)
    plot(mod, y=c("rate.Placebo", "rate.Vaccine", "rate.Placebo.het", "rate.Vaccine.het"),
         alpha = 0.8,
         ylim = c(0, 4.5),
         main = "Hazard",
         xlab = "days",
         ylab = "infections per 100 person yrs",
         legend = FALSE,
         col = c("blue", "green", "red", "orange"))
    legend("bottomright", legend = c("placebo, homogeneous risk", "vaccine, homogeneous risk", 
                                     "placebo, heterogeneous risk", "vaccine, heterogeneous risk"), 
           col = c("blue", "green", "red", "orange"), lwd = 2, cex = 0.9, bty = "n")
  })
}
