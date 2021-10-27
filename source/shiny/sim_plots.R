#------------------------------------------------------------------------------
# for creating cumulative infections plot
#------------------------------------------------------------------------------
createCumulativeInfectionsPlot <- function(output, reac) {
  output$CumulativeInfections <- renderPlot({
    
    mod <- runSim(reac)
    
    plot(mod, y = c("Ip", "total.Iph.Ipm.Ipl"), 
         alpha = 0.8, 
         main = "Cumulative infections",
         legend = FALSE,
         ylab = "infected",
         xlab = "days",
         col = c("blue", "red"))
    legend("bottomright", legend = c("homogeneous risk", "heterogeneous risk"), col = c("blue", "red"), lwd = 2, cex = 0.9, bty = "n")
  })
}

#------------------------------------------------------------------------------
# for creating placebo incidence plot
#------------------------------------------------------------------------------
createPlaceboRiskPlot <- function(output, reac) {
  output$PlaceboRiskPlot <- renderPlot({
    
    mod <- runSim(reac)
    
    plot(mod, y=c("rate.Placebo", "rate.Placebo.het"),
         alpha = 0.8,
         ylim = c(0, 6.0),
         main = "Instantaneous incidence in the placebo arm",
         xlab = "days",
         ylab = "infections per 100 person yrs",
         legend = FALSE,
         col = c("blue", "red"))
    legend("bottomright", legend = c("homogen. risk", "heterogen. risk"), col = c("blue", "red"), lwd = 2, cex = 0.9, bty = "n")
    
    par(mfrow=c(1,1))
  })
}

#------------------------------------------------------------------------------
# for creating placebo vaccine risk plot
#------------------------------------------------------------------------------
createPlaceboVaccineRiskPlot <- function(output, reac) {
  output$PlaceboVaccineRiskPlot <- renderPlot({
    mod <- runSim(reac)
    plot(mod, y=c("rate.Placebo", "rate.Vaccine", "rate.Placebo.het", "rate.Vaccine.het"),
         alpha = 0.8,
         ylim = c(0, 6.0),
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

#------------------------------------------------------------------------------
# for creating vaccine efficacy plot
#------------------------------------------------------------------------------
createVEPlot <- function(output, reac) {
  output$VEPlot <- renderPlot({
    mod <- runSim(reac)
    plot(mod, y=c("VE1.inst", "VE2.inst"),
       alpha = 0.8,
       main = "Vaccine efficacy",
       legend = FALSE, 
       xlab = "days",
       ylab = "Vaccine efficacy",
       col = c("blue", "red"))
    legend("topright", legend = c("VE, homogeneous risk", "VE, heterogeneous risk"), col = c("blue", "red"), lwd = 2, cex = 0.9, bty = "n")
  })
  
}

