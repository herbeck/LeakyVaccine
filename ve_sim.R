## Allen Roberts
## April 2021
## Investigate waning estimated vaccine efficacy under different levels of heterogeneity and initial conditions

#rm(list = ls())


library(deSolve)
library(tidyverse)
library(EpiModel)
library(survival)
library(EasyABC)
library(ggplot2)
library(viridis)

## Helper functions
source("ve_sim_fns.R")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Run simulation with different proportions in high risk group
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
runSimByPropHigh <- function(param) {
  
  beta <- param$beta        # transmission rate (per contact)
  c <- param$contactRate    # contact rate (contacts per day). Sets underlying risk in low risk group.
  prev <- param$prev        # needs some more consideration
  lambda <- beta*c*prev
  1-exp(-lambda*365)        # Annual risk 
  epsilon <- param$epsilon  # per contact vaccine efficacy
  n <- param$sampleSize     # Population (sample) size
  inc <- param$inc          # Overall estimated annual risk to calibrate to
  nsteps <- param$nsteps  
  
  ## Test different proportions in high risk group
  prop_highs <- seq(0.05, 0.5, by = 0.05)
  output <- data.frame(prop_high = rep(prop_highs, each = nsteps), step = 1:nsteps, cum_efficacy = NA, inst_efficacy = NA)
  
  for(prop_high in prop_highs) {
    
    ## Calculate risk multiplier
    risk <- uniroot(calc_inc, prop_high = prop_high, lambda = lambda, inc = inc, c(0, 1000), tol = 0.0001)$root
    print(prop_high)
    print(risk)
    
    init <- init.dcm(Sp = n, Ip = 0,
                     Sv = n, Iv = 0,
                     Sph = n*prop_high, Iph = 0,    #placebo, high risk
                     Spl = n*(1-prop_high), Ipl = 0,    #placebo, low risk
                     Svh = n*prop_high, Ivh = 0,    #vaccine
                     Svl = n*(1-prop_high), Ivl = 0,    #vaccine
                     SIp.flow = 0, SIv.flow = 0, 
                     SIph.flow = 0, SIpl.flow = 0,
                     SIvh.flow = 0, SIvl.flow = 0)
    control <- control.dcm(nsteps = nsteps, new.mod = si_ode)
    

    param <- param.dcm(lambda = lambda, epsilon = epsilon, inc = inc, prop_high = prop_high, risk = risk, n=n)
    print("before mod")
    mod <- dcm(param, init, control)
    print ("after mod")
    mod
    
    mod <- mod.manipulate(mod)
    
    output$cum_efficacy[output$prop_high == prop_high] <- mod$epi$VE2.cumul[, 1]
    output$inst_efficacy[output$prop_high == prop_high] <- mod$epi$VE2.inst[, 1]
    
  }
  
  #browser()
  
  out <- output %>%
    pivot_longer(cols = c("cum_efficacy", "inst_efficacy"), names_to = "metric")
    
  ve_by_prop_high <- ggplot(data = out, aes(x = step, y = value, group = prop_high)) +
    geom_line(aes(color = prop_high)) +
    geom_abline(intercept = epsilon, slope = 0, linetype = "dashed") +
    scale_color_viridis(name = "Proportion high risk") +
    labs(x = "Time (days)", y = "Estimated vaccine efficacy") +
    scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 1, by = 0.1)) +
    facet_wrap(~metric) +
    ggtitle(paste("Incidence = ", inc, "VE =", epsilon)) +
    theme_classic()
  
  return (ve_by_prop_high)
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Run simulation with different incidence levels
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
runSimByInc <- function(param) {
  
  beta <- param$beta        # transmission rate (per contact)
  c <- param$contactRate    # contact rate (contacts per day). Sets underlying risk in low risk group.
  prev <- param$prev        # needs some more consideration
  lambda <- beta*c*prev
  1-exp(-lambda*365)        # Annual risk 
  epsilon <- param$epsilon  # per contact vaccine efficacy
  n <- param$sampleSize     # Sample size
  inc <- param$inc          # Overall estimated annual risk to calibrate to
  nsteps <- param$nsteps  
  
  ## Test different incidence levels
  incs <- seq(0.01, 0.06, by = 0.005)
  prop_high <- 0.1
  
  output <- data.frame(inc = rep(incs, each = nsteps), step = 1:nsteps, cum_efficacy = NA, inst_efficacy = NA)
  
  for(inc in incs) {
    
    ## Calculate risk multiplier
    risk <- uniroot(calc_inc, prop_high = prop_high, lambda = lambda, inc = inc, c(0, 1000), tol = 0.0001)$root
    print(inc)
    print(risk)
    
    init <- init.dcm(Sp = n, Ip = 0,
                     Sv = n, Iv = 0,
                     Sph = n*prop_high, Iph = 0,    #placebo, high risk
                     Spl = n*(1-prop_high), Ipl = 0,    #placebo, low risk
                     Svh = n*prop_high, Ivh = 0,    #vaccine
                     Svl = n*(1-prop_high), Ivl = 0,    #vaccine
                     SIp.flow = 0, SIv.flow = 0, 
                     SIph.flow = 0, SIpl.flow = 0,
                     SIvh.flow = 0, SIvl.flow = 0)
    control <- control.dcm(nsteps = nsteps, new.mod = si_ode)
    
    param <- param.dcm(lambda = lambda, epsilon = epsilon, inc = inc, prop_high = prop_high, risk=risk, n=n)
    mod <- dcm(param, init, control)
    mod
    
    mod <- mod.manipulate(mod)
    
    
    
    output$cum_efficacy[output$inc == inc] <- mod$epi$VE2.cumul[, 1]
    output$inst_efficacy[output$inc == inc] <- mod$epi$VE2.inst[, 1]
    
  }
  
  out <- output %>%
    pivot_longer(cols = c("cum_efficacy", "inst_efficacy"), names_to = "metric")
  
  ve_by_inc <- ggplot(data = out, aes(x = step, y = value, group = inc)) +
    geom_line(aes(color = inc)) +
    geom_abline(intercept = epsilon, slope = 0, linetype = "dashed") +
    scale_color_viridis(name = "Incidence") +
    labs(x = "Time (days)", y = "Estimated vaccine efficacy") +
    scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 1, by = 0.1)) +
    facet_wrap(~metric) +
    ggtitle(paste("Prop high = ", prop_high, "VE =", epsilon)) +
    theme_classic()
  
  return (ve_by_inc)
  
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Run simulation with different vaccine efficacies
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
runSimByEpsilon <- function(param) {
  beta <- param$beta        # transmission rate (per contact)
  c <- param$contactRate    # contact rate (contacts per day). Sets underlying risk in low risk group.
  prev <- param$prev        # needs some more consideration
  lambda <- beta*c*prev
  1-exp(-lambda*365)        # Annual risk 
  epsilon <- param$epsilon  # per contact vaccine efficacy
  n <- param$sampleSize     # Sample size
  inc <- param$inc          # Overall estimated annual risk to calibrate to
  nsteps <- param$nsteps  
  
  ## Test different vaccine efficacies
  epsilons <- seq(0.1, 0.9, by = 0.1)
  prop_high <- 0.1
  inc <- 0.03
  
  output <- data.frame(epsilon = rep(epsilons, each = nsteps), step = 1:nsteps, cum_efficacy = NA, inst_efficacy = NA)
  
  for(epsilon in epsilons) {
    
    ## Calculate risk multiplier
    risk <- uniroot(calc_inc, prop_high = prop_high, lambda = lambda, inc = inc, c(0, 1000), tol = 0.0001)$root
    print(inc)
    print(risk)
    
    init <- init.dcm(Sp = n, Ip = 0,
                     Sv = n, Iv = 0,
                     Sph = n*prop_high, Iph = 0,    #placebo, high risk
                     Spl = n*(1-prop_high), Ipl = 0,    #placebo, low risk
                     Svh = n*prop_high, Ivh = 0,    #vaccine
                     Svl = n*(1-prop_high), Ivl = 0,    #vaccine
                     SIp.flow = 0, SIv.flow = 0, 
                     SIph.flow = 0, SIpl.flow = 0,
                     SIvh.flow = 0, SIvl.flow = 0)
    control <- control.dcm(nsteps = nsteps, new.mod = si_ode)
    
    param <- param.dcm(lambda = lambda, epsilon = epsilon, inc = inc, prop_high = prop_high, risk=risk, n=n)
    mod <- dcm(param, init, control)
    mod
    
    mod <- mod.manipulate(mod)
    
    output$cum_efficacy[output$epsilon == epsilon] <- mod$epi$VE2.cumul[, 1]
    output$inst_efficacy[output$epsilon == epsilon] <- mod$epi$VE2.inst[, 1]
    
  }
  
  out <- output %>%
    pivot_longer(cols = c("cum_efficacy", "inst_efficacy"), names_to = "metric")
  
  ve_by_epsilon <- ggplot(data = out, aes(x = step, y = value, group = epsilon)) +
    geom_line(aes(color = epsilon)) +
    geom_abline(aes(intercept = epsilon, slope = 0, colour = epsilon), linetype = "dashed") +
    scale_color_viridis(name = "Epsilon") +
    labs(x = "Time (days)", y = "Estimated vaccine efficacy") +
   # scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 1, by = 0.1)) +
    facet_wrap(~metric) +
    ggtitle(paste("Prop high = ", prop_high, "Incidence =", inc)) +
    theme_classic()
  
  ## Plot
  #browser()
  #pdf(file = "ve_sim_plots.pdf", height = 5, width = 5)
  # ve_by_prop_high
  # ve_by_inc
  # ve_by_epsilon
  # 
  return (ve_by_epsilon)
  
  dev.off()
}

# par(mar = c(3,3,2,1), mgp = c(2,1,0))
# plot(mod, y = c("Iv", "Ip", "total.Ivh.Ivl", "total.Iph.Ipl"), 
#      alpha = 0.8, 
#      main = "Cumulative infections",
#      legend = "full")
# plot(mod, y = c("SIv.flow", "SIp.flow", "SIvh.flow", "SIvl.flow", "SIph.flow", "SIpl.flow"), 
#      alpha = 0.8, 
#      main = "Daily infections, w/ all risk groups shown",
#      legend = "full")
# plot(mod, y = c("SIv.flow", "SIp.flow", "total.SIvh.SIvl.flow", "total.SIph.SIpl.flow"),
#      alpha = 0.8, 
#      main = "Daily infections, total mixed risk group",
#      legend = "full")
# plot(mod, y=c("rate.Vaccine", "rate.Placebo", "rate.Vaccine.het", "rate.Placebo.het"),
#      alpha = 0.8,
#      #ylim = c(0, 0.1),
#      main = "Instantaneous incidence rate",
#      legend = "full")
# plot(mod, y=c("cumul.rate.Vaccine", "cumul.rate.Placebo", "cumul.rate.Vaccine.het", "cumul.rate.Placebo.het"),
#      alpha = 0.8,
#      #ylim = c(0, 0.1),
#      main = "Cumulative incidence",
#      legend = "full")
# plot(mod, y=c("VE1.inst", "VE2.inst", "VE1.cumul", "VE2.cumul"),
#      alpha = 0.8,
#      main = "Vaccine efficacy",
#      legend = FALSE, 
#      col = 1:4)
# legend("topright", legend = c("Instantanteous VE, homogeneous risk", "Inst VE, heterogeneous risk", "Cumulative VE, homogeneous risk", "Cumul VE, heterogeneous risk"),
#        col = 1:4, lwd = 2, cex = 0.9, bty = "n")

