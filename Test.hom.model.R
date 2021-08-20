


library(deSolve)
library(tidyverse)
library(EpiModel)
library(survival)
library(EasyABC)

## Test effect of different incidence levels on vaccine efficacy
## in a homogeneous risk population

#incs <- seq(0.003, 0.03, by = 0.005)

betas <- seq(0.002, 0.02, by = 0.0005)   #transmission rate (per contact)
c <- 90/365    #contact rate (contacts per day)
prev <- 0.10    #needs some consideration
lambda <- beta*c*prev
epsilon <- 0.30  #per contact vaccine efficacy
risk <- 10.0   #risk multiplier
n <- 5000 #sample size
#prop_high <- 0.00  #proportion high exposure
#prop_low <- 0.25  #proportion low exposure

output <- data.frame(beta = rep(betas, each = nsteps), step = 1:nsteps, cum_efficacy = NA, inst_efficacy = NA)

for(beta in betas) {
  
  ## Calculate risk multiplier
  #risk <- uniroot(calc_inc, prop_high = prop_high, lambda = lambda, inc = inc, c(0, 1000), tol = 0.0001)$root
  print(prop_high)
  print(risk)
  
  init <- init.dcm(Sp = n, Ip = 0,  #placebo, homogeneous risk population
                   Sv = n, Iv = 0,  #vaccine, homogeneous risk population
                   #Sph = n*prop_high, Iph = 0,    #placebo, high risk
                   #Spm = n-(n*prop_high)-(n*prop_low), Ipm = 0,    #placebo, medium risk
                   #Spl = n*prop_low, Ipl = 0,    #placebo, low risk
                   #Svh = n*prop_high, Ivh = 0,    #vaccine, high risk
                   #Svm = n-(n*prop_high)-(n*prop_low), Ivm = 0,    #vaccine, medium risk
                   #Svl = n*prop_low, Ivl = 0,    #vaccine, low risk
                   SIp.flow = 0, SIv.flow = 0) 
                   #SIph.flow = 0, SIpm.flow = 0, SIpl.flow = 0,
                   #SIvh.flow = 0, SIvm.flow = 0, SIvl.flow = 0)
  
  control <- control.dcm(nsteps = nsteps, new.mod = si_ode)
  param <- param.dcm(lambda = lambda, epsilon = epsilon, inc = inc, prop_high = prop_high)
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

