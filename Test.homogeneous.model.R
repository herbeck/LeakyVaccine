


library(deSolve)
library(tidyverse)
library(EpiModel)
library(survival)
library(EasyABC)


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
  
  param <- param.dcm(lambda = lambda, epsilon = epsilon, inc = inc, prop_high = prop_high)
  mod <- dcm(param, init, control)
  mod
  
  mod <- mod.manipulate(mod)
  
  output$cum_efficacy[output$prop_high == prop_high] <- mod$epi$VE2.cumul[, 1]
  output$inst_efficacy[output$prop_high == prop_high] <- mod$epi$VE2.inst[, 1]
  
}
