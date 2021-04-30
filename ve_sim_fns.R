library(deSolve)
library(tidyverse)
library(EpiModel)
library(survival)
library(EasyABC)

## Calibrate risk multiplier to overall lambda
calc_inc <- function(risk, inc, prop_high, lambda) {
  
  return((1-exp(-lambda*365))*(1-prop_high) + (1-exp(-risk*lambda*365))*prop_high - inc)
  
}

si_ode <- function(times, init, param){
  with(as.list(c(init, param)), {
    #browser()
    
    # Flows
    # the number of people moving from S to I at each time step
    #Susceptible, Infected, placebo
    SIp.flow <- lambda*Sp
    SIv.flow <- lambda*(1-epsilon)*Sv
    
    #Susceptible, Infected, placebo, high, medium, low
    SIph.flow <- risk*lambda*Sph
    SIpl.flow <- lambda*Spl  
    
    #Susceptible, Infected, vaccine, high, medium, low
    SIvh.flow <- risk*lambda*(1-epsilon)*Svh
    SIvl.flow <- lambda*(1-epsilon)*Svl
    
    # ODEs
    # placebo; homogeneous risk
    dSp <- -SIp.flow
    dIp <- SIp.flow  #lambda*Sp
    
    # vaccine; homogeneous risk
    dSv <- -SIv.flow
    dIv <- SIv.flow  #lambda*epsilon*Sv
    # placebo; heterogeneous risk
    dSph <- -SIph.flow
    dIph <- SIph.flow  #risk*lambda*Sph
    dSpl <- -SIpl.flow
    dIpl <- SIpl.flow  #0*lambda*Spl
    
    # vaccine; heterogeneous risk
    dSvh <- -SIvh.flow
    dIvh <- SIvh.flow  #risk*lambda*(1-epsilon)*Svh
    dSvl <- -SIvl.flow
    dIvl <- SIvl.flow  #0*lambda*(1-epsilon)*Svl
    
    #Output
    list(c(dSp,dIp,
           dSv,dIv,
           dSph,dIph,
           dSpl,dIpl,
           dSvh,dIvh,
           dSvl,dIvl,
           SIp.flow,SIv.flow,
           SIph.flow,SIpl.flow,
           SIvh.flow,SIvl.flow))
  })
}

mod.manipulate <- function(mod){
  
  mod <- mutate_epi(mod, total.Svh.Svl = Svh + Svl) #all susceptible in heterogeneous risk vaccine pop
  mod <- mutate_epi(mod, total.Sph.Spl = Sph + Spl) #all susceptible in heterogeneous risk placebo pop
  mod <- mutate_epi(mod, total.Ivh.Ivl = Ivh + Ivl) #all infected in heterogeneous risk vaccine pop
  mod <- mutate_epi(mod, total.Iph.Ipl = Iph + Ipl) #all infected in heterogeneous risk placebo pop
  mod <- mutate_epi(mod, total.SIvh.SIvl.flow = SIvh.flow + SIvl.flow) #all infections per day in heterogeneous risk vaccine pop
  mod <- mutate_epi(mod, total.SIph.SIpl.flow = SIph.flow + SIpl.flow) #all infections in heterogeneous risk placebo pop
  #Incidence estimates, per 100 person years
  #Instantaneous incidence / hazard
  mod <- mutate_epi(mod, rate.Vaccine = (SIv.flow/Sv)*365*100)
  mod <- mutate_epi(mod, rate.Placebo = (SIp.flow/Sp)*365*100)
  mod <- mutate_epi(mod, rate.Vaccine.het = (total.SIvh.SIvl.flow/total.Svh.Svl)*365*100)
  mod <- mutate_epi(mod, rate.Placebo.het = (total.SIph.SIpl.flow/total.Sph.Spl)*365*100)
  #Cumulative incidence
  mod <- mutate_epi(mod, cumul.Iv = Iv)
  mod <- mutate_epi(mod, cumul.Ip = Ip)
  mod <- mutate_epi(mod, cumul.Ivh.Ivl = total.Ivh.Ivl)
  mod <- mutate_epi(mod, cumul.Iph.Ipl = total.Iph.Ipl)
  
  n <- mod$param$n
  
  mod <- mutate_epi(mod, cumul.rate.Vaccine = (cumul.Iv/n)*100)
  mod <- mutate_epi(mod, cumul.rate.Placebo = (cumul.Ip/n)*100)
  mod <- mutate_epi(mod, cumul.rate.Vaccine.het = (cumul.Ivh.Ivl/n)*100)
  mod <- mutate_epi(mod, cumul.rate.Placebo.het = (cumul.Iph.Ipl/n)*100)
  # mod <- mutate_epi(mod, cumul.Sv = cumsum(Sv)
  # mod <- mutate_epi(mod, cumul.Sp = cumsum(Sp))
  # mod <- mutate_epi(mod, cumul.Svh.Svl = cumsum(total.Svh.Svl))
  # mod <- mutate_epi(mod, cumul.Sph.Spl = cumsum(total.Sph.Spl))
  # mod <- mutate_epi(mod, cumul.rate.Vaccine = (Iv/cumul.Sv)*365*100)
  # mod <- mutate_epi(mod, cumul.rate.Placebo = (Ip/cumul.Sp)*365*100)
  # mod <- mutate_epi(mod, cumul.rate.Vaccine.het = (total.Ivh.Ivl/cumul.Svh.Svl)*365*100)
  # mod <- mutate_epi(mod, cumul.rate.Placebo.het = (total.Iph.Ipl/cumul.Sph.Spl)*365*100)
  #Vaccine efficacy (VE) estimates
  #VE <- 1 - Relative Risk; this is VE for instantaneous incidence / hazard
  mod <- mutate_epi(mod, VE1.inst = 1 - rate.Vaccine/rate.Placebo)
  mod <- mutate_epi(mod, VE2.inst = 1 - rate.Vaccine.het/rate.Placebo.het)
  #VE <- 1 - Relative Risk; this is VE from cumulative incidence
  mod <- mutate_epi(mod, VE1.cumul = 1 - cumul.rate.Vaccine/cumul.rate.Placebo)
  mod <- mutate_epi(mod, VE2.cumul = 1 - cumul.rate.Vaccine.het/cumul.rate.Placebo.het)
  #return(mod)
}
