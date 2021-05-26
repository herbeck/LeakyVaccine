
## MIDAS Vaccine trial model, with EpiModel

library("deSolve")
library("tidyverse")
library("EpiModel")

# Parameters
p <- 0.01   #transmission rate (per contact)
c <- 0.25   #contact rate (contacts per time)
prev <- 0.20

lambda <- p*c*prev
epsilon <- 0.50 #per act vaccine efficacy
risk <- 3.0 #risk multiplier

si_ode <- function(times, init, param){
  with(as.list(c(init, param)), {
    
    # Flows
    # the number of people moving from S to I at each time step
    SIp.flow <- lambda*Sp
    SIv.flow <- lambda*epsilon*Sv
    SIvh.flow <- risk*lambda*(1-epsilon)*Svh
    SIvl.flow <- lambda*(1-epsilon)*Svl
    
    # ODEs
    #placebo; homogeneous risk
    dSp <- -SIp.flow
    dIp <- SIp.flow  #lambda*Sp
    
    # vaccine; homogeneous risk
    dSv <- -SIv.flow
    dIv <- SIv.flow  #lambda*epsilon*Sv

    # vaccine; heterogeneous risk
    dSvh <- -SIvh.flow
    dIvh <- SIvh.flow  #risk*lambda*(1-epsilon)*Svh
    dSvl <- -SIvl.flow
    dIvl <- SIvl.flow  #lambda*(1-epsilon)*Svl

    #Output
    list(c(dSp,dIp,
           dSv,dIv,
           dSvh,dIvh,
           dSvl,dIvl,
           SIp.flow,SIv.flow,SIvh.flow,SIvl.flow))
  })
}

# EpiModel part

param <- param.dcm(lambda = lambda, epsilon = epsilon, risk = risk)
init <- init.dcm(Sp = 5000, Ip = 1,
                 Sv = 5000, Iv = 1,
                 Svh = 2500, Ivh = 1, 
                 Svl = 2500, Ivl = 1,
                 SIp.flow = 0, SIv.flow = 0, SIvh.flow = 0, SIvl.flow = 0)

control <- control.dcm(nsteps = 365*3, new.mod = si_ode)

mod <- dcm(param, init, control)
mod

par(mar = c(3,3,2,1), mgp = c(2,1,0))
plot(mod, y = c("Iv", "Ivh", "Ip", "Ivl"), 
     alpha = 0.8, 
     main = "Cumulative infections",
     legend = "full")
plot(mod, y = c("SIp.flow", "SIv.flow", "SIvh.flow", "SIvl.flow"), 
     alpha = 0.8, 
     main = "Incidence",
     legend = "full")

df <- as.data.frame(mod)
head(df, 25)




