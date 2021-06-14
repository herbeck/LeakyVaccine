
#------------------------------------------------------------------------------
# Right now this includes just two vaccine trial populations, each with a 
# vaccine arm and a placebo arm. One population has homogeneous exposure / risk
# of infection; the other population includes exposure heterogeneity, and this 
# heterogeneity is the same in both trial arms.
#------------------------------------------------------------------------------
si_ode <- function(times, init, param){
  with(as.list(c(init, param)), {
    
    # Flows
    # the number of people moving from S to I at each time step
    #Susceptible, Infected, placebo
    SIp.flow <- lambda*Sp
    SIv.flow <- lambda*(1-epsilon)*Sv
    
    #Susceptible, Infected, placebo, high, medium, low
    SIph.flow <- risk*lambda*Sph
    SIpm.flow <- lambda*Spl
    #SIpl.flow <- 0*lambda*Spl  #zero risk in the low risk group
    SIpl.flow <- (1/risk)*lambda*Spl  #inverse risk multiplier for low risk group
    
    #Susceptible, Infected, vaccine, high, medium, low
    SIvh.flow <- risk*lambda*(1-epsilon)*Svh
    SIvm.flow <- lambda*(1-epsilon)*Spl
    #SIvl.flow <- 0*lambda*(1-epsilon)*Svl  #zero risk in the low risk group 
    SIvl.flow <- (1/risk)*lambda*(1-epsilon)*Svl  #inverse risk 
    
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
    dSpm <- -SIpm.flow
    dIpm <- SIpm.flow  #lambda*Spm
    dSpl <- -SIpl.flow
    dIpl <- SIpl.flow  #0*lambda*Spl
    
    # vaccine; heterogeneous risk
    dSvh <- -SIvh.flow
    dIvh <- SIvh.flow  #risk*lambda*(1-epsilon)*Svh
    dSvm <- -SIvm.flow
    dIvm <- SIvm.flow  #lambda*Svm
    dSvl <- -SIvl.flow
    dIvl <- SIvl.flow  #0*lambda*(1-epsilon)*Svl
    
    #Output
    list(c(dSp,dIp,
           dSv,dIv,
           dSph,dIph,
           dSpm,dIpm,
           dSpl,dIpl,
           dSvh,dIvh,
           dSvm,dIvm,
           dSvl,dIvl,
           SIp.flow,SIv.flow,
           SIph.flow,SIpm.flow,SIpl.flow,
           SIvh.flow,SIvm.flow,SIvl.flow))
  })
}


#------------------------------------------------------------------------------
# Initial parameter settings
#------------------------------------------------------------------------------
# beta <- 0.004   #transmission rate (per contact)
# c <- 90/365    #contact rate (contacts per day)
# prev <- 0.10    #needs some consideration
# #prev <- 0.01
# lambda <- beta*c*prev
# #lambda <- 0.000008398975
# epsilon <- 0.30  #per contact vaccine efficacy
# risk <- 10.0   #risk multiplier


runSim <- function(param) {
  
  param <- param.dcm(lambda = param$beta * param$contactRate * param$prev, 
                     epsilon = param$epsilon, 
                     risk = param$risk)
  init <- init.dcm(Sp = 10000, Ip = 0,
                   Sv = 10000, Iv = 0,
                   Sph = 1000, Iph = 0,    #placebo, high risk
                   Spm = 7000, Ipm = 0,    #placebo, medium risk
                   Spl = 2500, Ipl = 0,    #placebo, low risk
                   Svh = 1000, Ivh = 0,    #vaccine, high risk
                   Svm = 7000, Ivm = 0,    #vaccine, medium risk
                   Svl = 2500, Ivl = 0,    #vaccine, low risk
                   SIp.flow = 0, SIv.flow = 0, 
                   SIph.flow = 0, SIpm.flow = 0, SIpl.flow = 0,
                   SIvh.flow = 0, SIvm.flow = 0, SIvl.flow = 0)
  
  control <- control.dcm(nsteps = 365*3, new.mod = si_ode)
  mod <- dcm(param, init, control)
  #mod
  
  mod <- mod.manipulate(mod)
  return (mod)
}


#------------------------------------------------------------------------------
# This function just takes the model output (a `mod` file) and uses the data to 
# create other data (e.g. incidence and VE estimates) for plotting.
#------------------------------------------------------------------------------
mod.manipulate <- function(mod){
  
  mod <- mutate_epi(mod, total.Svh.Svm.Svl = Svh + Svm + Svl) #all susceptible in heterogeneous risk vaccine pop
  mod <- mutate_epi(mod, total.Sph.Spm.Spl = Sph + Spm + Spl) #all susceptible in heterogeneous risk placebo pop
  mod <- mutate_epi(mod, total.Ivh.Ivm.Ivl = Ivh + Ivm + Ivl) #all infected in heterogeneous risk vaccine pop
  mod <- mutate_epi(mod, total.Iph.Ipm.Ipl = Iph + Ipm + Ipl) #all infected in heterogeneous risk placebo pop
  mod <- mutate_epi(mod, total.SIvh.SIvm.SIvl.flow = SIvh.flow + SIvm.flow + SIvl.flow) #all infections per day in heterogeneous risk vaccine pop
  mod <- mutate_epi(mod, total.SIph.SIpm.SIpl.flow = SIph.flow + SIpm.flow + SIpl.flow) #all infections in heterogeneous risk placebo pop
  
  #Instantaneous ncidence (hazard) estimates, per 100 person years
  #Instantaneous incidence / hazard
  mod <- mutate_epi(mod, rate.Vaccine = (SIv.flow/Sv)*365*100)
  mod <- mutate_epi(mod, rate.Placebo = (SIp.flow/Sp)*365*100)
  mod <- mutate_epi(mod, rate.Vaccine.het = (total.SIvh.SIvm.SIvl.flow/total.Svh.Svm.Svl)*365*100)
  mod <- mutate_epi(mod, rate.Placebo.het = (total.SIph.SIpm.SIpl.flow/total.Sph.Spm.Spl)*365*100)
  
  #Cumulative incidence
  mod <- mutate_epi(mod, cumul.Sv = cumsum(Sv))
  mod <- mutate_epi(mod, cumul.Sp = cumsum(Sp))
  mod <- mutate_epi(mod, cumul.Svh.Svm.Svl = cumsum(total.Svh.Svm.Svl))
  mod <- mutate_epi(mod, cumul.Sph.Spm.Spl = cumsum(total.Sph.Spm.Spl))
  mod <- mutate_epi(mod, cumul.rate.Vaccine = (Iv/cumul.Sv)*365*100)
  mod <- mutate_epi(mod, cumul.rate.Placebo = (Ip/cumul.Sp)*365*100)
  mod <- mutate_epi(mod, cumul.rate.Vaccine.het = (total.Ivh.Ivm.Ivl/cumul.Svh.Svm.Svl)*365*100)
  mod <- mutate_epi(mod, cumul.rate.Placebo.het = (total.Iph.Ipm.Ipl/cumul.Sph.Spm.Spl)*365*100)
  
  #Vaccine efficacy (VE) estimates
  #VE <- 1 - Relative Risk; this is VE for hazard
  mod <- mutate_epi(mod, VE1.inst = 1 - rate.Vaccine/rate.Placebo)
  mod <- mutate_epi(mod, VE2.inst = 1 - rate.Vaccine.het/rate.Placebo.het)
  
  #VE <- 1 - Relative Risk; this is VE from cumulative incidence
  mod <- mutate_epi(mod, VE1.cumul = 1 - cumul.rate.Vaccine/cumul.rate.Placebo)
  mod <- mutate_epi(mod, VE2.cumul = 1 - cumul.rate.Vaccine.het/cumul.rate.Placebo.het)
  
  return(mod)
}


