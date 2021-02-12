
## MIDAS Vaccine trial model, with EpiModel

# Original model ----------------------------------------------------------

library("deSolve")
library("tidyverse")

p = 0.01   #transmission rate (per contact)
c = 0.3   #contact rate (contacts per time)
Prev <- 0.20
lambda <- p*c*Prev
epsilon <- 0.50 #per act vaccine efficacy
risk <- 3.0 #risk multiplier

#parms <- c(lambda, epsilon, risk)
init <- c(Sv=1000,Iv=1,Sp=1000,Ip=1,Svh=500,Ivh=1,Svl=500,Ivl=1)
times <- seq(from=1, to=365*3, by=1)

sir_ode <- function(times, init, parms){
  with(as.list(c(init, parms)), {
    
    # Flows
    SIp.flow <- lambda*Sp
    SIv.flow <- lambda*epsilon*Sv
    SIvh.flow <- risk*lambda*(1-epsilon)*Svh
    SIvl.flow <- lambda*(1-epsilon)*Svl
    
    # ODEs
    #placebo
    #Np = Sp+Ip
    dSp <- -lambda*Sp
    dIp <- SIp.flow  #lambda*Sp
    
    # vaccine; homogeneous risk
    #Nv <- Sv+Iv
    dSv <- -lambda*epsilon*Sv
    dIv <- SIv.flow  #lambda*epsilon*Sv

    # vaccine; heterogeneous risk
    #Nv.risk <- Svh+Ivh+Svl+Ivl
    dSvh <- -risk*lambda*(1-epsilon)*Svh
    dIvh <-  SIvh.flow  #risk*lambda*(1-epsilon)*Svh
    dSvl <- -lambda*(1-epsilon)*Svl
    dIvl <- SIvl.flow  #lambda*(1-epsilon)*Svl

    #Output
    list(c(dSv,dIv,dSp,dIp,dSvh,dIvh,dSvl,dIvl,SIp.flow,SIv.flow,SIvh.flow,SIvl.flow))
    
  })
}

#sir_out <- lsoda(init, times, sir_ode, parms)

#sir_out_inf <- sir_out %>%
#  as.data.frame() %>%
#  gather("variable", "value", -time) %>%
#  filter(variable %in% c('Iv', 'Ip','Ivl','Ivh'))

#ggplot(sir_out_inf, aes(x = time/365, y = value, color = variable)) +
  #Add line
#  geom_line(lwd = 2) +
  #Add labels
#  xlab("Years") + ylab("Infections")


# EpiModel Conversion -----------------------------------------------------

library("EpiModel")

params <- param.dcm(lambda = lambda, epsilon = epsilon, risk = risk)
inits <- init.dcm(Sv = 1000, Iv = 1, Sp = 1000, Ip = 1,
                  Svh = 500, Ivh = 1, Svl = 500, Ivl = 1)
controls <- control.dcm(nsteps = 365*3, new.mod = sir_ode, print.mod = T)

mod <- dcm(params, inits, controls)
mod

par(mar = c(3,3,2,1), mgp = c(2,1,0))
plot(mod, y = c("Iv", "Ivh", "Ip", "Ivl"), 
     legend = "full")
plot(mod, y = c("Iv", "Ivh", "Ip", "Ivl"), 
     alpha = 0.8, 
     main = "Cumulative infections")

df <- as.data.frame(mod)
head(df, 25)



# Example of more sophisticated risk heterogeneity model ------------------

# This is an SIS (e.g., gonorrhea) model I use in ID modeling teaching that may
# be helpful... this represents a "core group" with higher partnership
# acquisition rates but with assortative mixing by those rates.

# This follows the SIS gonorrhea model in the textbook, in which mixing may
# range from assortative to proportional to disasortative along the continuous
# gradient defined by the Q statistic

# This model will be run in a large population of 20 million people, in which
# 2% are high-risk as defined by their contact rates (the core group) and the
# other 98% are low-risk
N <- 2e7
prop.high <- 0.02
prop.low <- 0.98

# Translate these proportions into population numbers
S.high <- N*prop.high
S.low <- N*prop.low

I.high <- 1
I.low <- 0

N.high <- S.high + I.high
N.low <- S.low + I.low

# For the contact rates, recall that we will specify an overall mean contact rate
# and also a contact rate in the low-risk group. Then we will calculate
c.mean <- 2
c.low <- 1.4

# To handle balancing, we will solve for the high-risk contact rate with these two
# c.mean*N <- c.high*N.high + c.low*N.low
c.high <- (c.mean*N - c.low*N.low)/N.high
c.high

# check against the original equation
c.mean*N == c.high*N.high + c.low*N.low

# Now, let's start with the mixing matrix. Recall this is our four-cell matrix
# of mixing between high and low risk persons.
#       h.j    l.j
#      ------ ------
# h.i | g.hh | g.hl |
#      ------ ------
# l.i | g.lh | g.ll |
#      ------ ------
# where h is high-risk and l is low-risk, and g is the proportion in each
# cell of the matrix.

# In proportional mixing, the proportion of partnerships of high-risk people
# that are with other high-risk people is the fraction of high-risk contacts
# available over the total number of contacts available, where the contact
# number is the product of the contact rate and the group size

g.hh <- (c.high*N.high)/(c.high*N.high + c.low*N.low)
g.hh

# Because the columns sum to 1, g.lh = 1-g.hh
g.lh <- 1 - g.hh
g.lh

# Next, we solve for g.hl. Because of balancing, we need g.hl*c.low*N.low
# must be equal to g.lh*c.high*N.high. Therefore, rearranging:
g.hl <- g.lh * ((c.high*N.high) / (c.low*N.low))
g.hl

# Finally, g.ll is just the inverse of g.hl
g.ll <- 1 - g.hl
g.ll

m <- matrix(c(g.hh, g.lh, g.hl, g.ll), ncol = 2, nrow = 2)
colnames(m) <- c("h.j", "l.j")
rownames(m) <- c("h.i", "l.i")
m

# Try a different value for the mean contact rate above and building
# the proportional mixing matrix again. Then put it back to 2.

# Next, let's input the Q statistic into this matrix. Recall, that when when
# Q = 1 there is purely assortative mixing, like this:

#       h.j    l.j
#      ------ ------
# h.i |  1   |  0   |
#      ------ ------
# l.i |  0   |  1   |
#      ------ ------

# When Q=0, we have proportional mixing. We will implement Q as a continuous
# number input into the g.hh cell of the matrix. When Q=0, we already calculated
# the value of g.hh:
g.hh <- (c.high*N.high)/(c.high*N.high + c.low*N.low)
g.hh

# Therefore:
Q <- (0.314 + 0.686 - 1)/(2 - 1)
Q

# When Q = 1 then g.hh = 1. With the denominator fixed, we need to update the
# fraction to shift more weight into the numerator. To make it 1, we need:
g.hh <- (c.high*N.high + c.low*N.low)/(c.high*N.high + c.low*N.low)
g.hh

# Therefore, Q will be implemented as a proportional element to scale
# c.low*N.low in the numerator:
Q <- 1
g.hh <- ((c.high*N.high) + (Q*c.low*N.low)) / ((c.high*N.high) + (c.low*N.low))
g.hh

# then the rest of our mixing matrix:
g.lh <- 1 - g.hh
g.hl <- g.lh * ((c.high*N.high) / (c.low*N.low))
g.ll <- 1 - g.hl
m <- matrix(c(g.hh, g.lh, g.hl, g.ll), ncol = 2, nrow = 2)
colnames(m) <- c("h.j", "l.j"); rownames(m) <- c("h.i", "l.i")
m

# Back to proportional mixing.
Q <- 0
g.hh <- ((c.high*N.high) + (Q*c.low*N.low)) / ((c.high*N.high) + (c.low*N.low))
g.lh <- 1 - g.hh
g.hl <- g.lh * ((c.high*N.high) / (c.low*N.low))
g.ll <- 1 - g.hl
m <- matrix(c(g.hh, g.lh, g.hl, g.ll), ncol = 2, nrow = 2)
colnames(m) <- c("h.j", "l.j"); rownames(m) <- c("h.i", "l.i")
m

# Voila!

# For disasortative mixing, if the contact rates between the two groups were
# closer or the population sizes were closer, we might get to a negative Q
# value approaching -1. But there are simply too many low risk people relative
# to high risk people in the population. The best we can do is set g.hh to 0
# and then continue with the mixing matrix:
g.hh <- 0
g.lh <- 1 - g.hh
g.hl <- g.lh * ((c.high*N.high) / (c.low*N.low))
g.ll <- 1 - g.hl
m <- matrix(c(g.hh, g.lh, g.hl, g.ll), ncol = 2, nrow = 2)
colnames(m) <- c("h.j", "l.j"); rownames(m) <- c("h.i", "l.i")
m

# All high-risk persons can mixing with low-risk persons, but there are too
# many low-risk persons (given their contact rate) to mix only with
# high-risk persons (there are excess implied partnerships). Therefore, they
# have to partner with some low-risk persons too. The Q value associated with
# this mixing matrix is:
Q <- (0 + 0.54 - 1)/(2 - 1)
Q

# That is the lower bound of Q for this particular set of contact rate parameters
# and population sizes. We can get to purely disasortative mixing by lowering
# the contact rate for the low-risk group to 1.02. Try that on your own.

# Moving on the model:
Qmod <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {

    ## Dynamic Calculations ##

    # Popsize
    N.high <- S.high + I.high
    N.low <- S.low + I.low

    # Contact rates
    c.high <- (c.mean*N - c.low*N.low)/N.high

    # mixing matrix calculations based on Q
    g.hh <- ((c.high*N.high) + (Q*c.low*N.low)) / ((c.high*N.high) + (c.low*N.low))
    g.lh <- 1 - g.hh
    g.hl <- (1 - g.hh) * ((c.high*N.high) / (c.low*N.low))
    g.ll <- 1 - g.hl

    # prob that p is infected, as a weighted average of the probability of
    # selecting a partner in each group and the prevalence in that group
    p.high <- (g.hh*I.high/N.high) + (g.lh*I.low/N.low)
    p.low <- (g.ll*I.low/N.low) + (g.hl*I.high/N.high)

    # lambda
    lambda.high <- tau * c.high * p.high
    lambda.low <- tau * c.low * p.low

    gamma <- 1/dur.inf

    ## Differential Equations ##
    dS.high <- -lambda.high*S.high + gamma*I.high
    dI.high <- lambda.high*S.high - gamma*I.high

    dS.low <- -lambda.low*S.low + gamma*I.low
    dI.low <- lambda.low*S.low - gamma*I.low


    ## Output ##
    list(c(dS.high, dI.high,
           dS.low, dI.low))

  })
}

# We parameterize this following the values suggested in the text. Note here
# for the first time we will use yearly time steps, so the model parameters
# are specified in terms of yearly rate. However, we will solve for the model
# in roughly weekly intervals by specifying dt = 1/50. We could do the same
# thing by converting the rates to weekly, then running the model for 5000
# time steps. The base model value for Q will be 0.
param <- param.dcm(c.mean = 2, c.low = 1.4, tau = 0.75, dur.inf = 2/12, Q = 0)
init <- init.dcm(S.high = S.high, I.high = I.high, S.low = S.low, I.low = I.low)
control <- control.dcm(nsteps = 100, dt = 0.02, new.mod = Qmod)

# Run the model and examine the output
mod <- dcm(param, init, control)
mod

# Add additional stats for the population sizes, and then the prevalence overall
# and within each risk group
mod <- mutate_epi(mod, N.high = S.high + I.high,
                  N.low = S.low + I.low,
                  N = S.high + I.high + S.low + I.low)
mod <- mutate_epi(mod, prev = (I.high + I.low)/N,
                  prev.h = I.high/N.high,
                  prev.l = I.low/N.low)
mod

# Plot the prevalence results together from this proportional mixing model
plot(mod, y = c("prev", "prev.h", "prev.l"), legend = "full", ylim = c(0, 0.4))

# Let's explore how mixing impacts disease transmission by varying Q over a
# broad range of values from disasortative to proportional to asortative.
param <- param.dcm(c.mean = 2, c.low = 1.4, tau = 0.75, dur.inf = 2/12,
                   Q = c(-0.46, 0, 0.33, 0.75, 1))
init <- init.dcm(S.high = S.high, I.high = I.high, S.low = S.low, I.low = I.low)
control <- control.dcm(nsteps = 100, dt = 0.02, new.mod = Qmod)
mod <- dcm(param, init, control)

# Add in the same summary statistics as before.
mod <- mutate_epi(mod, N.high = S.high + I.high,
                  N.low = S.low + I.low,
                  N = S.high + I.high + S.low + I.low)
mod <- mutate_epi(mod, prev = (I.high + I.low)/N,
                  prev.h = I.high/N.high,
                  prev.l = I.low/N.low)

# Plot the results. Some very interesting findings here. Let's think through how
# we arrived at them.
par(mfrow = c(1,3))
plot(mod, y = "prev", legend = "full", alpha = 0.75,
     main = "Overall Prev", ylim = c(0, 0.04), xlim = c(0, 20))
plot(mod, y = "prev.h", legend = "full", alpha = 0.75,
     main = "Prev High", ylim = c(0, 0.8), xlim = c(0, 20))
plot(mod, y = "prev.l", legend = "full", alpha = 0.75,
     main = "Prev Low", ylim = c(0, 0.03), xlim = c(0, 20))
