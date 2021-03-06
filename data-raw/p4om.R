# p4om.R - DESC
# mse/data-raw/p4om.R

# Copyright European Union, 2018
# Author: Iago Mosqueira (EC JRC) <iago.mosqueira@ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

library(mse)
library(FLa4a)
library(FLasher)
library(FLBRP)

data(ple4)
data(ple4.index)

stk <- ple4
idx <- FLIndices(SURVEY=ple4.index)

# Variables

it <- 50 # iterations
fy <- 2030 # final year
y0 <- range(stk)["minyear"] # initial data year
dy <- range(stk)["maxyear"] # final data year
iy <- 2008 # initial year of projection (also intermediate)
ny <- fy - iy + 1 # number of years to project from intial year
nsqy <- 3 # number of years to compute status quo metrics
vy <- ac(iy:fy) # vector of years to be projected

# Operating model conditioning

# fit stock assessment model
mcsave <- 200
mcmc <- mcsave*it
fit <- sca(stk, idx, fit="MCMC",
  mcmc = SCAMCMC(mcmc = mcmc, mcsave = mcsave, mcprobe = 0.4))

stk <- slim(stk + fit)

# Make SRRs

# average recruitment estimation sd
rv1 <- sqrt(mean(c(iterVars(log(rec(stk)))), na.rm=TRUE))

# average autocor lag1
ac1 <- mean(apply(window(rec(stk), end=2008)@.Data,6,function(x) c(acf(c(x), plot=FALSE, lag.max=1)$acf[2])))

# BevHolt
srbh <- fmle(as.FLSR(stk, model="bevholt"), method="L-BFGS-B", lower=c(1e-6, 1e-6), upper=c(max(rec(stk)) * 3, Inf))

# Deviances
devbh <- ar1rlnorm(rho=ac1, years=dy:fy, iters=it, margSD=rv1*2)
residuals(srbh) <- devbh

# Refpts

brp <- brp(FLBRP(stk, srbh))

# Set up future assumptions - means of 5 years
stk <- fwdWindow(stk, brp, end=fy)

# Fleet behaviour

fb <- mseCtrl(method=hyperstability.fb, args=list(beta=0.8))

# OM object
om <- FLom(stock=stk, sr=srbh, refpts=refpts(brp))

save(om, file="../data/p4om.RData", compress="xz")

# OEM

# Indices

# IEM

# save(om, oem, iem, file="../data/p4om.RData", compress="xz")
