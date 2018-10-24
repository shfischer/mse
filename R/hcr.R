# hcr.R - DESC
# mse/R/hcr.R

# Copyright European Union, 2018
# Author: Ernesto Jardim (EC JRC) <ernesto.jardim@ec.europa.eu>
#         Iago Mosqueira (EC JRC) <iago.mosqueira@ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

# Harvest Control Rule function: h()
#' Evaluate the chosen HCR function
#'
#' Evaluate the chosen HCR function using the current stock perception and a control
#' @param method Name of the chosen HCR function.
#' @param stk The perceived stock.
#' @param ay The current year. The management control (e.g. F or effort) will be set in ay+1.
#' @param EFF Effort array (if effort management is being used).
#' @param EFF0 Tracking array.
#' @param control The control object for the chosen HCR function. A list of parameters.
#h <- function(...){
#	args <- list(...)
#	method <- args$method
#	args$method <- NULL
#	# Check inputs
#	if(!is(args$stk,"FLStock")) stop("stk argument must be an FLStock")
#	# dispatch
#	ctrl <- do.call(method, args)
#	# check outputs
#	if(!is(ctrl, "fwdControl")) stop("The HCR must return and object of class fwdControl")	
#	# return
#	ctrl  
#}

# ices.hcr {{{

#' The typical HCR used by ICES
#'
#' The typical HCR used by ICES which sets a target F based on the SSB based on 4 parameters: blim, bsafe, fmin and ftrg.
#' F increases linearly between SSB = blim and SSB = bsafe, from F = fmin to F = ftrg.
#' If:
#' B < Blim, F = Fbycatch;
#' B > trigger, F = Fmsy;
#' B > Blim & B < trigger, F linear between Fbycatch and Fmsy;
#' F = ftrg is the maximum F, F = fmin is the minimum F.
#' F is set in year ay, based on SSB in year ay - ssb_lag.
#' The control argument is a list of parameters used by the HCR.
#' @param stk The perceived FLStock.
#' @param control A list with the elements fmin, ftrg, blim, bsafe and ssb_lag, all of which are numeric.
#' @param ay The year for which the target F is set, based on the SSB in year (ay - control$ssb_lag).
ices.hcr <- function(stk, fmin, ftrg, blim, bsafe, ssb_lag=1, ay, tracking){
	# rule
  ssb <- ssb(stk)[, ac(ay-ssb_lag)]
  fout <- FLQuant(fmin, dimnames=list(iter=dimnames(ssb)$iter))
  fout[ssb >= bsafe] <- ftrg
  inbetween <- (ssb < bsafe) & (ssb > blim)
  gradient <- (ftrg - fmin) / (bsafe - blim)
  fout[inbetween] <- (ssb[inbetween] - blim) * gradient + fmin
	# create control file
	ctrl <- getCtrl(c(fout), "f", ay+1, dim(fout)[6])
	# return
	list(ctrl=ctrl, tracking=tracking)
} # }}}

# fixedF.hcr {{{

#' A fixed target f
#'
#' No matter what get F = Ftarget
#' The control argument is a list of parameters used by the HCR.
#' @param stk The perceived FLStock.
#' @param control A list with the element ftrg (numeric).
fixedF.hcr <- function(stk, ftrg, ay, tracking){
	
  # rule 
	if(!is(ftrg, "FLQuant"))
    ftrg <- FLQuant(ftrg, dimnames=list(iter=dimnames(stk@catch)$iter))

	# create control file
	ctrl <- getCtrl(c(ftrg), "f", ay+1, dim(ftrg)[6])
	
  # return
	list(ctrl=ctrl, tracking=tracking)
} # }}}

# movingF.hcr {{{

movingF.hcr <- function(stk, hcrpars, ay, tracking){

	# rule 
	if(!is(hcrpars, "FLQuant"))
    hcrpars <- FLQuant(hcrpars, dimnames=list(iter=dimnames(stk@catch)$iter))
	
  # create control file
	ctrl <- getCtrl(c(hcrpars), "f", ay+1, dim(hcrpars)[6])
	
  # return
	list(ctrl=ctrl, tracking=tracking)
} # }}}

# catchSSB.hcr {{{

#' A HCR to set catch based on SSB
#'
# hcrparams=FLPar(dlimit=0.10, dtarget=0.40, lambda=1.0, dltac=0.15, dhtac=0.15),

#' @param dtarget
#' @param dlimit
#' @param lambda
#' @param MSY

catchSSB.hcr <- function(stk, dtarget=0.40, dlimit=0.10, lambda=1, MSY, ssb_lag=1,
  ay, tracking) {
  
  # COMPUTE depletion
  dep <- ssb(stk)[, ac(ay - ssb_lag)] / ssb(stk)[, 1]

  # RULE
  ca <- ifelse(dep <= dlimit, 0,
    ifelse(dep < dtarget, (lambda * MSY) / (dtarget - dlimit) * (dep - dlimit),
    MSY))
  
  # CONTROL
	ctrl <- getCtrl(c(ca), quant="catch", years=ay + 1, it=dim(ca)[6])

  # TAC limits
	
	return(list(ctrl=ctrl, tracking=tracking))

} # }}}

#' cpue.hcr
#'
#' @examples
#' data(ple4)

cpue.hcr <- function(stk, rule=~tac * (1 + lambda * slope), ay,
  lambda=1, tracking){
  
  slope <- tracking["cpue.est", ac(ay)]

  # TODO getCtrl
  ctrl <- fwdControl(quant="catch", value=eval(rule[[2]],
    list(tac=catch(stk)[, ac(ay-1)], lambda=lambda, slope=slope)), year=ay+1)
  
	return(list(ctrl=ctrl, tracking=tracking))
}

  
