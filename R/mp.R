# mp.R - DESC
# mse/R/mp.R

# Copyright European Union, 2018
# Author: Ernesto Jardim (EC JRC) <ernesto.jardim@ec.europa.eu>
#         Iago Mosqueira (EC JRC) <iago.mosqueira@ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

# mp {{{

mp <- function(om, oem=FLoem(), iem="missing", ctrl.mp, genArgs, 
  scenario="test", tracking="missing", verbose=TRUE){

	# PREPARE the om
	stk.om <- stock(om)	
	name(stk.om) <- scenario
	sr.om <- sr(om)
	sr.om.res <- residuals(sr.om)

  # VARIABLES
	fy <- genArgs$fy # final year
  # dy0 om$startyr | oem$args$
	y0 <- genArgs$y0 # initial data year
  # ay - oem$args$datalag
	dy <- genArgs$dy # final data year
  # TODO mlag
	iy <- genArgs$iy # initial year of projection (also intermediate)
  # om$its
  it <- genArgs$it
  #
	nsqy <- genArgs$nsqy # number of years to compute status quo metrics
  # freq

	vy <- ac(iy:fy) # vector of years to be projected

	# INIT tracking
  metric <- c("F.est", "B.est", "conv.est", "metric.hcr", "metric.is",
    "metric.iem", "metric.fb","F.om", "B.om", "C.om")
	
  if (!missing(tracking))
    metric <- c(metric, tracking)

	tracking <- FLQuant(NA, dimnames=list(metric=metric, year=c(iy-1,vy), iter=1:it))

	# GET historical
  # TODO CHECK
	tracking["metric.is", ac(iy)] <- catch(stk.om)[,ac(iy)]

	# SET seed
	if (!is.null(genArgs$seed)) set.seed(genArgs$seed)
  
	# --- GO FISH!
  for(i in vy[-length(vy)]) {

		gc()

    if(verbose)
  		cat(i, " > ")
		
    ay <- an(i)
    vy0 <- 1:(ay-y0) # data years (positions vector) - one less than current year
		sqy <- ac((ay-1):(ay-nsqy)) # years for status quo computations 
		
    # TRACK om
		tracking["F.om", ac(ay-1)] <- fbar(stk.om)[,ac(ay-1)]    
		tracking["B.om", ac(ay-1)] <- ssb(stk.om)[,ac(ay-1)]    
		tracking["C.om", ac(ay-1)] <- catch(stk.om)[,ac(ay-1)]    
		
		# --- OEM
		ctrl.oem <- args(oem)
		ctrl.oem$method <- method(oem)
		ctrl.oem$deviances <- deviances(oem)
		ctrl.oem$observations <- observations(oem)
		ctrl.oem$stk <- stk.om
		ctrl.oem$vy0 <- vy0
		ctrl.oem$ay <- ay
		ctrl.oem$tracking <- tracking
		ctrl.oem$ioval <- list(iv=list(t1=flsval), ov=list(t1=flsval, t2=flival))
		o.out <- do.call("mpDispatch", ctrl.oem)
		stk0 <- o.out$stk
		idx0 <- o.out$idx
		observations(oem) <- o.out$observations
    tracking <- o.out$tracking

		# --- EST, estimator of stock statistics

		if (!is.null(ctrl.mp$ctrl.est)){
			ctrl.est <- args(ctrl.mp$ctrl.est)
			ctrl.est$method <- method(ctrl.mp$ctrl.est)
			ctrl.est$stk <- stk0
			ctrl.est$idx <- idx0
			ctrl.est$ay <- ay
			ctrl.est$tracking <- tracking
			ctrl.est$ioval <- list(iv=list(t1=flsval, t2=flival), ov=list(t1=flsval))
			out.assess <- do.call("mpDispatch", ctrl.est)
			stk0 <- out.assess$stk
			tracking <- out.assess$tracking
		}
		tracking["F.est",ac(ay)] <- fbar(stk0)[,ac(ay-1)]
		# tracking["B.est",ac(ay)] <- ssb(stk0)[,ac(ay-1)]
	

		# --- PHCR, HCR parameterization
		
    if (!is.null(ctrl.mp$ctrl.phcr)){
			ctrl.phcr <- args(ctrl.mp$ctrl.phcr)
			ctrl.phcr$method <- method(ctrl.mp$ctrl.phcr) 
			ctrl.phcr$stk <- stk0
			ctrl.phcr$ay <- ay
			ctrl.phcr$iy <- iy
			ctrl.phcr$tracking <- tracking
			if(exists("hcrpars")) ctrl.phcr$hcrpars <- hcrpars
			ctrl.phcr$ioval <- list(iv=list(t1=flsval), ov=list(t1=flpval))
			out <- do.call("mpDispatch", ctrl.phcr)
			hcrpars <- out$hcrpars
			tracking <- out$tracking
		}

		# --- HCR
		
    if (!is.null(ctrl.mp$ctrl.hcr)){
			ctrl.hcr <- args(ctrl.mp$ctrl.hcr)
			ctrl.hcr$method <- method(ctrl.mp$ctrl.hcr)
			ctrl.hcr$stk <- stk0
			ctrl.hcr$ay <- ay
			ctrl.hcr$tracking <- tracking
			if(exists("hcrpars")) ctrl.hcr$hcrpars <- hcrpars
			ctrl.hcr$ioval <- list(iv=list(t1=flsval), ov=list(t1=flfval))
			out <- do.call("mpDispatch", ctrl.hcr)
			ctrl <- out$ctrl
			tracking <- out$tracking
		} else {
			ctrl <- getCtrl(yearMeans(fbar(stk0)[,sqy]), "f", ay+1, it)
		}

    # TODO HANDLE multiple targets
  	tracking["metric.hcr", ac(ay)] <- na.exclude(c(ctrl[ac(ay+1),]$value))
		
		# --- IS, implementation system
		
    if (!is.null(ctrl.mp$ctrl.is)){
			ctrl.is <- args(ctrl.mp$ctrl.is)
			ctrl.is$method <- method(ctrl.mp$ctrl.is)
			ctrl.is$ctrl <- ctrl
			ctrl.is$stk <- stk0
			ctrl.is$ay <- ay
			ctrl.is$tracking <- tracking
			ctrl.is$ioval <- list(iv=list(t1=flsval, t2=flfval), ov=list(t1=flfval))
			out <- do.call("mpDispatch", ctrl.is)
			ctrl <- out$ctrl
			tracking <- out$tracking
		  tracking["metric.is", ac(ay)] <- ctrl[ac(ay+1),]$value
		} else {
			tracking["metric.is", ac(ay)] <- tracking["metric.hcr", ac(ay+1)]
		}

		# --- TM, technical measures

		if (!is.null(ctrl.mp$ctrl.tm)){
			ctrl.tm <- args(ctrl.mp$ctrl.tm)
			ctrl.tm$method <- method(ctrl.mp$ctrl.tm)
			ctrl.tm$stk <- stk0
			ctrl.tm$sqy <- sqy
			ctrl.tm$tracking <- tracking
			ctrl.tm$ioval <- list(iv=list(t1=flsval), ov=list(t1=flqval))
			out <- do.call("mpDispatch", ctrl.tm)
			attr(ctrl, "snew") <- out$flq
			tracking <- out$tracking
		}

		# --- IEM, implementation error
		
    if(!missing(iem)){
			ctrl.iem <- args(iem)
			ctrl.iem$method <- method(iem)
			ctrl.iem$ctrl <- ctrl
			ctrl.iem$tracking <- tracking
			ctrl.iem$ioval <- list(iv=list(t1=flfval), ov=list(t1=flfval))
			out <- do.call("mpDispatch", ctrl.iem)
			ctrl <- out$ctrl
			tracking <- out$tracking
		}
  	tracking["metric.iem", ac(ay)] <- na.exclude(c(ctrl[ac(ay+1),]$value))

		# --- FB, fleet dynamics/behaviour

		if (exists(fleetBehaviour(om))){
			ctrl.fb <- args(fleetBehaviour(om))
			ctrl.fb$method <- method(fleetBehaviour(om))
			ctrl.fb$tracking <- tracking
			ctrl.fb$ctrl <- ctrl
			ctrl.fb$ioval <- list(iv=list(t1=flfval), ov=list(t1=flfval))
			out <- do.call("mpDispatch", ctrl.fb)
			ctrl <- out$ctrl
			tracking <- out$tracking
		}
	  
  	tracking["metric.fb", ac(ay)] <- na.exclude(c(ctrl[ac(ay+1),]$value))

    # --- FWD
		
    # NEW selectivity?
    if(!is.null(attr(ctrl, "snew")))
      harvest(stk.om)[, ac(ay+1)] <- attr(ctrl, "snew")
    stk.om <- fwd(stk.om, control=ctrl, sr=sr.om, deviances = sr.om.res,
      effort_max=3)
	}
    if(verbose)
      cat("\n")

  # --- OUTPUT
    res <- as(om, "FLmse")
    stock(res) <- window(stk.om, start=iy, end=fy)
    tracking(res) <- window(tracking, end=fy)
    genArgs(res) <- genArgs
    # TODO accessors
    res@oem <- oem
    res@control <- ctrl.mp

	return(res)
}

