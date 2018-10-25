# tune.R - DESC
# /tune.R

# Copyright European Union, 2018
# Author: Iago Mosqueira (EC JRC) <iago.mosqueira@ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

# BISECTION tuning {{{

#' 
#' @examples

tunebisect <- function(om, control, metrics, indicator, tune=list(ftrg=c(0.20, 1.20)),
  mpargs, prob=0.5, tol=0.01, maxit=12, verbose=TRUE, ...) {
  
  # RUN at min, ...
  cmin <- control
  cmin$ctrl.hcr@args[names(tune)] <- lapply(tune, '[', 1)
  rmin <- mp(om, ctrl.mp=cmin, genArgs=mpargs, scenario=paste0("min"), ...)
  pmin <- performance(metrics(stock(rmin), metrics=metrics), indicator=indicator,
    refpts=refpts(om), probs=NULL, years=mpargs$vy)
  obmin <- mean(pmin$data) - prob

  # CHECK cmin result
  if(isTRUE(all.equal(obmin, 0, tolerance=tol)))
    return(rmin)

  # ... & max
  cmax <- control
  cmax$ctrl.hcr@args[names(tune)] <- lapply(tune, '[', 2)
  rmax <- mp(om, ctrl.mp=cmax, genArgs=mpargs, scenario=paste0("max"), ...)
  pmax <- performance(metrics(stock(rmax), metrics=metrics), indicator=indicator,
    refpts=refpts(om), probs=NULL, years=mpargs$vy)
  obmax <- mean(pmax$data) - prob

  # CHECK cmax result
  if(isTRUE(all.equal(obmax, 0, tolerance=tol)))
    return(rmax)

  # CHECK range includes 0
  if((obmin * obmax) > 0)
    stop("Range of hcr param(s) cannot achieve requested tuning objective probability")

  # LOOP bisecting
  count <- 0
  while(count <= maxit) {

    # RUN at mid
    cmid <- control
    cmid$ctrl.hcr@args[[names(tune)]] <-
      (cmin$ctrl.hcr@args[[names(tune)]] + cmax$ctrl.hcr@args[[names(tune)]]) / 2
    rmid <- mp(om, ctrl.mp=cmid, genArgs=mpargs, scenario=paste0("mid"), ...)
    pmid <- performance(metrics(stock(rmid), metrics=metrics), indicator=indicator,
      refpts=refpts(om), probs=NULL, years=mpargs$vy)
    obmid <- mean(pmid$data) - prob

    if(verbose)
      print(paste0("ob: ", obmid, "; ", names(tune), ": ",
        unlist(cmid$ctrl.hcr@args[names(tune)])))
  
    # CHECK and RETURN cmid result
    if(isTRUE(all.equal(obmid, 0, tolerance=tol))) {
      return(rmid)
    }

    # TEST LEFT
    if((obmin * obmid) < 0) {

      # SET max as new mid
      cmax <- cmid
      obmax <- obmid

    } else {
      
      # SET min as new mid
      cmin <- cmid
      obmin <- obmid
    }

    count <- count + 1
  }

  stop("Solution not found!")

} # }}}
