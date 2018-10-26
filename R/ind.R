# ind.R - DESC
# /ind.R

# Copyright European Union, 2018
# Author: Iago Mosqueira (EC JRC) <iago.mosqueira@ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

# cpue.ind {{{

cpue.ind <- function(stk, idx, nyears=5, ayears=3, ay, tracking) {
  
  # INDEX slot
  idx <- index(idx[[1]])

  # SUBSET last nyears
  ind <- idx[, seq(dim(idx)[2] - nyears, dim(idx)[2])]

  # SLOPE by iter
  dat <- data.table(as.data.frame(ind))
  slope <- dat[, .(slope=coef(lm(log(data)~year))[2]), by=iter]

  # WEIGHTED average index of last ayears
  mind <- yearSums(tail(ind, ayears) * 
   c(0.50 * seq(1, ayears - 1) / sum(seq(1, ayears - 1)), 0.50))
                   
  # OUTPUT
  tracking["cpue.slope", ac(ay)] <- c(slope[,slope])
  tracking["cpue.mean", ac(ay)] <- mind

  list(stk=stk, tracking=tracking)
} # }}}

# perfectB.ind {{{

perfectB.ind <- function(stk, idx, nyears=5, ay, tracking) {

  # INDEX slot
  idx <- index(idx[[1]])

  # SUBSET last nyears
  ind <- idx[, seq(dim(idx)[2] - nyears, dim(idx)[2])]

  # SLOPE by iter
  dat <- data.table(as.data.frame(ind))
  slope <- dat[, .(slope=coef(lm(log(data)~year))[2]), by=iter]
	
  # OUTPUT
  tracking["cpue.ind", ac(ay)] <- c(slope[,slope])

  # ADD ind to tracking$B.est
  list(stk=stk, tracking=tracking)
} # }}}
