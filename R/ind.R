# ind.R - DESC
# /ind.R

# Copyright European Union, 2018
# Author: Iago Mosqueira (EC JRC) <iago.mosqueira@ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

# cpue.ind {{{

cpue.ind <- function(stk, idx, nyears=5, ay, tracking) {

  # INDEX slot
  idx <- index(idx[[1]])

  # SUBSET last nyears
  ind <- idx[, seq(dim(idx)[2] - nyears, dim(idx)[2])]

  # SLOPE by iter
  dat <- data.table(as.data.frame(ind))
  slope <- dat[, .(slope=coef(lm(log(data)~year))[2]), by=iter]
	
  # OUTPUT
  tracking["cpue.est", ac(ay)] <- c(slope[,slope])

  # ADD ind to tracking$B.est
  list(stk=stk, tracking=tracking)
} # }}}
