# fwd.R - DESC
# /fwd.R

# Copyright European Union, 2018
# Author: Iago Mosqueira (EC JRC) <iago.mosqueira@ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

# getCtrl


getCtrl <- function(values, quant, years, it, rel.years = NA) {

  fwdControl(list(year = years, quant=quant, value = values, rel.year=rel.years))

}
