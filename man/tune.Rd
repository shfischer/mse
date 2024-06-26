% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run.R
\name{tune}
\alias{tune}
\title{tune}
\usage{
tune(mp, grid, statistics, refpts, ...)
}
\arguments{
\item{mp}{A function executing a projection applying a given MP, see \code{\link{mseBasic}} for an example}

\item{grid}{A name list of \emph{mp} argument values to loop along}

\item{statistics}{A list of performance statistics}

\item{refpts}{The reference points needed to compute the statistics, \emph{FLPar}}

\item{...}{Any other arguments to be passed on to \emph{mp}}
}
\value{
A list or aggregatecd FLR object, depending on the output of \emph{mp}
}
\description{
Carry out multiple runs of an MP for a given dataset over a grid of values
for the MP/HCR paramaters, in order to find the parameter combination(s) that
give the best performance over the chosen statistics, a.k.a. \emph{tuning}
}
\details{
DETAILS
}
\seealso{
\code{\link[mse]{mp}},\code{\link[mse]{performance}}
}
