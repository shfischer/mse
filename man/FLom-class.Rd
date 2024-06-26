% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FLom-class.R, R/FLmse-class.R
\docType{class}
\name{FLom}
\alias{FLom}
\alias{FLom-class}
\alias{FLom-methods}
\alias{stock,FLom-method}
\alias{stock<-,FLom,FLStock-method}
\alias{sr,FLom-method}
\alias{sr<-,FLom,FLSR-method}
\alias{combine,FLom,FLom-method}
\alias{om<-,FLmse,FLo-method}
\alias{tracking<-,FLmse,FLQuants-method}
\alias{control<-,FLmse,mpCtrl-method}
\alias{oem<-,FLmse,FLoem-method}
\title{A class for an operating model (OM)}
\usage{
FLom(...)

FLom(...)

\S4method{stock}{FLom}(object)

\S4method{stock}{FLom,FLStock}(object) <- value

\S4method{sr}{FLom}(object)

\S4method{sr}{FLom,FLSR}(object) <- value

\S4method{combine}{FLom,FLom}(x, y, ...)

\S4method{om}{FLmse,FLo}(object) <- value

\S4method{tracking}{FLmse,FLQuants}(object) <- value

\S4method{control}{FLmse,mpCtrl}(object) <- value

\S4method{oem}{FLmse,FLoem}(object) <- value
}
\arguments{
\item{...}{additional argument list that might never be used}

\item{object}{object of relevant class (see signature of method)}

\item{value}{Object to assign in slot}
}
\description{
Lorem ipsum dolor sit amet, consectetur adipiscing elit. Pellentesque eleifend
odio ac rutrum luctus. Aenean placerat porttitor commodo. Pellentesque eget porta
libero. Pellentesque molestie mi sed orci feugiat, non mollis enim tristique.
Suspendisse eu sapien vitae arcu lobortis ultrices vitae ac velit. Curabitur id
}
\section{Slots}{

\describe{
\item{\code{stock}}{The population and catch history, \code{FLStock}.}

\item{\code{sr}}{The stock-recruitment relationship, \code{FLSR}.}

\item{\code{refpts}}{The estimated reference points, \code{FLPar}.}

\item{\code{fleetBehaviour}}{Dynamics of the fishing fleet to be used in projections, \code{mseCtrl}.}
}}

\section{Validity}{


\describe{
\item{stock and sr dimensions}{Dimensions 2:6 of the \code{stock} and \code{sr} slots must match.}
\item{rec age}{Stock and stock recruitment residuals must use the recruitment age.}
}
You can inspect the class validity function by using
\code{getValidity(getClassDef('FLom'))}
}

\section{Accessors}{

All slots in the class have accessor and replacement methods defined that
allow retrieving and substituting individual slots.

The values passed for replacement need to be of the class of that slot.
A numeric vector can also be used when replacing FLQuant slots, and the
vector will be used to substitute the values in the slot, but not its other
attributes.
}

\section{Constructor}{

A construction method exists for this class that can take named arguments for
any of its slots. All unspecified slots are then created to match the
requirements of the class validity function.
}

\section{Methods}{

Methods exist for various calculations based on values stored in the class:

\describe{
\item{METHOD}{Neque porro quisquam est qui dolorem ipsum.}
}
}

\examples{
data(sol274)
comb <- combine(iter(om, 1:50), iter(om, 51:100))
all.equal(om, comb)
}
\seealso{
\link{FLComp}
}
\author{
The FLR Team
}
\keyword{classes}
