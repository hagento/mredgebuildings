% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcHDDCDD.R
\name{compHDDCDDFactors}
\alias{compHDDCDDFactors}
\title{Calculate HDD/CDD values for different ambient/limit temperature combinations}
\usage{
compHDDCDDFactors(tlow, tup, tlim, tambStd = 5, tlimStd = 5)
}
\arguments{
\item{tlow}{lower temperature boundary}

\item{tup}{upper temperature boundary}

\item{tlim}{named list of limit temperature sequences for \code{HDD} and \code{CDD}}

\item{tambStd}{std of ambient temperature}

\item{tlimStd}{std of limit temperature}
}
\value{
data frame of HDD/CDD values
}
\description{
HDD/CDD values are pre-calculated for an interval \code{tlow}-\code{tup} and
for a set of limit temperatures \code{tlim} with a temperature resolution of
0.1C.
}
\details{
The respective heating/cooling degree days are calculated as the difference
between the assumed ambient and a limit temperature, aggregated to a full day.
The latter defines a threshold above/below which cooling/heating is assumed to
be initiated.

To account for heterogenity in heating/cooling behavior, the ambient and limit
temperature, \code{tamb} and \code{tlim}, are assumed to be normally distributed.
This changes the calculation of a degree day to a double integration of
\code{tLimit - T_ambient_day} with integration boundaries set at 3 standard
deviations, \code{tambStd} and \code{tlimStd}, from \code{tamb} and \code{tlim}
respectively.

As consequence, the ramp function of \code{HDD_day = max(0, tLimit - T_ambient_day)}
changes to a curved function that is above zero even if the mean of \code{T_ambient_day}
is above the mean of \code{tLimit}.
}
