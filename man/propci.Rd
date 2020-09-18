\name{propci}
\alias{propci}
\title{Approximate and Exact Confidence Intervals for Proportions}
\description{
Produces confidence intervals for proportions through normal approximation and using exact methods based on the Beta distribution. }
\usage{
propci(x, n=NULL, conf.level=.95)
}
\arguments{
  \item{x}{An integer value representing the number of ones (successes) or a variable with only zeros, ones and (potentially) NAs.}
  \item{n}{If \code{x} is an integer, then \code{n} is the total number of observations from which the number of successes comes. Otherwise, \code{n} will be disregarded.}
  \item{conf.level}{Level of confidence for the interval}
}
\value{
A matrix giving the normal approximation and exact confidence intervals.
}
\examples{
data(alberta)
propci(alberta$sex-1)
}
