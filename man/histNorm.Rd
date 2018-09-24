\name{histNorm}
\alias{histNorm}
\title{Histogram with Density Curves Super-imposed}
\description{Histogram with density curves super-imposed}
\usage{
histNorm(data, normCurv=TRUE, densCurve=FALSE, ...)
}
\arguments{
  \item{data}{Values to be plotted}
  \item{normCurve}{Whether or not to impose a normal curve with same mean and SD as plotted variable (solid line). }
  \item{densCurve}{Whether or not to impose a smoothe density estimate of the plotted variable (dashed line). }
  \item{\dots}{Other arguments to be passed down to \code{\link{hist}}. }
}
\value{
A graph
}
