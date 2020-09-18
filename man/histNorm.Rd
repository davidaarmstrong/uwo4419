\name{histNorm}
\alias{histNorm}
\title{Histogram with Density Curves Super-imposed}
\description{Histogram with density curves super-imposed}
\usage{
histNorm(x, data, normCurve=TRUE, densCurve=FALSE, bins=30)
}
\arguments{
  \item{x}{Character string giving the name of the variable to be plotted}
  \item{data}{A data frame in which \code{x} can be found.}
  \item{normCurve}{Whether or not to impose a normal curve with same mean and SD as plotted variable (solid line). }
  \item{densCurve}{Whether or not to impose a smoothe density estimate of the plotted variable (dashed line). }
  \item{bins}{Number of bins to use in the histogram}
}
\value{
A \code{ggplot}
}
\examples{
data(alberta)
alberta$age <- ifelse(alberta$age > 94, NA, alberta$age)
histNorm("age", alberta)
}


