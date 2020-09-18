\name{confidenceInterval}
\alias{confidenceInterval}
\title{Confidence Interval}
\description{Make a confidence interval for the mean of a vector of values.}
\usage{
confidenceInterval(x, confidence = 0.95, na.rm = TRUE, distr = c("normal", "t"))
}
\arguments{
  \item{x}{A vector of values for which the mean and its confidence interval will be calculated.}
  \item{confidence}{Desired level of confidence for the confidence interval.  Defaults to 0.95.}
  \item{na.rm}{Logical indicating whether missing observations should be removed.  Defaults to \code{TRUE}.}
  \item{distr}{Distribution used to find the critical value.  Default is \code{"normal"}.}
}
\value{
A vector of values giving the mean, confidence interval and standard error.}
\examples{
data(alberta)
confidenceInterval(alberta$k7)
}
