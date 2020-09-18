\name{plotCIgroup}
\alias{plotCIgroup}
\title{Plot Confidence Intervals by Groups}
\description{Plot confidence intervals by group.}
\usage{
plotCIgroup(form, data, includeOverall = TRUE, ...)
}
\arguments{
  \item{form}{Formula giving dependent variable and grouping variable.}
  \item{data}{Data used to calculate confidence intervals}
  \item{includeOverall}{Logical indicating whether the overall mean and confidence interval for all observations should be included.}
  \item{\dots}{Other arguments to be passed down to \code{\link{confidenceInterval}}. }
}
\value{
A ggplot.
}
\examples{
data(alberta)
alberta$k6group <- factor(alberta$k6group,
   labels = attr(attr(alberta$k6group, "labels"), "names"))
plotCIgroup(k3c ~ k6group, data=alberta)
}

