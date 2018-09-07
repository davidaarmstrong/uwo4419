\name{plotCIgroup}
\alias{plotCIgroup}
\title{Plot Confidence Intervals by Groups}
\description{Plot confidence intervals by group.}
\usage{
plotCIgroup(form, data, horiz = FALSE, arrowlen = 0, includeOverall = TRUE,
  distr = c("normal", "t"), conflevel = 0.95, las = 2, ...)
}
\arguments{
  \item{form}{Formula giving dependent variable and grouping variable.}
  \item{data}{Data used to calculate confidence intervals}
  \item{horiz}{Logical indicating whether the plotted intervals should be vertical (default) or horizontal (if \code{TRUE}). }
  \item{arrowlen}{Length of whiskers, defaults to zero. }
  \item{includeOverall}{Logical indicating whether the overall mean and confidence interval for all observations should be included.}
  \item{distr}{Distribution to be used for the confidence interval.}
  \item{conflevel}{Confidence level to be used, defaults to 0.95.}
  \item{las}{Graphical parameter to be passed down to the graphic function controling whether labels print parallel or perpendicular to axes.}
  \item{\dots}{Other arguments to be passed down to \code{\link{plot}}. }
}
\value{
A graph
}
