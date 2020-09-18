\name{GKGamma}
\alias{GKGamma}
\alias{print.gkg}
\title{Goodman-Kruskall Gamma}
\description{
Calculates Goodman and Kruskall's Gamma along with its confidence interval
}
\usage{
GKGamma(x, y = NULL, conf.level = NA, ...)
}
\arguments{
  \item{x}{Either a matrix providing a cross-tabulation of two variables (if \code{y} is \code{NULL}) or one variable that will be cross-tabulated with \code{y}.}
  \item{y}{Either \code{NULL} if \code{x} is a matrix or a variable that will be cross-tabulated with \code{x}.}
  \item{conf.level}{The confidence level of the desired confidence interval.  If this is \code{NA} (the default), then no confidence intervall will be produced.}
  \item{\dots}{Other arguments to be passed down to the \code{\link{table}} function.}
}
\value{
A test statistic and optional confidence interval.
}
\examples{
data(alberta)
GKGamma(alberta$d2, alberta$k6group)
}
