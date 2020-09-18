\name{plotStdRes}
\alias{plotStdRes}
\title{Plot Standardized Residuals from a Cross-tabulation}
\description{
Uses \code{ggplot2} to visualize the standardized residuals from a cross-tabulation}
\usage{
plotStdRes(x)
}
\arguments{
  \item{x}{A cross-tabulation of two variables}
}
\value{
A ggplot
}
\examples{
alberta$k6group <- rio::factorize(alberta$k6group)
alberta$d2 <- rio::factorize(alberta$d2)
tab <- with(alberta, table(k6group, d2))
plotStdRes(tab)
}


