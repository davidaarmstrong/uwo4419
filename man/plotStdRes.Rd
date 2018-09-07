\name{plotStdRes}
\alias{plotStdRes}
\title{Plot Standardized Residuals from a Cross-tabulation}
\description{
Uses \code{levelplot} to visualize the standardized residuals from a cross-tabulation}
\usage{
plotStdRes(x, col = RColorBrewer::brewer.pal(10, "RdBu"), ...)
}
\arguments{
  \item{x}{A cross-tabulation of two variables}
  \item{col}{A vector of colrs.  Defaults to the \code{RdBu} palette from the \pkg{RColorBrewer} package.}
  \item{\dots}{Other arguments to be passed down to \code{\link{levelplot}}.}
}
\value{
A lattice graph
}
