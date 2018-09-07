\name{pmc}
\alias{pmc}
\title{Plot of Pairwise Correlations}
\description{
Produces a plot of pairwise correlations of terms from the model matrix of a \code{lm} or \code{glm}.}
\usage{
pmc(obj, col = RColorBrewer::brewer.pal(6, "RdBu"), ...)
}
\arguments{
  \item{obj}{A model of class \code{glm} or \code{lm}. }
  \item{col}{A vector of colors, defaults to the \code{RdBu} palette from the \pkg{RColorBrewer} package.}
  \item{\dots}{Other arguments to be passed down to \code{\link{levelplot}}.}
}
\value{
A lattice graph.
}
