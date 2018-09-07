\name{pwCorrMat}
\alias{pwCorrMat}
\alias{print.pwc}
\alias{sig.cor}
\title{Calculate Pairwise Correlations and Significance}
\description{
Calculates pairwise Pearson correlations and their significance through one of three ways.}
\usage{
pwCorrMat(X, method = c("z", "t", "sim"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A matrix of variables whose correlation is to be calculated.}
  \item{method}{Method for calculating significance. The \code{z} and \code{t} options use an approximation where \code{sim} uses Monte Carlo simulation.}
  \item{\dots}{Other arguments to be passed down to \code{sig.cor}.  In particular, these would be \code{n.sim}, the number of simulations used and \code{two.sided} a logical indicating whether the significance calculation should be one- or two-sided.  You could also use this to pass other arguments down to the \code{\link{cor}} function.}
}
\value{
A list with three elements:
\describe{
  \item{rSig}{A matrix indicating which pairwise correlations are significant}
  \item{r}{A matrix giving the pairwise correlations}
  \item{p}{A matrix giving the p-values for the pairwise correlations}
}
}
