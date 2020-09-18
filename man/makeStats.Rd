\name{makeStats}
\alias{makeStats}
\alias{lambda}
\alias{phi}
\alias{V}
\alias{ord.somers.d}
\alias{ord.gamma}
\alias{tau.b}
\alias{simtable}
\alias{simrho}
\alias{print.ktb}
\alias{concordant}
\alias{discordant}
\title{Measures of Association}
\description{
Calculate measures of association for a pair of \code{factor} variables. }
\usage{
makeStats(x, y, chisq = FALSE, phi = FALSE, cramersV = FALSE, lambda = FALSE,
  gamma = FALSE, d = FALSE, taub = FALSE, rho = FALSE, n = 1000)
}
\arguments{
  \item{x}{Either a matrix providing a cross-tabulation of two variables (if \code{y} is \code{NULL}) or one variable that will be cross-tabulated with \code{y}.}
  \item{y}{Either \code{NULL} if \code{x} is a matrix or a variable that will be cross-tabulated with \code{x}.}
  \item{chisq}{Logical indicating whether Chi-squared should be calculated.}
  \item{phi}{Logical indicating whether phi should be calculated.}
  \item{cramersV}{Logical indicating whether Cramer's Vshould be calculated.}
  \item{lambda}{Logical indicating whether Lambda should be calculated.}
  \item{gamma}{Logical indicating whether Gamma should be calculated.}
  \item{d}{Logical indicating whether Somer's D should be calculated.}
  \item{taub}{Logical indicating whether Kendall's Tau b should be calculated.}
  \item{rho}{Logical indicating whether Spearman's Rho should be calculated.}
  \item{n}{Number of simulated values used to calculate the p-value}
}
\value{
  A matrix of statistics and simulated p-values.
}
\examples{
data(alberta)
makeStats(alberta$k6group, alberta$d2, chisq=TRUE, taub=TRUE)
}
