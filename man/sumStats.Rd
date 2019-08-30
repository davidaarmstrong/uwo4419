\name{sumStats}
\alias{sumStats}
\title{Summary Statistics}
\description{
Producessummary statistics for a set of variables.}
\usage{
sumStats(data, vars, byvar=NULL, convertFactors=FALSE)
}
\arguments{
  \item{data}{A \code{data.frame} or \code{tibble}.}
  \item{vars}{A vector of character strings giving the names of variables to be summarized}
  \item{byvar}{A character string giving the name of a grouping variable in \code{data}.}
  \item{convertFactors}{Logical indicating whether factors should be converted to numberic and then summarized.}
}
\value{
A matrix giving summary statistics
}
