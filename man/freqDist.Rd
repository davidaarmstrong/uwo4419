\name{freqDist}
\alias{freqDist}
\title{Frequency Distribution}
\description{
Produces a frequency distribution, relative frequency distribution and cumulative distribution from a factor vector. }
\usage{
freqDist(x)
}
\arguments{
  \item{x}{A \code{factor} or a numeric variable with relatively few unique vales.}
}
\value{
A matrix with the following columns:
\describe{
  \item{Freq}{The raw frequency count}
  \item{Percent}{The relative frequency or percentage of observations in each group}
  \item{CuPct}{The cumulative percentage of observations in each group}
}
}
\examples{
data(alberta)
freqDist(alberta$k6group)
}
