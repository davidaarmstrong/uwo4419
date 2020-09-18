\name{histDiscrete}
\alias{histDiscrete}
\title{Create a histogram from a discrete variable}
\description{A histogram for a discrete variable is basically a bar chart with no space between the bars.}
\usage{
histDiscrete(x, data, ...)
}
\arguments{
  \item{x}{A character string giving the name of the variable to be plotted.}
  \item{data}{A data frame in which \code{x} can be found. }
  \item{...}{Not implemented}
}
\value{
A \code{ggplot}.
}
\examples{
data(alberta)
alberta$k3c <- ifelse(alberta$k3c > 11, NA, alberta$k3c)
histDiscrete("k3c", alberta)
}


