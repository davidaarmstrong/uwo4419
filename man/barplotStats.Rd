\name{barplotStats}
\alias{barplotStats}
\title{Barplot of Summary Statistics}
\description{
Produces a barplot of summary statistics.}
\usage{
barplotStats(x, y, data, stat="sum", includeN=FALSE, ...)
}
\arguments{
  \item{x}{Character string giving the name of the grouping variable. }
  \item{y}{Character string giving the name of the numeric variable to be summarized. }
  \item{data}{Data frame in which \code{x} and \code{y} can be found. }
  \item{stat}{A string given the statistic that will be used to summarize the data.}
  \item{includeN}{Logical indicating whether the number of observations in the group should be appended to the variable name for plotting purposes.}
  \item{...}{other arguments that will be passed down to \code{stat}.}
}
\value{
A plot
}
\examples{
data(alberta)
barplotStats("k6group", "k3c", alberta, stat="mean")
}
