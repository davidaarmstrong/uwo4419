\name{barplotStats}
\alias{barplotStats}
\title{Barplot of Summary Statistics}
\description{
Produces a barplot of summary statistics.}
\usage{
barplotStats(x, y=NULL, stat="sum", ...)
}
\arguments{
  \item{x}{A \code{factor} or a numeric variable with relatively few unique vales that will be used as categories on the \code{x}-axis of the graph.}
  \item{y}{A numeric variable that will be summarized by each group of \code{x}.  By default, when \code{y} is \code{NULL}, this is a vector of ones that will be summed to give the count of observations in each group of \code{x}.}
  \item{stat}{A string given the statistic that will be used to summarize the data.}
  \item{...}{other arguments that will be passed down to the \code{\link{barplot}} function.}
}
\value{
A plot
}
