\name{tTest}
\alias{tTest}
\alias{print.tTest}
\title{Difference of Means Test}
\description{
Produces a more helpful summary of the existing \code{t.test} function.
}
\usage{
tTest(x, y, data, ...)
}
\arguments{
  \item{x}{Character string giving name of grouping variable.}
  \item{y}{Character string giving name of continuous variable. }
  \item{data}{Data where both \code{x} and \code{y} can be found.}
  \item{...}{Other arguments to be passed down to the \code{t.test} function.}
}
\value{
A list of class \code{tTest} with two elements:
\describe{
  \item{sum}{Means, n, and standard errors for both groups and the difference.}
  \item{tt}{The result produced by \code{t.test}.}
}
}
\examples{
data(alberta)
alberta$sex <- rio::factorize(alberta$sex)
tTest("sex", "k7", data=alberta)
}
