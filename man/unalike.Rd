\name{unalike}
\alias{unalike}
\title{Unalikability}
\description{Unalikability is a measure of variation for categorical variables.  It essentially measures the probability of drawing two non-equal values at random. }
\usage{
unalike(x)
}
\arguments{
  \item{x}{A factor or other variable with relatively few levels}
}
\value{
The unalikability score.
}
\references{
Gary D. Kader and Mike Perry (2007)  \sQuote{Variability for Categorical Variables} Journal of Statistics Education Volume 15, Number 2. \url{http://jse.amstat.org/v15n2/kader.html}
}
\examples{
data(Alberta)
unalike(alberta$k6group)
}
