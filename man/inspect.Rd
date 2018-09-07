\name{inspect}
\alias{inspect}
\alias{inspect.data.frame}
\alias{inspect.tbl_df}
\title{
  Inspect a variable in a dataset.
}
\description{Inspect a variable in a dataset. }
\usage{
inspect(data, x, ...)
}
\arguments{
  \item{data}{A \code{data.frame} or \code{tibble} that contains the variable to be inspected.}
  \item{x}{A character string giving the name of the variable to be inspected}
  \item{\dots}{Currently un-implemented.}
}
\value{
A list with three elements:
\describe{
  \item{variable.label}{The variables descriptive label (if available).}
  \item{value.labels}{The labels of the variables numeric values (if available).}
  \item{freq_dist}{Frequency distribution and relative frequency distribution of responses.}
}
}
