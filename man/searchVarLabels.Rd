\name{searchVarLabels}
\alias{searchVarLabels}
\alias{searchVarLabels.data.frame}
\alias{searchVarLabels.tbl_df}
\title{Search variable labels}
\description{Searches through variable labels and variable names for a given string. }
\usage{
searchVarLabels(dat, str)
}
\arguments{
  \item{dat}{A \code{data.frame} or \code{tibble} to be searched}
  \item{str}{A string to search for.}
}
\value{
A matrix giving the column number, variable name and variable label of the matching elements.
}
