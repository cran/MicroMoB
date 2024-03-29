% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mosquito_RM.R
\name{output_mosquitoes.RM}
\alias{output_mosquitoes.RM}
\title{Get output for Ross-Macdonald mosquito populations}
\usage{
\method{output_mosquitoes}{RM}(model)
}
\arguments{
\item{model}{an object from \link{make_MicroMoB}}
}
\value{
a \link{data.frame} with columns \code{M} (all adult mosquitoes), \code{Y} (infected mosquitoes), and \code{Z} (infectious mosquitoes), and rows
correspond to places.
}
\description{
Return a \link{data.frame}.
}
