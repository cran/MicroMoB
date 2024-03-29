% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mosquito_interface.R
\name{compute_Z}
\alias{compute_Z}
\title{Compute density of infective mosquitoes (\eqn{Z})}
\usage{
compute_Z(model)
}
\arguments{
\item{model}{an object from \link{make_MicroMoB}}
}
\value{
a vector of length \code{p} giving the density of infected and infectious mosquitoes in each patch
}
\description{
This method dispatches on the type of \code{model$mosquito}. \eqn{Z}
is also known as the "sporozoite rate" in malariology.
}
