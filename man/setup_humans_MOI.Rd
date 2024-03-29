% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/humans_MOI.R
\name{setup_humans_MOI}
\alias{setup_humans_MOI}
\title{Setup humans with MOI (multiplicity of infection) pathogen model}
\usage{
setup_humans_MOI(
  model,
  stochastic,
  theta,
  wf = NULL,
  H,
  MOI,
  b = 0.55,
  c = 0.15,
  r = 1/200,
  sigma = 1
)
}
\arguments{
\item{model}{an object from \link{make_MicroMoB}}

\item{stochastic}{should the model update deterministically or stochastically?}

\item{theta}{a time spent matrix}

\item{wf}{biting weights}

\item{H}{vector of strata population sizes}

\item{MOI}{a matrix giving the distribution of persons across strata (columns) and
multiplicity of infection (rows).}

\item{b}{transmission efficiency (mosquito to human)}

\item{c}{transmission efficiency (human to mosquito)}

\item{r}{recovery rate (inverse of infectious duration)}

\item{sigma}{control non-independence of pathogen clearance; \code{sigma > 1} indicates competition
(clearance is faster than independent) and \code{sigma < 1} indicates facilitation (clearance is slower than independent).}
}
\value{
no return value
}
\description{
This is a queueing model (M/M/inf) of superinfection in humans.
}
\note{
The \link{step_humans} method for the MOI model will grow the \code{MOI}
matrix (add rows) if an individual's MOI exceeds the size of the matrix; therefore
it's a good idea to pad the input matrix with extra empty rows to avoid
reallocating memory during the simulation as much as possible.
}
