% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visitor_trace.R
\name{setup_visitor_trace}
\alias{setup_visitor_trace}
\title{Setup trace driven visitors}
\usage{
setup_visitor_trace(model, Wd = NULL, xd = NULL)
}
\arguments{
\item{model}{an object from \link{make_MicroMoB}}

\item{Wd}{a time varying trace of visitor host availability passed to \link{time_patch_varying_parameter}
or \code{NULL} to set to \code{0} (no visitors)}

\item{xd}{a time varying trace of visitor net infectiousness passed to \link{time_patch_varying_parameter}
or \code{NULL} to set to \code{0} (no visitors)}
}
\value{
no return value
}
\description{
This model complies with the visitors component interface. It adds
a named list \code{model$visitor}.
}
