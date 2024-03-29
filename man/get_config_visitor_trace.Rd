% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visitor_trace.R
\name{get_config_visitor_trace}
\alias{get_config_visitor_trace}
\title{Get parameters for trace driven visitors}
\usage{
get_config_visitor_trace(path)
}
\arguments{
\item{path}{a file path to a JSON file}
}
\value{
a named \link{list}
}
\description{
The JSON config file should have two entries:
\itemize{
\item Wd: vector or matrix (see \link{time_patch_varying_parameter} for valid dimensions)
\item xd: vector or matrix (see \link{time_patch_varying_parameter} for valid dimensions)
}

For interpretation of the entries, please read \link{setup_visitor_trace}.
}
\examples{
# to see an example of proper JSON input, run the following
library(jsonlite)
par <- list(
 "Wd" = rep(1, 5),
 "xd" = rep(0.01, 365)
)
toJSON(par, pretty = TRUE)
}
