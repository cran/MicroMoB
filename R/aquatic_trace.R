# classes and methods to deal with immature mosquitoes and oviposition

#' @title Setup aquatic (immature) mosquito model with trace (forced) emergence
#' @description Emergence is passed as a (possibly time varying) parameter which is
#' decoupled from the adult mosquito dynamics.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param lambda daily emergence of mosquitoes, may be time and patch varying, see [MicroMoB::time_patch_varying_parameter]
#' @param stochastic should the model update deterministically or stochastically?
#' @return no return value
#' @export
setup_aqua_trace <- function(model, lambda, stochastic) {
  stopifnot(inherits(model, "MicroMoB"))

  tmax <- model$global$tmax
  p <- model$global$p

  stopifnot(is.finite(lambda))
  stopifnot(lambda >= 0)

  lambda_mat <- time_patch_varying_parameter(param = lambda, p = p, tmax = tmax)

  aqua_class <- c("trace")
  if (stochastic) {
    aqua_class <- c(aqua_class, "trace_stochastic")
  } else {
    aqua_class <- c(aqua_class, "trace_deterministic")
  }

  model$aqua <- structure(list(), class = aqua_class)
  model$aqua$lambda <- lambda_mat

}


# step function

#' @title Update aquatic (immature) mosquito populations for forced emergence
#' @description This function does nothing as trace models are do not have
#' endogenous dynamics.
#' @inheritParams step_aqua
#' @return no return value
#' @export
step_aqua.trace <- function(model) {invisible()}


# get emerging adults

#' @title Compute number of newly emerging adults from forcing term
#' @description This function dispatches on the second class attribute of `model$aqua`
#' for stochastic or deterministic behavior.
#' @inheritParams compute_emergents
#' @details see [MicroMoB::compute_emergents.trace_deterministic] and [MicroMoB::compute_emergents.trace_stochastic]
#' @return no return value
#' @export
compute_emergents.trace <- function(model) {
  NextMethod()
}

#' @title Compute number of newly emerging adults from forcing term (deterministic)
#' @description Return the column of the lambda matrix for this day.
#' @inheritParams compute_emergents
#' @return a vector of length `p` giving the number of newly emerging adult in each patch
#' @export
compute_emergents.trace_deterministic <- function(model) {
  return(model$aqua$lambda[, model$global$tnow])
}

#' @title Compute number of newly emerging adults from forcing term (stochastic)
#' @description Draw a Poisson distributed number of emerging adults with mean parameter
#' from the column of the trace matrix for this day.
#' @inheritParams compute_emergents
#' @return a vector of length `p` giving the number of newly emerging adult in each patch
#' @importFrom stats rpois
#' @export
compute_emergents.trace_stochastic <- function(model) {
  return(rpois(n = model$global$p, lambda = model$aqua$lambda[, model$global$tnow]))
}