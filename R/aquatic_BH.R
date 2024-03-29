# a Beverton-Holt style model of aquatic stages

#' @title Setup aquatic (immature) mosquito model with Beverton-Holt dynamics
#' @description A single compartment for all aquatic stages is modeled which
#' suffers density dependent mortality like the Beverton-Holt model.
#' @details All parameters can be passed either as a vector of length equal to `l`, a matrix with `l` rows
#' and `tmax` columns, or a matrix with `l` rows and `365` columns.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param stochastic should the model update deterministically or stochastically?
#' @param molt proportion of immature stages which will mature and emerge as adults each day (may be time and patch varying see [MicroMoB::time_patch_varying_parameter])
#' @param surv daily survival probability (may be time and patch varying see [MicroMoB::time_patch_varying_parameter])
#' @param K carrying capacity (may be time and patch varying see [MicroMoB::time_patch_varying_parameter])
#' @param L initial number of immature mosquitoes
#' @return no return value
#' @export
setup_aqua_BH <- function(model, stochastic, molt, surv, K, L) {
  stopifnot(inherits(model, "MicroMoB"))

  tmax <- model$global$tmax
  l <- model$global$l

  molt_mat <- time_patch_varying_parameter(param = molt, p = l, tmax = tmax)
  surv_mat <- time_patch_varying_parameter(param = surv, p = l, tmax = tmax)
  K_mat <- time_patch_varying_parameter(param = K, p = l, tmax = tmax)

  stopifnot(is.finite(molt_mat))
  stopifnot(molt_mat >= 0)
  stopifnot(is.finite(surv_mat))
  stopifnot(surv_mat >= 0)
  stopifnot(is.finite(K_mat))
  stopifnot(K_mat >= 0)

  stopifnot(length(L) == l)
  A <- rep(0, l)
  eggs <- rep(0, l)

  aqua_class <- c("BH")
  if (stochastic) {
    aqua_class <- c(aqua_class, "BH_stochastic")
    storage.mode(L) <- "integer"
    storage.mode(A) <- "integer"
    storage.mode(eggs) <- "integer"
  } else {
    aqua_class <- c(aqua_class, "BH_deterministic")
  }

  model$aqua <- structure(list(), class = aqua_class)
  model$aqua$molt <- molt_mat
  model$aqua$surv <- surv_mat
  model$aqua$K <- K_mat

  model$aqua$L <- L
  model$aqua$A <- A
  model$aqua$eggs <- eggs

}


#' @title Get parameters for aquatic (immature) model with Beverton-Holt dynamics
#' @description The JSON config file should have two entries:
#'  * stochastic: a boolean value
#'  * molt: a scalar, vector, or matrix (row major)
#'  * surv: a scalar, vector, or matrix (row major)
#'  * K: a scalar, vector, or matrix (row major)
#'  * L: a vector
#'
#' Please see [MicroMoB::time_patch_varying_parameter] for allowed dimensions of entries
#' `molt`, `surv`, and `K`. `L` should be of length equal to the number of patches.
#' For interpretation of the entries, please read [MicroMoB::setup_aqua_BH].
#' @param path a file path to a JSON file
#' @return a named [list]
#' @importFrom jsonlite read_json
#' @examples
#' # to see an example of proper JSON input, run the following
#' library(jsonlite)
#' p <- 5 # number of patches
#' t <- 10 # number of days to simulate
#' par <- list(
#'  "stochastic" = FALSE,
#'  "molt" = 0.3,
#'  "surv" = rep(0.5, 365),
#'  "K" = matrix(rpois(n = t * p, lambda = 100), nrow = p, ncol = t),
#'  "L" = rep(10, p)
#' )
#' toJSON(par, pretty = TRUE)
#' @export
get_config_aqua_BH <- function(path) {
  pars <- read_json(path = file.path(path), simplifyVector = TRUE)

  stopifnot(length(pars) == 5L)
  stopifnot(is.logical(pars$stochastic))
  stopifnot(length(pars$stochastic) == 1L)

  stopifnot(is.numeric(pars$molt))
  stopifnot(is.vector(pars$molt) | is.matrix(pars$molt))

  stopifnot(is.numeric(pars$surv))
  stopifnot(is.vector(pars$surv) | is.matrix(pars$surv))

  stopifnot(is.numeric(pars$K))
  stopifnot(is.vector(pars$K) | is.matrix(pars$K))

  stopifnot(is.numeric(pars$molt))
  stopifnot(is.vector(pars$molt) | is.matrix(pars$molt))

  stopifnot(is.numeric(pars$L))
  return(pars)
}


# output

#' @title Get output for aquatic (immature) mosquito populations with Beverton-Holt dynamics
#' @description Return a [data.frame].
#' @inheritParams output_aqua
#' @return a [data.frame] with columns `L` (immature) and `A` (emerging pupae)
#' @export
output_aqua.BH <- function(model) {
  data.frame(L = model$aqua$L, A = model$aqua$A)
}


# step function

#' @title Update aquatic (immature) mosquito populations for Beverton-Holt dynamics
#' @description This function dispatches on the second class attribute of `model$aqua`
#' for stochastic or deterministic behavior.
#' @inheritParams step_aqua
#' @return no return value
#' @export
step_aqua.BH <- function(model) {
  NextMethod()
}

#' @title Update aquatic (immature) mosquito populations for deterministic Beverton-Holt dynamics
#' @description Run a deterministic state update.
#' @inheritParams step_aqua
#' @return no return value
#' @export
step_aqua.BH_deterministic <- function(model) {

  tnow <- model$global$tnow

  molt <- model$aqua$molt[, tnow]
  surv <- model$aqua$surv[, tnow]
  K <- model$aqua$K[, tnow]
  L <- model$aqua$L

  # get eggs from adult mosquitoes
  eggs <- compute_oviposit(model)

  model$aqua$L = eggs + (1-molt)*surv*L*(K / (L + K))
  model$aqua$A = molt*surv*L*(K / (L + K))
}

#' @title Update aquatic (immature) mosquito populations for stochastic Beverton-Holt dynamics
#' @description Run a stochastic state update.
#' @inheritParams step_aqua
#' @return no return value
#' @importFrom stats rbinom
#' @export
step_aqua.BH_stochastic <- function(model) {

  tnow <- model$global$tnow
  l <- model$global$l

  molt <- model$aqua$molt[, tnow]
  surv <- model$aqua$surv[, tnow]
  K <- model$aqua$K[, tnow]
  L <- model$aqua$L

  # get eggs from adult mosquitoes
  eggs <- compute_oviposit(model)

  # survivors
  survived <- rbinom(n = l, size = L, prob = surv*(K / (L + K)))
  emerging <- rbinom(n = l, size = survived, prob = molt)

  model$aqua$L = eggs + survived - emerging
  model$aqua$A = emerging
}


# get emerging adults

#' @title Compute number of newly emerging adults from Beverton-Holt dynamics
#' @description This function dispatches on the second class attribute of `model$aqua`
#' for stochastic or deterministic behavior.
#' @inheritParams compute_emergents
#' @return a vector of length `l` giving the number of newly emerging adult in each patch
#' @export
compute_emergents.BH <- function(model) {
  model$aqua$A
}
