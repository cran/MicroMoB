## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(MicroMoB)
library(ggplot2)
library(data.table)
library(parallel)

## -----------------------------------------------------------------------------
# mosquito parameters
f <- 0.3
q <- 1
eip <- 14
EIP <- eip + 1

lifespan <- 20
g <- 1/lifespan

# human parameters
b <- 0.55
c <- 0.15
r <- 1/200

S <- 1e3
I <- 300
N <- S + I

# transmission parameters
kappa <- (I/N)*c

# equilibrium solutions
Z <- (r*I*N) / (b*f*q*S)
Y <- Z / exp(-g*EIP)
M <- (Z*(g + (f*q*kappa))) / (f*q*kappa*exp(-g*EIP))
lambda <- g*M

## -----------------------------------------------------------------------------
patches <- 1
nstrata <- 1
tmax <- 365 * 2

# human parameters
theta <- diag(nstrata)
H <- N
X <- I

# mosquito parameters
p <- 1 - 1/lifespan
psi <- diag(patches)

# incubating mosquitoes
pp <- p^(eip:1)
pp <- pp / sum(pp)
ZZ <- (Y-Z)*pp
ZZ <- as.matrix(ZZ)

## -----------------------------------------------------------------------------
mod <- make_MicroMoB(tmax = tmax, p = patches)
setup_humans_SIS(mod, stochastic = FALSE, theta = theta, H = H, X = X, b = b, c = c, r = r)
setup_aqua_trace(mod, stochastic = FALSE, lambda = lambda)
setup_mosquito_RM(mod, stochastic = FALSE, f = f, q = q, eip = eip, p = p, psi = psi, M = M, Y = Y, Z = Z)
setup_alternative_trace(mod)
setup_visitor_trace(mod)

mod$mosquito$ZZ <- ZZ

## -----------------------------------------------------------------------------
# matrices to hold output
mosquito_out <- matrix(data = 0, nrow = tmax, ncol = 4, dimnames = list(NULL, c("Day", "M", "Y", "Z")))
mosquito_out[, "Day"] <- 1:tmax
human_out <- matrix(data = 0, nrow = tmax, ncol = 3, dimnames = list(NULL, c("Day", "S", "I")))
human_out[, "Day"] <- 1:tmax

# run it
while (mod$global$tnow <= tmax) {
  compute_bloodmeal(model = mod)
  step_aqua(model = mod)
  step_mosquitoes(model = mod)
  step_humans(model = mod)
  mosquito_out[mod$global$tnow, 2:4] <- c(mod$mosquito$M, mod$mosquito$Y, mod$mosquito$Z)
  human_out[mod$global$tnow, 2:3] <- c(mod$human$H - mod$human$X, mod$human$X)
  mod$global$tnow <- mod$global$tnow + 1L
}

mosquito_out <- as.data.table(mosquito_out)
mosquito_out <- melt(mosquito_out, id.vars = "Day", variable.name = "Compartment", value.name = "Count")
mosquito_out[, "Species" := "Mosquito"]
human_out <- as.data.table(human_out)
human_out <- melt(human_out, id.vars = "Day", variable.name = "Compartment", value.name = "Count")
human_out[, "Species" := "Human"]

det_out <- rbind(mosquito_out, human_out)

## -----------------------------------------------------------------------------
sto_out <- mclapply(X = 1:10, FUN = function(runid) {
  
  mod <- make_MicroMoB(tmax = tmax, p = patches)
  setup_humans_SIS(mod, stochastic = TRUE, theta = theta, H = H, X = X, b = b, c = c, r = r)
  setup_aqua_trace(mod, stochastic = TRUE, lambda = lambda)
  setup_mosquito_RM(mod, stochastic = TRUE, f = f, q = q, eip = eip, p = p, psi = psi, M = M, Y = Y, Z = Z)
  setup_alternative_trace(mod)
  setup_visitor_trace(mod)
  
  # matrices to hold output
  mosquito_out <- matrix(data = 0, nrow = tmax, ncol = 4, dimnames = list(NULL, c("Day", "M", "Y", "Z")))
  mosquito_out[, "Day"] <- 1:tmax
  human_out <- matrix(data = 0, nrow = tmax, ncol = 3, dimnames = list(NULL, c("Day", "S", "I")))
  human_out[, "Day"] <- 1:tmax
  
  # run it
  while (mod$global$tnow <= tmax) {
    compute_bloodmeal(model = mod)
    step_aqua(model = mod)
    step_mosquitoes(model = mod)
    step_humans(model = mod)
    mosquito_out[mod$global$tnow, 2:4] <- c(mod$mosquito$M, mod$mosquito$Y, mod$mosquito$Z)
    human_out[mod$global$tnow, 2:3] <- c(mod$human$H - mod$human$X, mod$human$X)
    mod$global$tnow <- mod$global$tnow + 1L
  }
  
  mosquito_out <- as.data.table(mosquito_out)
  mosquito_out <- melt(mosquito_out, id.vars = "Day", variable.name = "Compartment", value.name = "Count")
  mosquito_out[, "Species" := "Mosquito"]
  human_out <- as.data.table(human_out)
  human_out <- melt(human_out, id.vars = "Day", variable.name = "Compartment", value.name = "Count")
  human_out[, "Species" := "Human"]
  
  out <- rbind(mosquito_out, human_out)
  out[, "Run" := as.integer(runid)]
  
  return(out)
})

sto_out <- do.call(rbind, sto_out)

ggplot(sto_out) +
    geom_line(aes(x = Day, y = Count, color = Compartment, group = interaction(Run, Compartment)), alpha = 0.3) +
    geom_line(data = det_out, aes(x = Day, y = Count, color = Compartment)) +
    facet_wrap(. ~ Species, scales = "free")

