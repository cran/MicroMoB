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

lifespan <- 20
g <- 1/lifespan
p <- 1- g

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
Y <- Z / (p^eip)
M <- (Z*(g + (f*q*p*kappa))) / (f*q*p*kappa*(p^eip))
lambda <- g*M

## -----------------------------------------------------------------------------
patches <- 1
nstrata <- 1
tmax <- 365 * 2
theta <- diag(nstrata)
psi <- diag(patches)

## -----------------------------------------------------------------------------
mod <- make_MicroMoB(tmax = tmax, p = patches)
setup_humans_SIS(mod, stochastic = FALSE, theta = theta, H = N, X = I, b = b, c = c, r = r)
setup_aqua_trace(mod, stochastic = FALSE, lambda = lambda)
setup_mosquito_RM(mod, stochastic = FALSE, f = f, q = q, eip = eip, p = p, psi = psi, M = M, Y = Y, Z = Z)

## -----------------------------------------------------------------------------
mosy_out <- data.table::CJ(day = 1:tmax, state = c('M', 'Y', 'Z'), value = NaN, species = "mosquito")
mosy_out <- mosy_out[c('M', 'Y', 'Z'), on="state"]
data.table::setkey(mosy_out, day)

human_out <- data.table::CJ(day = 1:tmax, state = c('S', 'I'), value = NaN, species = "human")
human_out <- human_out[c('S', 'I'), on="state"]
data.table::setkey(human_out, day)

while (get_tnow(mod) <= tmax) {
  
  compute_bloodmeal_simple(model = mod)
  step_aqua(model = mod)
  step_mosquitoes(model = mod)
  step_humans(model = mod)
  
  mosy_out[day == get_tnow(mod) & state == 'M', value := mod$mosquito$M]
  mosy_out[day == get_tnow(mod) & state == 'Y', value := mod$mosquito$Y]
  mosy_out[day == get_tnow(mod) & state == 'Z', value := mod$mosquito$Z]
  
  human_out[day == get_tnow(mod) & state == 'S', value := mod$human$H - mod$human$X]
  human_out[day == get_tnow(mod) & state == 'I', value := mod$human$X]

  mod$global$tnow <- mod$global$tnow + 1L
}

det_out <- rbind(mosy_out, human_out)

## -----------------------------------------------------------------------------
sto_out <- mclapply(X = 1:10, FUN = function(runid) {
  
  mod <- make_MicroMoB(tmax = tmax, p = patches)
  setup_humans_SIS(mod, stochastic = TRUE, theta = theta, H = N, X = I, b = b, c = c, r = r)
  setup_aqua_trace(mod, stochastic = TRUE, lambda = lambda)
  setup_mosquito_RM(mod, stochastic = TRUE, f = f, q = q, eip = eip, p = p, psi = psi, M = M, Y = Y, Z = Z)
  
  mosy_out <- data.table::CJ(day = 1:tmax, state = c('M', 'Y', 'Z'), value = NaN, species = "mosquito")
  mosy_out <- mosy_out[c('M', 'Y', 'Z'), on="state"]
  data.table::setkey(mosy_out, day)
  
  human_out <- data.table::CJ(day = 1:tmax, state = c('S', 'I'), value = NaN, species = "human")
  human_out <- human_out[c('S', 'I'), on="state"]
  data.table::setkey(human_out, day)
    
  while (get_tnow(mod) <= tmax) {
    
    compute_bloodmeal_simple(model = mod)
    step_aqua(model = mod)
    step_mosquitoes(model = mod)
    step_humans(model = mod)
    
    mosy_out[day == get_tnow(mod) & state == 'M', value := mod$mosquito$M]
    mosy_out[day == get_tnow(mod) & state == 'Y', value := mod$mosquito$Y]
    mosy_out[day == get_tnow(mod) & state == 'Z', value := mod$mosquito$Z]
    
    human_out[day == get_tnow(mod) & state == 'S', value := mod$human$H - mod$human$X]
    human_out[day == get_tnow(mod) & state == 'I', value := mod$human$X]
  
    mod$global$tnow <- mod$global$tnow + 1L
  }
  
  out <- rbind(mosy_out, human_out)
  out[, "run" := as.integer(runid)]
  
  return(out)
})

sto_out <- data.table::rbindlist(sto_out)

ggplot(sto_out) +
    geom_line(aes(x = day, y = value, color = state, group = run), alpha = 0.3) +
    geom_line(data = det_out, aes(x = day, y = value, color = state)) +
    facet_wrap(species ~ state, scales = "free")

