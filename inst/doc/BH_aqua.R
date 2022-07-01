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
patches <- 1
tmax <- 1e2

M <- 120
p <- 0.9
lambda <- M*(1-p)

nu <- 25
f <- 0.3
eggs <- nu * f * M

# static pars
molt <-  0.1
surv <- 0.9

# solve L
L <- lambda * ((1/molt) - 1) + eggs
K <- - (lambda * L) / (lambda - L*molt*surv)

## -----------------------------------------------------------------------------
# deterministic run
mod <- make_MicroMoB(tmax = tmax, p = patches)
setup_aqua_BH(model = mod, stochastic = FALSE, molt = molt, surv = surv, K = K, L = L)
setup_mosquito_RM(model = mod, stochastic = FALSE, f = f, q = 0.9, eip = 10, p = p, psi = diag(1), nu = nu, M = M, Y = 0, Z = 0)

out_det <- data.table::CJ(day = 1:tmax, state = c('L', 'A', 'M'), value = NaN)
out_det <- out_det[c('L', 'A', 'M'), on="state"]
data.table::setkey(out_det, day)

while (get_tnow(mod) <= tmax) {
  
  step_aqua(model = mod)
  step_mosquitoes(model = mod)
  
  out_det[day == get_tnow(mod) & state == 'L', value := mod$aqua$L]
  out_det[day == get_tnow(mod) & state == 'A', value := mod$aqua$A]
  out_det[day == get_tnow(mod) & state == 'M', value := mod$mosquito$M]

  mod$global$tnow <- mod$global$tnow + 1L
}

## -----------------------------------------------------------------------------
# stochastic runs
out_sto <- mclapply(X = 1:10, FUN = function(runid) {
  
  mod <- make_MicroMoB(tmax = tmax, p = patches)
  setup_aqua_BH(model = mod, stochastic = TRUE, molt = molt, surv = surv, K = K, L = L)
  setup_mosquito_RM(model = mod, stochastic = TRUE, f = f, q = 0.9, eip = 10, p = p, psi = diag(1), nu = nu, M = M, Y = 0, Z = 0)
  
  out <- data.table::CJ(day = 1:tmax, state = c('L', 'A', 'M'), value = NaN)
  out <- out[c('L', 'A', 'M'), on="state"]
  data.table::setkey(out, day)
  
  while (get_tnow(mod) <= tmax) {
    
    step_aqua(model = mod)
    step_mosquitoes(model = mod)
    
    out[day == get_tnow(mod) & state == 'L', value := mod$aqua$L]
    out[day == get_tnow(mod) & state == 'A', value := mod$aqua$A]
    out[day == get_tnow(mod) & state == 'M', value := mod$mosquito$M]
  
    mod$global$tnow <- mod$global$tnow + 1L
  }
  
  out[, 'run' := as.integer(runid)] 
  return(out)
})

out_sto <- data.table::rbindlist(out_sto)

## -----------------------------------------------------------------------------
ggplot(out_sto) +
    geom_line(aes(x = day, y = value, color = state, group = run), alpha = 0.35) +
    geom_line(data = out_det, aes(x = day, y = value, color = state)) +
    facet_wrap(. ~ state, scales = "free")

