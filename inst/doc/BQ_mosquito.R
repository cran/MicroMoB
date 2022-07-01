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
p <- l <- 1
tmax <- 1e2

M <- 120
pB <- 0.8
pQ <- 0.95
PsiB <- 0.5
PsiQ <- 0.85

B <- (M - (M*pQ*(1-PsiQ))) / ((pB*PsiB) - (pQ*(1-PsiQ)) + 1)
Q <- (M*pB*PsiB) / ((pB*PsiB) - (pQ*(1-PsiQ)) + 1)

lambda <- B - (pB*(1-PsiB)*B) - (pQ*PsiQ*Q)

nu <- 25
eggs <- nu * PsiQ * Q

# static pars
molt <-  0.1
surv <- 0.9

# solve L
L <- lambda * ((1/molt) - 1) + eggs
K <- - (lambda * L) / (lambda - L*molt*surv)

## -----------------------------------------------------------------------------
# deterministic run
mod <- make_MicroMoB(tmax = tmax, p = p, l = l)
setup_aqua_BH(model = mod, stochastic = FALSE, molt = molt, surv = surv, K = K, L = L)
setup_mosquito_BQ(model = mod, stochastic = FALSE, eip = 5, pB = pB, pQ = pQ, psiQ = PsiQ, Psi_bb = matrix(1), Psi_bq = matrix(1), Psi_qb = matrix(1), Psi_qq = matrix(1), nu = nu, M = c(B, Q), Y = matrix(0, nrow = 2, ncol = 6))

out_det <- data.table::CJ(day = 1:tmax, state = c('L', 'A', 'B', 'Q'), value = NaN)
out_det <- out_det[c('L', 'A', 'B', 'Q'), on="state"]
data.table::setkey(out_det, day)

mod$mosquito$q <- 0.3
mod$mosquito$f <- log(1 - PsiB) / -0.3

while (get_tnow(mod) <= tmax) {
  
  step_aqua(model = mod)
  step_mosquitoes(model = mod)
  
  out_det[day == get_tnow(mod) & state == 'L', value := mod$aqua$L]
  out_det[day == get_tnow(mod) & state == 'A', value := mod$aqua$A]
  out_det[day == get_tnow(mod) & state == 'B', value := mod$mosquito$M[1]]
  out_det[day == get_tnow(mod) & state == 'Q', value := mod$mosquito$M[2]]
  
  mod$global$tnow <- mod$global$tnow + 1L
}

## -----------------------------------------------------------------------------
# stochastic runs
out_sto <- mclapply(X = 1:10, FUN = function(runid) {
  
  mod <- make_MicroMoB(tmax = tmax, p = p, l = l)
  setup_aqua_BH(model = mod, stochastic = TRUE, molt = molt, surv = surv, K = K, L = L)
  setup_mosquito_BQ(model = mod, stochastic = TRUE, eip = 5, pB = pB, pQ = pQ, psiQ = PsiQ, Psi_bb = matrix(1), Psi_bq = matrix(1), Psi_qb = matrix(1), Psi_qq = matrix(1), nu = nu, M = c(B, Q), Y = matrix(0, nrow = 2, ncol = 6))
  
  out <- data.table::CJ(day = 1:tmax, state = c('L', 'A', 'B', 'Q'), value = NaN)
  out <- out[c('L', 'A', 'B', 'Q'), on="state"]
  data.table::setkey(out, day)

  mod$mosquito$q <- 0.3
  mod$mosquito$f <- log(1 - PsiB) / -0.3
  
  while (get_tnow(mod) <= tmax) {
    step_aqua(model = mod)
    step_mosquitoes(model = mod)
    
    out[day == get_tnow(mod) & state == 'L', value := mod$aqua$L]
    out[day == get_tnow(mod) & state == 'A', value := mod$aqua$A]
    out[day == get_tnow(mod) & state == 'B', value := mod$mosquito$M[1]]
    out[day == get_tnow(mod) & state == 'Q', value := mod$mosquito$M[2]]
    
    mod$global$tnow <- mod$global$tnow + 1L
  }
  
  out[, 'run' := as.integer(runid)]
  return(out)
})

## -----------------------------------------------------------------------------
out_sto <- data.table::rbindlist(out_sto)

ggplot(data = out_sto) +
  geom_line(aes(x = day, y = value, color = state, group = run), alpha = 0.35) +
  geom_line(data = out_det, aes(x = day, y = value, color = state)) +
  facet_wrap(. ~ state, scales = "free")

