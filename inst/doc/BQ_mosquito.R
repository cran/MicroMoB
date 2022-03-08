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

out_det <- matrix(data = NaN, nrow = tmax + 1, ncol = 4)
out_det[1L, ] <- c(mod$aqua$L, mod$aqua$A, mod$mosquito$M)

mod$mosquito$q <- 0.3
mod$mosquito$f <- log(1 - PsiB) / -0.3

while (mod$global$tnow <= tmax) {
  step_aqua(model = mod)
  step_mosquitoes(model = mod)
  out_det[mod$global$tnow + 1L, ] <- c(mod$aqua$L, mod$aqua$A, mod$mosquito$M)
  mod$global$tnow <- mod$global$tnow + 1L
}

## -----------------------------------------------------------------------------
# stochastic runs
out_sto <- mclapply(X = 1:10, FUN = function(runid) {
  
  mod <- make_MicroMoB(tmax = tmax, p = p, l = l)
  setup_aqua_BH(model = mod, stochastic = TRUE, molt = molt, surv = surv, K = K, L = L)
  setup_mosquito_BQ(model = mod, stochastic = TRUE, eip = 5, pB = pB, pQ = pQ, psiQ = PsiQ, Psi_bb = matrix(1), Psi_bq = matrix(1), Psi_qb = matrix(1), Psi_qq = matrix(1), nu = nu, M = c(B, Q), Y = matrix(0, nrow = 2, ncol = 6))
  
  out_run <- matrix(data = NaN, nrow = tmax + 1, ncol = 4)
  out_run[1L, ] <- c(mod$aqua$L, mod$aqua$A, mod$mosquito$M)
  
  mod$mosquito$q <- 0.3
  mod$mosquito$f <- log(1 - PsiB) / -0.3
  
  while (mod$global$tnow <= tmax) {
    step_aqua(model = mod)
    step_mosquitoes(model = mod)
    out_run[mod$global$tnow + 1L, ] <- c(mod$aqua$L, mod$aqua$A, mod$mosquito$M)
    mod$global$tnow <- mod$global$tnow + 1L
  }
  
  out_run <- as.data.frame(out_run)
  out_run$run <- as.integer(runid)
  return(out_run)
})

## -----------------------------------------------------------------------------
out_det <- as.data.table(out_det)
out_det[, "Day" := 0:tmax]
out_det <- melt(out_det, id.vars = "Day", variable.name = "Stage", value.name = "Count")
levels(out_det$Stage) <- c("L", "A", "B", "Q")

out_sto <- do.call(rbind, out_sto)
out_sto <- as.data.table(out_sto)
out_sto <- melt(out_sto, id.vars = "run", variable.name = "Stage", value.name = "Count")
out_sto[, "Day" := 0:tmax, by = c("run", "Stage")]
levels(out_sto$Stage) <- c("L", "A", "B", "Q")

ggplot(data = out_sto) +
    geom_line(aes(x = Day, y = Count, color = Stage, group = run), alpha = 0.35) +
    geom_line(data = out_det, mapping = aes(x = Day, y = Count, color = Stage)) +
    facet_wrap(. ~ Stage, scales = "free")

