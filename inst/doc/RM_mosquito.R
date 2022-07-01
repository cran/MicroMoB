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

## ---- echo=FALSE--------------------------------------------------------------
pardefault <- par()$mar

## ---- fig.width=8, fig.height=6-----------------------------------------------
tmax <- 365 * 3
p <- 3

lambda <- dnorm(x = 1:365, mean = 180, sd = 90)
lambda <- lambda * (100/max(lambda))
lambda <- t(replicate(p, lambda))

psurv <- (sin((1:365)/365*2*pi) + 1.01)/2 * 0.9

EIP <- 5

f <- 0.3
q <- 1
psi <- matrix(
  c(
    0.9, 0.05, 0.05,
    0.05, 0.9, 0.05,
    0.05, 0.05, 0.9
  ), nrow = 3, ncol = 3, 
  byrow = TRUE
)


par(mar = c(5,5,2,5))
plot(lambda[1, ], type = "l", col = "red", xlab = "Day", ylab = "Lambda (red)")
par(new = TRUE)
plot(psurv, type = "l", axes = F, xlab = NA, ylab = NA, col = "blue")
axis(side = 4)
mtext(side = 4, line = 3, 'Survival Probability (blue)')

## ---- echo=FALSE--------------------------------------------------------------
par(mar = pardefault)

## ---- fig.width=10, fig.height=8----------------------------------------------
M <- c(100, 100, 100)
Y <- c(0, 0, 0)
Z <- c(0, 0, 0)

mod <- make_MicroMoB(tmax = tmax, p = 3)
setup_mosquito_RM(mod, stochastic = FALSE, f = f, q = q, eip = EIP, p = psurv, psi = psi, M = M, Y = Y, Z = Z)
setup_aqua_trace(model = mod, lambda = lambda, stochastic = FALSE)

det_out <- data.table::CJ(day = 1:tmax, state = c('M', 'Y', 'Z'), patch = 1:3, value = NaN)
det_out <- det_out[c('M', 'Y', 'Z'), on="state"]
data.table::setkey(det_out, day)

# run it
while(get_tnow(mod) <= tmax) {
  
  mod$mosquito$kappa <- rep(0.05, 3)
  step_aqua(model = mod)
  step_mosquitoes(model = mod)
  
  det_out[day == get_tnow(mod) & state == 'M', value := mod$mosquito$M]
  det_out[day == get_tnow(mod) & state == 'Y', value := mod$mosquito$Y]
  det_out[day == get_tnow(mod) & state == 'Z', value := mod$mosquito$Z]
  
  mod$global$tnow <- mod$global$tnow + 1L
}

ggplot(det_out) +
    geom_line(aes(x = day, y = value, color = state)) +
    facet_grid(patch ~ .)

## ---- fig.width=10, fig.height=8----------------------------------------------
sto_out <- mclapply(X = 1:10, FUN = function(runid) {
  
  mod <- make_MicroMoB(tmax = tmax, p = 3)
  setup_mosquito_RM(mod, stochastic = TRUE, f = f, q = q, eip = EIP, p = psurv, psi = psi, M = M, Y = Y, Z = Z)
  setup_aqua_trace(model = mod, lambda = lambda, stochastic = TRUE)
  
  out <- data.table::CJ(day = 1:tmax, state = c('M', 'Y', 'Z'), patch = 1:3, value = NaN)
  out <- out[c('M', 'Y', 'Z'), on="state"]
  data.table::setkey(out, day)
  
  while(get_tnow(mod) <= tmax) {
    
    mod$mosquito$kappa <- rep(0.05, 3)
    step_aqua(model = mod)
    step_mosquitoes(model = mod)
    
    out[day == get_tnow(mod) & state == 'M', value := mod$mosquito$M]
    out[day == get_tnow(mod) & state == 'Y', value := mod$mosquito$Y]
    out[day == get_tnow(mod) & state == 'Z', value := mod$mosquito$Z]
    
    mod$global$tnow <- mod$global$tnow + 1L
  }
  
  out[, 'run' := as.integer(runid)]
  return(out)
})

sto_out <- data.table::rbindlist(sto_out)

ggplot(sto_out) +
    geom_line(aes(x = day, y = value, color = state, group = interaction(run, state)), alpha = 0.35) +
    facet_grid(patch ~ .)

