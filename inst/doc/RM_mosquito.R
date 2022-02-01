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

M_det <- matrix(data = 0, nrow = tmax, ncol = 3)
Y_det <- matrix(data = 0, nrow = tmax, ncol = 3)
Z_det <- matrix(data = 0, nrow = tmax, ncol = 3)

# run it
while(mod$global$tnow <= tmax) {
  mod$mosquito$kappa <- rep(0.05, 3)
  step_aqua(model = mod)
  step_mosquitoes(model = mod)
  M_det[mod$global$tnow, ] <- mod$mosquito$M
  Y_det[mod$global$tnow, ] <- mod$mosquito$Y
  Z_det[mod$global$tnow, ] <- mod$mosquito$Z
  mod$global$tnow <- mod$global$tnow + 1L
}

det_out <- as.data.table(rbind(M_det, Y_det, Z_det))

det_out[, "Day" := as.integer(rep(1:tmax, 3))]
det_out[, "Compartment" := rep(c("M", "Y", "Z"), times = rep(tmax, 3))]
det_out <- melt(det_out, id.vars = c("Day", "Compartment"), variable.name = "Patch", value.name = "Count")
det_out[, "Patch" := as.integer(Patch)]

ggplot(det_out) +
    geom_line(aes(x = Day, y = Count, color = Compartment)) +
    facet_grid(Patch ~ .)

## ---- fig.width=10, fig.height=8----------------------------------------------
sto_out <- mclapply(X = 1:10, FUN = function(runid) {
  
  mod <- make_MicroMoB(tmax = tmax, p = 3)
  setup_mosquito_RM(mod, stochastic = TRUE, f = f, q = q, eip = EIP, p = psurv, psi = psi, M = M, Y = Y, Z = Z)
  setup_aqua_trace(model = mod, lambda = lambda, stochastic = TRUE)
  
  M_out <- as.data.frame(matrix(data = 0, nrow = tmax, ncol = 3))
  Y_out <- as.data.frame(matrix(data = 0, nrow = tmax, ncol = 3))
  Z_out <- as.data.frame(matrix(data = 0, nrow = tmax, ncol = 3))
  
  # run it
  while(mod$global$tnow <= tmax) {
    mod$mosquito$kappa <- rep(0.05, 3)
    step_aqua(model = mod)
    step_mosquitoes(model = mod)
    M_out[mod$global$tnow, ] <- mod$mosquito$M
    Y_out[mod$global$tnow, ] <- mod$mosquito$Y
    Z_out[mod$global$tnow, ] <- mod$mosquito$Z
    mod$global$tnow <- mod$global$tnow + 1L
  }
  
  out <- as.data.table(rbind(M_out, Y_out, Z_out))
  out[, "Day" := as.integer(rep(1:tmax, 3))]
  out[, "Compartment" := rep(c("M", "Y", "Z"), times = rep(tmax, 3))]
  out <- melt(out, id.vars = c("Day", "Compartment"), variable.name = "Patch", value.name = "Count")
  out[, "Patch" := as.integer(Patch)]
  out[, "Run" := as.integer(runid)]
  
  return(out)
})

sto_out <- do.call(rbind, sto_out)

ggplot(sto_out) +
    geom_line(aes(x = Day, y = Count, color = Compartment, group = interaction(Run, Compartment)), alpha = 0.35) +
    facet_grid(Patch ~ .)

