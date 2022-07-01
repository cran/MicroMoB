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
h <- 0.025
r <- 1/200
b <- 0.55

EIR <- -log(1 - h) / b

n <- 1
tmax <- 1e3
MOI_init <- matrix(data = c(1e5, rep(0, 1e2)), nrow = 101, ncol = n)

mod <- make_MicroMoB(tmax = tmax, p = 1)
setup_humans_MOI(model = mod, stochastic = TRUE, theta = matrix(1, nrow = n, ncol = 1), H = colSums(MOI_init), MOI = MOI_init, r = r, b = b)

human_out <- data.table::CJ(day = 1:tmax, MOI = 0:(nrow(MOI_init)-1), value = NaN)
data.table::setkey(human_out, day)

while (get_tnow(mod) <= mod$global$tmax) {
  mod$human$EIR <- EIR
  step_humans(model = mod)
  human_out[day==get_tnow(mod), value := as.vector(mod$human$MOI)]
  mod$global$tnow <- mod$global$tnow + 1L
}

## -----------------------------------------------------------------------------
weighted.mean(x = 0:100, w = human_out[day == tmax, value])

## -----------------------------------------------------------------------------
ggplot(data = human_out) +
  geom_line(aes(x = day, y = value, group = MOI, color = MOI))

## -----------------------------------------------------------------------------
human_final <- human_out[day == tmax, ]
human_final[, 'theoretical' := dpois(x = MOI, lambda = h/r)]
human_final[, 'empirical' := value / sum(value)]

ggplot(human_final, aes(MOI, empirical)) +
    geom_bar(stat = 'identity') +
    geom_line(aes(x = MOI, y = theoretical), color = "blue")

