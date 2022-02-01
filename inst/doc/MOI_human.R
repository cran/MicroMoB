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

MOI_out <- matrix(data = NaN, nrow = nrow(MOI_init), ncol = tmax + 1)
MOI_out[, 1L] <- MOI_init

while (mod$global$tnow <= mod$global$tmax) {
  mod$human$EIR <- EIR
  step_humans(model = mod)
  MOI_out[, mod$global$tnow + 1L] <- mod$human$MOI
  mod$global$tnow <- mod$global$tnow + 1L
}

## -----------------------------------------------------------------------------
weighted.mean(x = 0:100, w = MOI_out[, tmax])

## -----------------------------------------------------------------------------
# plot output
MOI_out_dt <- as.data.table(t(MOI_out))
MOI_out_dt <- suppressWarnings(melt(MOI_out_dt))
setnames(MOI_out_dt, new = c("MOI", "Count"))
levels(MOI_out_dt$MOI) <- 0:100
MOI_out_dt[, "Day" := 0:tmax, by = MOI]

ggplot(data = MOI_out_dt) +
  geom_line(aes(x = Day, y = Count, group = MOI, color = as.numeric(MOI)))

## -----------------------------------------------------------------------------
MOI_final <- MOI_out_dt[Day == tmax, ]
MOI_final[, "prop" := Count / sum(Count)]
MOI_final[, "MOI" := 0:100]
MOI_final[, "theoretical" := dpois(x = MOI, lambda = h/r)]

ggplot(MOI_final, aes(MOI, prop)) +
    geom_bar(stat = 'identity') +
    geom_line(aes(x = MOI, y = theoretical), color = "blue")

