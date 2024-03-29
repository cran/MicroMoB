test_that("RM model setup is working", {
  tmax <- 20

  f <- 0.3
  q <- 1
  a <- f * q
  psi <- matrix(
    c(0.9, 0.025, 0.075,
      0.15, 0.6, 0.25,
      0.01, 0.04, 0.95), nrow = 3, ncol = 3, byrow = TRUE
  )
  M <- c(100, 150, 120)
  Y <- c(10, 15, 8)
  Z <- c(8, 10, 4)

  mod <- make_MicroMoB(tmax = tmax, p = 3)

  # basic errors with scalar values for p, eip
  expect_error(setup_mosquito_RM(mod, stochastic = FALSE, f = f, q = q, eip = 5, p = -2, psi = psi, M = M, Y = Y, Z = Z))
  expect_error(setup_mosquito_RM(mod, stochastic = FALSE, f = f, q = q, eip = 5, p = 2, psi = psi, M = M, Y = Y, Z = Z))
  expect_error(setup_mosquito_RM(mod, stochastic = FALSE, f = f, q = q, eip = 5, p = numeric(0), psi = psi, M = M, Y = Y, Z = Z))
  expect_error(setup_mosquito_RM(mod, stochastic = FALSE, f = f, q = q, eip = -5, p = 0.9, psi = psi, M = M, Y = Y, Z = Z))
  expect_error(setup_mosquito_RM(mod, stochastic = FALSE, f = f, q = q, eip = numeric(0), p = 0.9, psi = psi, M = M, Y = Y, Z = Z))

  # eip tests
  setup_mosquito_RM(mod, stochastic = FALSE, f = f, q = q, eip = 5, p = 0.9, psi = psi, M = M, Y = Y, Z = Z)
  expect_equal(mod$mosquito$eip, rep(5, tmax))
  expect_equal(mod$mosquito$p, matrix(0.9, 3, tmax))

  mod <- make_MicroMoB(tmax = tmax, p = 3)
  eip <- 1:365
  setup_mosquito_RM(mod, stochastic = FALSE, f = f, q = q, eip = eip, p = 0.9, psi = psi, M = M, Y = Y, Z = Z)
  expect_equal(mod$mosquito$eip, 1:tmax)
  expect_equal(mod$mosquito$p, matrix(0.9, 3, tmax))

  mod <- make_MicroMoB(tmax = tmax, p = 3)
  eip <- 1:tmax
  setup_mosquito_RM(mod, stochastic = FALSE, f = f, q = q, eip = eip, p = 0.9, psi = psi, M = M, Y = Y, Z = Z)
  expect_equal(mod$mosquito$eip, 1:tmax)
  expect_equal(mod$mosquito$p, matrix(0.9, 3, tmax))

  # p tests
  mod <- make_MicroMoB(tmax = tmax, p = 3)
  p <- seq(from = 0.01, to = 0.99, length.out = 365)
  setup_mosquito_RM(mod, stochastic = FALSE, f = f, q = q, eip = 5, p = p, psi = psi, M = M, Y = Y, Z = Z)
  expect_equal(mod$mosquito$p, t(replicate(3, p[1:tmax])))

  mod <- make_MicroMoB(tmax = tmax, p = 3)
  p <- seq(from = 0.01, to = 0.99, length.out = tmax)
  setup_mosquito_RM(mod, stochastic = FALSE, f = f, q = q, eip = 5, p = p, psi = psi, M = M, Y = Y, Z = Z)
  expect_equal(mod$mosquito$p, t(replicate(3, p[1:tmax])))

  # oviposit tests
  expected_nu_det <- mod$mosquito$nu * f * M
  expect_equal(compute_oviposit(mod), expected_nu_det)

  mod <- make_MicroMoB(tmax = tmax, p = 3)
  setup_mosquito_RM(mod, stochastic = TRUE, f = f, q = q, eip = 5, p = 0.9, psi = psi, M = M, Y = Y, Z = Z)
  expect_true(all(compute_oviposit(mod) > 0))

  # other objs
  expect_equal(length(mod$mosquito$kappa), 3)
  expect_equal(length(mod$mosquito$M), 3)
  expect_equal(length(mod$mosquito$Y), 3)
  expect_equal(length(mod$mosquito$Z), 3)
  expect_equal(dim(mod$mosquito$ZZ), c(5, 3))

  # compute results
  expect_equal(compute_Z(mod), Z)
  expect_equal(compute_f(mod, rep(1 ,3)), rep(f, 3))
})


test_that("deterministic RM step is working with pulse of infection, no dispersal", {
  tmax <- 20

  f <- 0.3
  q <- 1
  a <- f * q
  psi <- diag(3)
  M <- c(100, 150, 120)
  Y <- c(0, 0, 0)
  Z <- c(0, 0, 0)

  mod <- make_MicroMoB(tmax = tmax, p = 3)
  setup_mosquito_RM(mod, stochastic = FALSE, f = f, q = q, eip = 3, p = 0.9, psi = psi, M = M, Y = Y, Z = Z)
  setup_aqua_trace(model = mod, lambda = c(1e1, 1e2, 1e3), stochastic = FALSE)

  expect_equal(mod$mosquito$Y, Y)
  expect_equal(mod$mosquito$Z, Z)

  # time = 1
  mod$mosquito$kappa <- rep(1, 3)
  step_mosquitoes(model = mod)

  expect_equal(mod$mosquito$M, (M*0.9) + c(1e1,1e2,1e3))
  expect_true(all(mod$mosquito$Y > 0))
  expect_true(all(mod$mosquito$Z == 0))
  expect_true(all(mod$mosquito$ZZ[1:2, ] == 0))
  expect_equal(mod$mosquito$ZZ[3, ], mod$mosquito$Y)

  # time = 2
  mod$mosquito$kappa <- rep(0, 3)
  mod$global$tnow <- 2
  step_mosquitoes(model = mod)

  expect_true(all(mod$mosquito$Y > 0))
  expect_true(all(mod$mosquito$Z == 0))
  expect_equal(mod$mosquito$ZZ[2, ], mod$mosquito$Y)
  expect_true(all(mod$mosquito$ZZ[-2, ] == 0))

  # time = 3
  mod$global$tnow <- 3
  step_mosquitoes(model = mod)

  expect_true(all(mod$mosquito$Y > 0))
  expect_true(all(mod$mosquito$Z == 0))
  expect_equal(mod$mosquito$ZZ[1, ], mod$mosquito$Y)
  expect_true(all(mod$mosquito$ZZ[-1, ] == 0))

  # time = 4 (expect Z mosquitoes)
  mod$global$tnow <- 4
  step_mosquitoes(model = mod)

  expect_true(all(mod$mosquito$Y > 0))
  expect_equal(mod$mosquito$Y, mod$mosquito$Z)
  expect_true(all(mod$mosquito$ZZ == 0))

  # by hand
  expect_equal((M * a) * (0.9^4), mod$mosquito$Z)

  out <- output_mosquitoes(mod)
  expect_equal(out$M, mod$mosquito$M)
  expect_equal(out$Y, mod$mosquito$Y)
  expect_equal(out$Z, mod$mosquito$Z)

})


test_that("deterministic RM step is working with pulse of infection, no dispersal, 2 EIP vals", {
  tmax <- 20

  f <- 0.3
  q <- 1
  a <- f * q
  psi <- diag(3)
  M <- c(100, 150, 120)
  Y <- c(0, 0, 0)
  Z <- c(0, 0, 0)

  mod <- make_MicroMoB(tmax = tmax, p = 3)
  setup_mosquito_RM(mod, stochastic = FALSE, f = f, q = q, eip = rep(c(4, 1), times = c(1, 19)), p = 0.9, psi = psi, M = M, Y = Y, Z = Z)
  setup_aqua_trace(model = mod, lambda = c(1e1, 1e2, 1e3), stochastic = FALSE)

  expect_equal(mod$mosquito$Y, Y)
  expect_equal(mod$mosquito$Z, Z)

  # time = 1
  mod$mosquito$kappa <- rep(1, 3)
  step_mosquitoes(model = mod)

  expect_true(all(mod$mosquito$Y > 0))
  expect_true(all(mod$mosquito$Z == 0))
  expect_true(all(mod$mosquito$ZZ[-4, ] == 0))
  expect_equal(mod$mosquito$ZZ[4, ], mod$mosquito$Y)

  # time = 2
  mod$mosquito$kappa <- rep(1, 3)
  mod$global$tnow <- 2
  step_mosquitoes(model = mod)

  expect_true(all(mod$mosquito$Y > 0))
  expect_true(all(mod$mosquito$Z == 0))
  expect_equal(colSums(mod$mosquito$ZZ), mod$mosquito$Y)
  expect_true(all(mod$mosquito$ZZ[c(2, 4), ] == 0))

  # time = 3 (Z mosquitoes from those infected on t=2)
  mod$mosquito$kappa <- rep(0, 3)
  mod$global$tnow <- 3
  step_mosquitoes(model = mod)

  expect_true(all(mod$mosquito$Y > 0))
  expect_true(all(mod$mosquito$Z > 0))
  expect_equal(mod$mosquito$Y - colSums(mod$mosquito$ZZ), mod$mosquito$Z)
  expect_true(all(mod$mosquito$ZZ[-2, ] == 0))

  # time = 4
  mod$global$tnow <- 4
  step_mosquitoes(model = mod)

  expect_true(all(mod$mosquito$Y > 0))
  expect_true(all(mod$mosquito$Z > 0))
  expect_equal(mod$mosquito$Y - colSums(mod$mosquito$ZZ), mod$mosquito$Z)
  expect_true(all(mod$mosquito$ZZ[-1, ] == 0))

  # time = 5 (Z mosquitoes from those infected on t=1)
  mod$global$tnow <- 5
  step_mosquitoes(model = mod)

  expect_true(all(mod$mosquito$Y > 0))
  expect_true(all(mod$mosquito$Z > 0))
  expect_equal(mod$mosquito$Y, mod$mosquito$Z)
  expect_true(all(mod$mosquito$ZZ == 0))

})

test_that("deterministic RM step is working with pulse of infection, with dispersal", {
  tmax <- 20

  f <- 0.3
  q <- 1
  a <- f * q
  psi <- matrix(
    c(0.9, 0.025, 0.075,
      0.1, 0.7, 0.2,
      0.05, 0.15, 0.8),
    nrow = 3, ncol = 3, byrow = TRUE
  )
  M <- c(100, 150, 120)
  Y <- c(0, 0, 0)
  Z <- c(0, 0, 0)

  mod <- make_MicroMoB(tmax = tmax, p = 3)
  setup_mosquito_RM(mod, stochastic = FALSE, f = f, q = q, eip = 2, p = 0.9, psi = psi, M = M, Y = Y, Z = Z)
  setup_aqua_trace(model = mod, lambda = c(1e1, 1e2, 1e3), stochastic = FALSE)

  expect_equal(mod$mosquito$Y, Y)
  expect_equal(mod$mosquito$Z, Z)

  # time = 1
  mod$mosquito$kappa <- rep(1, 3)
  step_mosquitoes(model = mod)

  expect_equal(mod$mosquito$M, as.vector(((M*0.9)) %*% psi) + c(1e1,1e2,1e3))
  expect_true(all(mod$mosquito$Y > 0))
  expect_true(all(mod$mosquito$Z == 0))
  expect_true(all(mod$mosquito$ZZ[1, ] == 0))
  expect_equal(mod$mosquito$ZZ[2, ], mod$mosquito$Y)

  # time = 2
  mod$mosquito$kappa <- rep(0, 3)
  mod$global$tnow <- 2
  step_mosquitoes(model = mod)

  expect_true(all(mod$mosquito$Y > 0))
  expect_true(all(mod$mosquito$Z == 0))
  expect_equal(mod$mosquito$ZZ[1, ], mod$mosquito$Y)
  expect_true(all(mod$mosquito$ZZ[2, ] == 0))

  # time = 3 (expect Z mosquitoes)
  mod$global$tnow <- 3
  step_mosquitoes(model = mod)

  expect_true(all(mod$mosquito$Y > 0))
  expect_true(all(mod$mosquito$Z > 0))
  expect_true(all(mod$mosquito$ZZ == 0))

  # by hand
  expected_Z <- ((M * a) * (0.9^3)) %*% psi %*% psi %*% psi
  expect_equal(as.vector(expected_Z), mod$mosquito$Z)

})



test_that("stochastic RM step is working with pulse of infection, no dispersal", {
  tmax <- 20

  f <- 0.3
  q <- 1
  a <- f * q
  psi <- diag(3)
  M <- c(1e5, 5e5, 2e5)
  Y <- c(0, 0, 0)
  Z <- c(0, 0, 0)

  mod <- make_MicroMoB(tmax = tmax, p = 3)
  setup_mosquito_RM(mod, stochastic = TRUE, f = f, q = q, eip = 3, p = 0.9, psi = psi, M = M, Y = Y, Z = Z)
  setup_aqua_trace(model = mod, lambda = c(1e1, 1e2, 1e3), stochastic = FALSE)

  expect_equal(mod$mosquito$Y, Y)
  expect_equal(mod$mosquito$Z, Z)

  # time = 1
  mod$mosquito$kappa <- rep(1, 3)
  step_mosquitoes(model = mod)

  expect_true(all(mod$mosquito$Y > 0))
  expect_true(all(mod$mosquito$Z == 0))
  expect_true(all(mod$mosquito$ZZ[1:2, ] == 0))
  expect_equal(mod$mosquito$ZZ[3, ], mod$mosquito$Y)

  # time = 2
  mod$mosquito$kappa <- rep(0, 3)
  mod$global$tnow <- 2
  step_mosquitoes(model = mod)

  expect_true(all(mod$mosquito$Y > 0))
  expect_true(all(mod$mosquito$Z == 0))
  expect_equal(mod$mosquito$ZZ[2, ], mod$mosquito$Y)
  expect_true(all(mod$mosquito$ZZ[-2, ] == 0))

  # time = 3
  mod$global$tnow <- 3
  step_mosquitoes(model = mod)

  expect_true(all(mod$mosquito$Y > 0))
  expect_true(all(mod$mosquito$Z == 0))
  expect_equal(mod$mosquito$ZZ[1, ], mod$mosquito$Y)
  expect_true(all(mod$mosquito$ZZ[-1, ] == 0))

  # time = 4 (expect Z mosquitoes)
  mod$global$tnow <- 4
  step_mosquitoes(model = mod)

  expect_true(all(mod$mosquito$Y > 0))
  expect_equal(mod$mosquito$Y, mod$mosquito$Z)
  expect_true(all(mod$mosquito$ZZ == 0))

  expected <-  (M * a) * (0.9^4)
  expect_true(all(abs(mod$mosquito$Z - expected) / expected < 0.025))

  out <- output_mosquitoes(mod)
  expect_equal(out$M, mod$mosquito$M)
  expect_equal(out$Y, mod$mosquito$Y)
  expect_equal(out$Z, mod$mosquito$Z)

})


test_that("stochastic RM step is working with pulse of infection, no dispersal, 2 EIP vals", {
  tmax <- 20

  f <- 0.3
  q <- 1
  a <- f * q
  psi <- diag(3)
  M <- c(100, 150, 120)
  Y <- c(0, 0, 0)
  Z <- c(0, 0, 0)

  mod <- make_MicroMoB(tmax = tmax, p = 3)
  setup_mosquito_RM(mod, stochastic = TRUE, f = f, q = q, eip = rep(c(4, 1), times = c(1, 19)), p = 0.9, psi = psi, M = M, Y = Y, Z = Z)
  setup_aqua_trace(model = mod, lambda = c(1e1, 1e2, 1e3), stochastic = FALSE)

  expect_equal(mod$mosquito$Y, Y)
  expect_equal(mod$mosquito$Z, Z)

  # time = 1
  mod$mosquito$kappa <- rep(1, 3)
  step_mosquitoes(model = mod)

  expect_true(all(mod$mosquito$Y > 0))
  expect_true(all(mod$mosquito$Z == 0))
  expect_true(all(mod$mosquito$ZZ[-4, ] == 0))
  expect_equal(mod$mosquito$ZZ[4, ], mod$mosquito$Y)

  # time = 2
  mod$mosquito$kappa <- rep(1, 3)
  mod$global$tnow <- 2
  step_mosquitoes(model = mod)

  expect_true(all(mod$mosquito$Y > 0))
  expect_true(all(mod$mosquito$Z == 0))
  expect_equal(colSums(mod$mosquito$ZZ), mod$mosquito$Y)
  expect_true(all(mod$mosquito$ZZ[c(2, 4), ] == 0))

  # time = 3 (Z mosquitoes from those infected on t=2)
  mod$mosquito$kappa <- rep(0, 3)
  mod$global$tnow <- 3
  step_mosquitoes(model = mod)

  expect_true(all(mod$mosquito$Y > 0))
  expect_true(all(mod$mosquito$Z > 0))
  expect_equal(mod$mosquito$Y - colSums(mod$mosquito$ZZ), mod$mosquito$Z)
  expect_true(all(mod$mosquito$ZZ[-2, ] == 0))

  # time = 4
  mod$global$tnow <- 4
  step_mosquitoes(model = mod)

  expect_true(all(mod$mosquito$Y > 0))
  expect_true(all(mod$mosquito$Z > 0))
  expect_equal(mod$mosquito$Y - colSums(mod$mosquito$ZZ), mod$mosquito$Z)
  expect_true(all(mod$mosquito$ZZ[-1, ] == 0))

  # time = 5 (Z mosquitoes from those infected on t=1)
  mod$global$tnow <- 5
  step_mosquitoes(model = mod)

  expect_true(all(mod$mosquito$Y > 0))
  expect_true(all(mod$mosquito$Z > 0))
  expect_equal(mod$mosquito$Y, mod$mosquito$Z)
  expect_true(all(mod$mosquito$ZZ == 0))

})


test_that("stochastic RM step is working with pulse of infection, with dispersal", {
  tmax <- 20

  f <- 0.3
  q <- 1
  a <- f * q
  psi <- matrix(
    c(0.9, 0.025, 0.075,
      0.1, 0.7, 0.2,
      0.05, 0.15, 0.8),
    nrow = 3, ncol = 3, byrow = TRUE
  )
  M <- c(1e5, 5e5, 2e5)
  Y <- c(0, 0, 0)
  Z <- c(0, 0, 0)

  mod <- make_MicroMoB(tmax = tmax, p = 3)
  setup_mosquito_RM(mod, stochastic = TRUE, f = f, q = q, eip = 2, p = 0.9, psi = psi, M = M, Y = Y, Z = Z)
  setup_aqua_trace(model = mod, lambda = c(1e1, 1e2, 1e3), stochastic = FALSE)

  expect_equal(mod$mosquito$Y, Y)
  expect_equal(mod$mosquito$Z, Z)

  # time = 1
  mod$mosquito$kappa <- rep(1, 3)
  step_mosquitoes(model = mod)

  # expect_equal(mod$mosquito$M, as.vector(((M*0.9) + c(1e1,1e2,1e3)) %*% psi))
  expect_true(all(mod$mosquito$Y > 0))
  expect_true(all(mod$mosquito$Z == 0))
  expect_true(all(mod$mosquito$ZZ[1, ] == 0))
  expect_equal(colSums(mod$mosquito$ZZ), mod$mosquito$Y)

  # time = 2
  mod$mosquito$kappa <- rep(0, 3)
  mod$global$tnow <- 2
  step_mosquitoes(model = mod)

  expect_true(all(mod$mosquito$Y > 0))
  expect_true(all(mod$mosquito$Z == 0))
  expect_equal(colSums(mod$mosquito$ZZ), mod$mosquito$Y)
  expect_true(all(mod$mosquito$ZZ[2, ] == 0))

  # time = 3 (expect Z mosquitoes)
  mod$global$tnow <- 3
  step_mosquitoes(model = mod)

  expect_true(all(mod$mosquito$Y > 0))
  expect_true(all(mod$mosquito$Z > 0))
  expect_true(all(mod$mosquito$ZZ == 0))

  # by hand
  expected_Z <- ((M * a) * (0.9^3)) %*% psi %*% psi %*% psi
  expected_Z <- as.vector(expected_Z)
  expect_true(all(abs(mod$mosquito$Z - expected_Z) / expected_Z < 0.025))

})


test_that("test JSON config working", {

  library(jsonlite)

  t <- 10 # days to simulate
  p <- 2 # number of patches

  EIP <-  rep(5, t)
  p_surv <- 0.95
  psi <- matrix(rexp(p^2), nrow = p, ncol = p)
  psi <- psi / rowSums(psi)

  # sending to JSON does not change R type when read back in
  par <- list(
    "stochastic" = FALSE,
    "f" = 0.3,
    "q" = 0.9,
    "eip" = EIP,
    "p" = p_surv,
    "psi" = psi,
    "nu" = 20,
    "M" = rep(100, p),
    "Y" = rep(20, p),
    "Z" = rep(5, p)
  )

  json_path <- tempfile(pattern = "mosquito_par", fileext = ".json")
  write_json(x = par, path = json_path, digits = NA)
  par_in <- get_config_mosquito_RM(path = json_path)
  expect_true(all.equal(par, par_in))

  # reject obviously bad input
  par <- list(
    "stochastic" = FALSE,
    "f" = 0.3,
    "q" = 0.9,
    "eip" = EIP,
    "p" = p_surv,
    "psi" = NULL,
    "nu" = 20,
    "M" = rep(100, p),
    "Y" = rep(20, p),
    "Z" = rep(5, p)
  )

  json_path <- tempfile(pattern = "mosquito_par", fileext = ".json")
  write_json(x = par, path = json_path, digits = NA)
  expect_error(get_config_mosquito_RM(path = json_path))

  unlink(x = json_path)

})


test_that("JSON parameters can read in", {
  path <- system.file("extdata", "mosquito_RM.json", package = "MicroMoB")
  pars <- get_config_mosquito_RM(path = path)

  expect_true(length(pars) == 10L)
  expect_true(is.logical(pars$stochastic))
  expect_true(length(pars$stochastic) == 1L)

  expect_true(is.numeric(pars$f))
  expect_true(is.numeric(pars$q))

  expect_true(is.numeric(pars$eip))
  expect_true(is.vector(pars$eip))

  expect_true(is.numeric(pars$p))
  expect_true(is.vector(pars$p))

  expect_true(is.numeric(pars$psi))
  expect_true(is.matrix(pars$psi))

  expect_true(is.numeric(pars$nu))

  expect_true(is.numeric(pars$M))
  expect_true(is.numeric(pars$Y))
  expect_true(is.numeric(pars$Z))
  expect_true(is.vector(pars$M))
  expect_true(is.vector(pars$Y))
  expect_true(is.vector(pars$Z))

})
