library(testthat)
library(AMRSpreadModel)

test_that("create_model works", {
  results <- create_model(S0 = 0.9, Ds0 = 0.05, As0 = 0.025, Ts0 = 0.0125,
                          DR0 = 0.005, AR0 = 0.0025, TR0 = 0.001,
                          a = 0.3, b_h = 0.1, m = 0.7, Phi = 0.7,
                          fT = 0.2, rD = 0.05, rA = 0.03, rTs = 0.02,
                          rTR = 0.01, Sv0 = 0.9, Ev_s0 = 0, Iv_s0 = 0,
                          Ev_r0 = 0, Iv_r0 = 0,
                          e = 0.1, mu = 0.1, n = 10,
                          Lambda_v_s = 0.005, Lambda_v_r = 0.003,
                          EIR_s = NULL, EIR_R = NULL)
  model_results <- results$model$run(0:500)
  expect_true(!is.null(model_results))
  expect_true(length(model_results) > 0)

  expect_equal(results$FOI$EIR_s, 0.7 * 0.3 * 0)
  expect_equal(results$FOI$EIR_R, 0.7 * 0.3 * 0)

  EIR_s <- calculate_eir(m = 0.7, a = 0.3, Iv = 0)
  EIR_R <- calculate_eir(m = 0.7, a = 0.3, Iv = 0)
  expect_equal(results$EIR_s, EIR_s)
  expect_equal(results$EIR_R, EIR_R)

  expect_equal(results$Prevalence, 0.05 + 0.025 + 0.0125 + 0.005 + 0.0025 + 0.001)
  expect_equal(results$Prevalence_R, (0.005 + 0.0025 + 0.001) / (0.005 + 0.0025 + 0.001 + 0.025 + 0.05 + 0.0125))
})

test_that("run_model_for_ranges works", {
  results <- run_model_for_ranges(fT_range = seq(0.1, 0.9, by = 0.1),
                                  EIR_range = seq(10, 200, by = 10),
                                  S0 = 0.9, Ds0 = 0.05, As0 = 0.025, Ts0 = 0.0125,
                                  DR0 = 0.005, AR0 = 0.0025, TR0 = 0.001,
                                  b_h = 0.1, Phi = 0.7,
                                  rD = 0.05, rA = 0.03, rTs = 0.02, rTR = 0.01,
                                  Sv0 = 0.9, Ev_s0 = 0, Iv_s0 = 0, Ev_r0 = 0, Iv_r0 = 0,
                                  e = 0.1, mu = 0.1, n = 10, Lambda_v_s = 0.005, Lambda_v_r = 0.003)
  expect_true(!is.null(results))
  expect_true(nrow(results) > 0)
})

test_that("calculate_eir works", {
  EIR <- calculate_eir(m = 0.7, a = 0.3, Iv = 0.05)
  expect_equal(EIR, 0.7 * 0.3 * 0.05 * 365)
})
