test_that("svlk() returns the correct covariance matrices and derivatives", {

  # Testing the GxExM case =====================================================
  LK <- matrix(c(1.0, 0.8, 0.7,
                 0.8, 1.0, 0.9,
                 0.7, 0.9, 1.0),
               nrow = 3, byrow = TRUE)

  order <- 9
  kappa <- c(2, 3, 4, 0.2, 0.4, 0.5)
  vf <- svlk(C = LK)
  result <- vf(order = order, kappa = kappa)

  ## Covariance matrix =========================================================
  expectation <- rbind(cbind(kappa[1] * LK,    sqrt(kappa[1]) * sqrt(kappa[2]) * kappa[4] * LK,    sqrt(kappa[1]) * sqrt(kappa[3]) * kappa[5] * LK),
                       cbind(sqrt(kappa[1]) * sqrt(kappa[2]) * kappa[4] * LK,    kappa[2] * LK,    sqrt(kappa[2]) * sqrt(kappa[3]) * kappa[6] * LK),
                       cbind(sqrt(kappa[1]) * sqrt(kappa[3]) * kappa[5] * LK,    sqrt(kappa[2]) * sqrt(kappa[3]) * kappa[6] * LK,    kappa[3] * LK))

  expect_equal(result[[1]], expectation)

  ## Derivative wrt variance 1 =================================================
  expectation <- rbind(cbind(LK,    0.5 * sqrt(kappa[2]) * kappa[4] * LK / sqrt(kappa[1]),    0.5 * sqrt(kappa[3]) * kappa[5] * LK / sqrt(kappa[1])),
                       cbind(0.5 * sqrt(kappa[2]) * kappa[4] * LK / sqrt(kappa[1]),    0 * LK,    0 * LK),
                       cbind(0.5 * sqrt(kappa[3]) * kappa[5] * LK / sqrt(kappa[1]),    0 * LK,    0 * LK))

  expect_equal(result[[2]], expectation)

  ## Derivative wrt variance 2 =================================================
  expectation <- rbind(cbind(0 * LK,    0.5 * sqrt(kappa[1]) * kappa[4] * LK / sqrt(kappa[2]),    0 * LK),
                       cbind(0.5 * sqrt(kappa[1]) * kappa[4] * LK / sqrt(kappa[2]),    LK,    0.5 * sqrt(kappa[3]) * kappa[6] * LK / sqrt(kappa[2])),
                       cbind(0 * LK,    0.5 * sqrt(kappa[3]) * kappa[6] * LK / sqrt(kappa[2]),    0 * LK))

  expect_equal(result[[3]], expectation)

  ## Derivative wrt variance 3 =================================================
  expectation <- rbind(cbind(0 * LK,    0 * LK,    0.5 * sqrt(kappa[1]) * kappa[5] * LK / sqrt(kappa[3])),
                       cbind(0 * LK,    0 * LK,    0.5 * sqrt(kappa[2]) * kappa[6] * LK / sqrt(kappa[3])),
                       cbind(0.5 * sqrt(kappa[1]) * kappa[5] * LK / sqrt(kappa[3]),    0.5 * sqrt(kappa[2]) * kappa[6] * LK / sqrt(kappa[3]),    LK))

  expect_equal(result[[4]], expectation)

  ## Derivative wrt correlation 1-2 ============================================
  expectation <- rbind(cbind(0 * LK,    sqrt(kappa[1]) * sqrt(kappa[2]) * LK,    0 * LK),
                       cbind(sqrt(kappa[1]) * sqrt(kappa[2]) * LK,    0 * LK,    0 * LK),
                       cbind(0 * LK,    0 * LK,    0 * LK))

  expect_equal(result[[5]], expectation)

  ## Derivative wrt correlation 1-3 ============================================
  expectation <- rbind(cbind(0 * LK,    0 * LK,    sqrt(kappa[1]) * sqrt(kappa[3]) * LK),
                       cbind(0 * LK,    0 * LK,    0 * LK),
                       cbind(sqrt(kappa[1]) * sqrt(kappa[3]) * LK,    0 * LK,    0 * LK))

  expect_equal(result[[6]], expectation)

  ## Derivative wrt correlation 2-3 ============================================
  expectation <- rbind(cbind(0 * LK,    0 * LK,    0 * LK),
                       cbind(0 * LK,    0 * LK,    sqrt(kappa[2]) * sqrt(kappa[3]) * LK),
                       cbind(0 * LK,    sqrt(kappa[2]) * sqrt(kappa[3]) * LK,    0 * LK))

  expect_equal(result[[7]], expectation)

  # Testing the GxE case =======================================================
  LK <- matrix(c(1.0, 0.8, 0.7,
                 0.8, 1.0, 0.9,
                 0.7, 0.9, 1.0),
               nrow = 3, byrow = TRUE)

  order <- 3
  kappa <- c(2)
  vf <- svlk(C = LK)
  result <- vf(order = order, kappa = kappa)

  ## Covariance matrix =========================================================
  expectation <- kappa[1] * LK
  expect_equal(result[[1]], expectation)

  ## Derivative wrt variance 1 =================================================
  expectation <- LK
  expect_equal(result[[2]], expectation)
})


