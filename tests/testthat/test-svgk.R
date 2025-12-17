test_that("svgk() returns the correct covariance matrices and derivatives", {

  # Testing the GxExM case =====================================================
  LK <- matrix(c(1.0, 0.8, 0.7,
                 0.8, 1.0, 0.9,
                 0.7, 0.9, 1.0),
               nrow = 3, byrow = TRUE)
  DM <- log(LK) / -0.15

  order <- 9
  kappa <- c(2, 3, 4, 0.2, 0.4, 0.5, 0.1)
  vf <- svgk(D = DM)
  result <- vf(order = order, kappa = kappa)

  ## Covariance matrix =========================================================
  expectation <- rbind(cbind(kappa[1] * exp(-kappa[7] * DM),    sqrt(kappa[1]) * sqrt(kappa[2]) * kappa[4] * exp(-kappa[7] * DM),    sqrt(kappa[1]) * sqrt(kappa[3]) * kappa[5] * exp(-kappa[7] * DM)),
                       cbind(sqrt(kappa[1]) * sqrt(kappa[2]) * kappa[4] * exp(-kappa[7] * DM),    kappa[2] * exp(-kappa[7] * DM),    sqrt(kappa[2]) * sqrt(kappa[3]) * kappa[6] * exp(-kappa[7] * DM)),
                       cbind(sqrt(kappa[1]) * sqrt(kappa[3]) * kappa[5] * exp(-kappa[7] * DM),    sqrt(kappa[2]) * sqrt(kappa[3]) * kappa[6] * exp(-kappa[7] * DM),    kappa[3] * exp(-kappa[7] * DM)))

  expect_equal(result[[1]], expectation)

  ## Derivative wrt variance 1 =================================================
  expectation <- rbind(cbind(exp(-kappa[7] * DM),    0.5 * sqrt(kappa[2]) * kappa[4] * exp(-kappa[7] * DM) / sqrt(kappa[1]),    0.5 * sqrt(kappa[3]) * kappa[5] * exp(-kappa[7] * DM) / sqrt(kappa[1])),
                       cbind(0.5 * sqrt(kappa[2]) * kappa[4] * exp(-kappa[7] * DM) / sqrt(kappa[1]),    0 * exp(-kappa[7] * DM),    0 * exp(-kappa[7] * DM)),
                       cbind(0.5 * sqrt(kappa[3]) * kappa[5] * exp(-kappa[7] * DM) / sqrt(kappa[1]),    0 * exp(-kappa[7] * DM),    0 * exp(-kappa[7] * DM)))

  expect_equal(result[[2]], expectation)

  ## Derivative wrt variance 2 =================================================
  expectation <- rbind(cbind(0 * exp(-kappa[7] * DM),    0.5 * sqrt(kappa[1]) * kappa[4] * exp(-kappa[7] * DM) / sqrt(kappa[2]),    0 * exp(-kappa[7] * DM)),
                       cbind(0.5 * sqrt(kappa[1]) * kappa[4] * exp(-kappa[7] * DM) / sqrt(kappa[2]),    exp(-kappa[7] * DM),    0.5 * sqrt(kappa[3]) * kappa[6] * exp(-kappa[7] * DM) / sqrt(kappa[2])),
                       cbind(0 * exp(-kappa[7] * DM),    0.5 * sqrt(kappa[3]) * kappa[6] * exp(-kappa[7] * DM) / sqrt(kappa[2]),    0 * exp(-kappa[7] * DM)))

  expect_equal(result[[3]], expectation)

  ## Derivative wrt variance 3 =================================================
  expectation <- rbind(cbind(0 * exp(-kappa[7] * DM),    0 * exp(-kappa[7] * DM),    0.5 * sqrt(kappa[1]) * kappa[5] * exp(-kappa[7] * DM) / sqrt(kappa[3])),
                       cbind(0 * exp(-kappa[7] * DM),    0 * exp(-kappa[7] * DM),    0.5 * sqrt(kappa[2]) * kappa[6] * exp(-kappa[7] * DM) / sqrt(kappa[3])),
                       cbind(0.5 * sqrt(kappa[1]) * kappa[5] * exp(-kappa[7] * DM) / sqrt(kappa[3]),    0.5 * sqrt(kappa[2]) * kappa[6] * exp(-kappa[7] * DM) / sqrt(kappa[3]),    exp(-kappa[7] * DM)))

  expect_equal(result[[4]], expectation)

  ## Derivative wrt correlation 1-2 ============================================
  expectation <- rbind(cbind(0 * exp(-kappa[7] * DM),    sqrt(kappa[1]) * sqrt(kappa[2]) * exp(-kappa[7] * DM),    0 * exp(-kappa[7] * DM)),
                       cbind(sqrt(kappa[1]) * sqrt(kappa[2]) * exp(-kappa[7] * DM),    0 * exp(-kappa[7] * DM),    0 * exp(-kappa[7] * DM)),
                       cbind(0 * exp(-kappa[7] * DM),    0 * exp(-kappa[7] * DM),    0 * exp(-kappa[7] * DM)))

  expect_equal(result[[5]], expectation)

  ## Derivative wrt correlation 1-3 ============================================
  expectation <- rbind(cbind(0 * exp(-kappa[7] * DM),    0 * exp(-kappa[7] * DM),    sqrt(kappa[1]) * sqrt(kappa[3]) * exp(-kappa[7] * DM)),
                       cbind(0 * exp(-kappa[7] * DM),    0 * exp(-kappa[7] * DM),    0 * exp(-kappa[7] * DM)),
                       cbind(sqrt(kappa[1]) * sqrt(kappa[3]) * exp(-kappa[7] * DM),    0 * exp(-kappa[7] * DM),    0 * exp(-kappa[7] * DM)))

  expect_equal(result[[6]], expectation)

  ## Derivative wrt correlation 2-3 ============================================
  expectation <- rbind(cbind(0 * exp(-kappa[7] * DM),    0 * exp(-kappa[7] * DM),    0 * exp(-kappa[7] * DM)),
                       cbind(0 * exp(-kappa[7] * DM),    0 * exp(-kappa[7] * DM),    sqrt(kappa[2]) * sqrt(kappa[3]) * exp(-kappa[7] * DM)),
                       cbind(0 * exp(-kappa[7] * DM),    sqrt(kappa[2]) * sqrt(kappa[3]) * exp(-kappa[7] * DM),    0 * exp(-kappa[7] * DM)))

  expect_equal(result[[7]], expectation)

  ## Derivative wrt the bandwidth ==============================================
  expectation <- rbind(cbind(-DM * kappa[1] * exp(-kappa[7] * DM),    -DM * sqrt(kappa[1]) * sqrt(kappa[2]) * kappa[4] * exp(-kappa[7] * DM),    -DM * sqrt(kappa[1]) * sqrt(kappa[3]) * kappa[5] * exp(-kappa[7] * DM)),
                       cbind(-DM * sqrt(kappa[1]) * sqrt(kappa[2]) * kappa[4] * exp(-kappa[7] * DM),    -DM * kappa[2] * exp(-kappa[7] * DM),    -DM * sqrt(kappa[2]) * sqrt(kappa[3]) * kappa[6] * exp(-kappa[7] * DM)),
                       cbind(-DM * sqrt(kappa[1]) * sqrt(kappa[3]) * kappa[5] * exp(-kappa[7] * DM),    -DM * sqrt(kappa[2]) * sqrt(kappa[3]) * kappa[6] * exp(-kappa[7] * DM),    -DM * kappa[3] * exp(-kappa[7] * DM)))

  expect_equal(result[[8]], expectation)

  # Testing the GxE case =======================================================
  LK <- matrix(c(1.0, 0.8, 0.7,
                 0.8, 1.0, 0.9,
                 0.7, 0.9, 1.0),
               nrow = 3, byrow = TRUE)
  DM <- log(LK) / -0.15

  order <- 3
  kappa <- c(2, 0.1)
  vf <- svgk(D = DM)
  result <- vf(order = order, kappa = kappa)

  ## Covariance matrix =========================================================
  expectation <- kappa[1] * exp(-kappa[2] * DM)
  expect_equal(result[[1]], expectation)

  ## Derivative wrt variance 1 =================================================
  expectation <- exp(-kappa[2] * DM)
  expect_equal(result[[2]], expectation)

  ## Derivative wrt the bandwidth ==============================================
  expectation <- -DM * kappa[1] * exp(-kappa[2] * DM)
  expect_equal(result[[3]], expectation)
})


