test_that("mvgk() returns the correct covariance matrices and derivatives", {

  # Testing the GxExM case =====================================================
  LK <- matrix(c(1.0, 0.8, 0.7,
                 0.8, 1.0, 0.9,
                 0.7, 0.9, 1.0),
               nrow = 3, byrow = TRUE)
  DM <- log(LK) / -0.15

  order <- 9
  kappa <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 0.2, 0.4, 0.5, 0.1)
  vf <- mvgk(D = DM)
  result <- vf(order = order, kappa = kappa)

  ## Covariance matrix =========================================================
  expectation <- rbind(cbind(sqrt(kappa[1] * kappa[1]) * exp(-kappa[13] * DM)[1,1],    sqrt(kappa[1] * kappa[2]) * exp(-kappa[13] * DM)[1,2],    sqrt(kappa[1] * kappa[3]) * exp(-kappa[13] * DM)[1,3],    sqrt(kappa[1] * kappa[4]) * kappa[10] * exp(-kappa[13] * DM)[1,1],    sqrt(kappa[1] * kappa[5]) * kappa[10] * exp(-kappa[13] * DM)[1,2],    sqrt(kappa[1] * kappa[6]) * kappa[10] * exp(-kappa[13] * DM)[1,3],    sqrt(kappa[1] * kappa[7]) * kappa[11] * exp(-kappa[13] * DM)[1,1],    sqrt(kappa[1] * kappa[8]) * kappa[11] * exp(-kappa[13] * DM)[1,2],    sqrt(kappa[1] * kappa[9]) * kappa[11] * exp(-kappa[13] * DM)[1,3]),
                       cbind(sqrt(kappa[2] * kappa[1]) * exp(-kappa[13] * DM)[2,1],    sqrt(kappa[2] * kappa[2]) * exp(-kappa[13] * DM)[2,2],    sqrt(kappa[2] * kappa[3]) * exp(-kappa[13] * DM)[2,3],    sqrt(kappa[2] * kappa[4]) * kappa[10] * exp(-kappa[13] * DM)[2,1],    sqrt(kappa[2] * kappa[5]) * kappa[10] * exp(-kappa[13] * DM)[2,2],    sqrt(kappa[2] * kappa[6]) * kappa[10] * exp(-kappa[13] * DM)[2,3],    sqrt(kappa[2] * kappa[7]) * kappa[11] * exp(-kappa[13] * DM)[2,1],    sqrt(kappa[2] * kappa[8]) * kappa[11] * exp(-kappa[13] * DM)[2,2],    sqrt(kappa[2] * kappa[9]) * kappa[11] * exp(-kappa[13] * DM)[2,3]),
                       cbind(sqrt(kappa[3] * kappa[1]) * exp(-kappa[13] * DM)[3,1],    sqrt(kappa[3] * kappa[2]) * exp(-kappa[13] * DM)[3,2],    sqrt(kappa[3] * kappa[3]) * exp(-kappa[13] * DM)[3,3],    sqrt(kappa[3] * kappa[4]) * kappa[10] * exp(-kappa[13] * DM)[3,1],    sqrt(kappa[3] * kappa[5]) * kappa[10] * exp(-kappa[13] * DM)[3,2],    sqrt(kappa[3] * kappa[6]) * kappa[10] * exp(-kappa[13] * DM)[3,3],    sqrt(kappa[3] * kappa[7]) * kappa[11] * exp(-kappa[13] * DM)[3,1],    sqrt(kappa[3] * kappa[8]) * kappa[11] * exp(-kappa[13] * DM)[3,2],    sqrt(kappa[3] * kappa[9]) * kappa[11] * exp(-kappa[13] * DM)[3,3]),

                       cbind(sqrt(kappa[4] * kappa[1]) * kappa[10] * exp(-kappa[13] * DM)[1,1],    sqrt(kappa[4] * kappa[2]) * kappa[10] * exp(-kappa[13] * DM)[1,2],    sqrt(kappa[4] * kappa[3]) * kappa[10] * exp(-kappa[13] * DM)[1,3],    sqrt(kappa[4] * kappa[4]) * exp(-kappa[13] * DM)[1,1],    sqrt(kappa[4] * kappa[5]) * exp(-kappa[13] * DM)[1,2],    sqrt(kappa[4] * kappa[6]) * exp(-kappa[13] * DM)[1,3],    sqrt(kappa[4] * kappa[7]) * kappa[12] * exp(-kappa[13] * DM)[1,1],    sqrt(kappa[4] * kappa[8]) * kappa[12] * exp(-kappa[13] * DM)[1,2],    sqrt(kappa[4] * kappa[9]) * kappa[12] * exp(-kappa[13] * DM)[1,3]),
                       cbind(sqrt(kappa[5] * kappa[1]) * kappa[10] * exp(-kappa[13] * DM)[2,1],    sqrt(kappa[5] * kappa[2]) * kappa[10] * exp(-kappa[13] * DM)[2,2],    sqrt(kappa[5] * kappa[3]) * kappa[10] * exp(-kappa[13] * DM)[2,3],    sqrt(kappa[5] * kappa[4]) * exp(-kappa[13] * DM)[2,1],    sqrt(kappa[5] * kappa[5]) * exp(-kappa[13] * DM)[2,2],    sqrt(kappa[5] * kappa[6]) * exp(-kappa[13] * DM)[2,3],    sqrt(kappa[5] * kappa[7]) * kappa[12] * exp(-kappa[13] * DM)[2,1],    sqrt(kappa[5] * kappa[8]) * kappa[12] * exp(-kappa[13] * DM)[2,2],    sqrt(kappa[5] * kappa[9]) * kappa[12] * exp(-kappa[13] * DM)[2,3]),
                       cbind(sqrt(kappa[6] * kappa[1]) * kappa[10] * exp(-kappa[13] * DM)[3,1],    sqrt(kappa[6] * kappa[2]) * kappa[10] * exp(-kappa[13] * DM)[3,2],    sqrt(kappa[6] * kappa[3]) * kappa[10] * exp(-kappa[13] * DM)[3,3],    sqrt(kappa[6] * kappa[4]) * exp(-kappa[13] * DM)[3,1],    sqrt(kappa[6] * kappa[5]) * exp(-kappa[13] * DM)[3,2],    sqrt(kappa[6] * kappa[6]) * exp(-kappa[13] * DM)[3,3],    sqrt(kappa[6] * kappa[7]) * kappa[12] * exp(-kappa[13] * DM)[3,1],    sqrt(kappa[6] * kappa[8]) * kappa[12] * exp(-kappa[13] * DM)[3,2],    sqrt(kappa[6] * kappa[9]) * kappa[12] * exp(-kappa[13] * DM)[3,3]),

                       cbind(sqrt(kappa[7] * kappa[1]) * kappa[11] * exp(-kappa[13] * DM)[1,1],    sqrt(kappa[7] * kappa[2]) * kappa[11] * exp(-kappa[13] * DM)[1,2],    sqrt(kappa[7] * kappa[3]) * kappa[11] * exp(-kappa[13] * DM)[1,3],    sqrt(kappa[7] * kappa[4]) * kappa[12] * exp(-kappa[13] * DM)[1,1],    sqrt(kappa[7] * kappa[5]) * kappa[12] * exp(-kappa[13] * DM)[1,2],    sqrt(kappa[7] * kappa[6]) * kappa[12] * exp(-kappa[13] * DM)[1,3],    sqrt(kappa[7] * kappa[7]) * exp(-kappa[13] * DM)[1,1],    sqrt(kappa[7] * kappa[8]) * exp(-kappa[13] * DM)[1,2],    sqrt(kappa[7] * kappa[9]) * exp(-kappa[13] * DM)[1,3]),
                       cbind(sqrt(kappa[8] * kappa[1]) * kappa[11] * exp(-kappa[13] * DM)[2,1],    sqrt(kappa[8] * kappa[2]) * kappa[11] * exp(-kappa[13] * DM)[2,2],    sqrt(kappa[8] * kappa[3]) * kappa[11] * exp(-kappa[13] * DM)[2,3],    sqrt(kappa[8] * kappa[4]) * kappa[12] * exp(-kappa[13] * DM)[2,1],    sqrt(kappa[8] * kappa[5]) * kappa[12] * exp(-kappa[13] * DM)[2,2],    sqrt(kappa[8] * kappa[6]) * kappa[12] * exp(-kappa[13] * DM)[2,3],    sqrt(kappa[8] * kappa[7]) * exp(-kappa[13] * DM)[2,1],    sqrt(kappa[8] * kappa[8]) * exp(-kappa[13] * DM)[2,2],    sqrt(kappa[8] * kappa[9]) * exp(-kappa[13] * DM)[2,3]),
                       cbind(sqrt(kappa[9] * kappa[1]) * kappa[11] * exp(-kappa[13] * DM)[3,1],    sqrt(kappa[9] * kappa[2]) * kappa[11] * exp(-kappa[13] * DM)[3,2],    sqrt(kappa[9] * kappa[3]) * kappa[11] * exp(-kappa[13] * DM)[3,3],    sqrt(kappa[9] * kappa[4]) * kappa[12] * exp(-kappa[13] * DM)[3,1],    sqrt(kappa[9] * kappa[5]) * kappa[12] * exp(-kappa[13] * DM)[3,2],    sqrt(kappa[9] * kappa[6]) * kappa[12] * exp(-kappa[13] * DM)[3,3],    sqrt(kappa[9] * kappa[7]) * exp(-kappa[13] * DM)[3,1],    sqrt(kappa[9] * kappa[8]) * exp(-kappa[13] * DM)[3,2],    sqrt(kappa[9] * kappa[9]) * exp(-kappa[13] * DM)[3,3]))

  expect_equal(result[[1]], expectation)

  ## Derivative wrt variance 1 =================================================
  expectation <- rbind(cbind(exp(-kappa[13] * DM)[1,1],    0.5 * sqrt(kappa[2]) * exp(-kappa[13] * DM)[1,2] / sqrt(kappa[1]),   0.5 * sqrt(kappa[3]) * exp(-kappa[13] * DM)[1,3] / sqrt(kappa[1]),    0.5 * sqrt(kappa[4]) * kappa[10] * exp(-kappa[13] * DM)[1,1] / sqrt(kappa[1]),    0.5 * sqrt(kappa[5]) * kappa[10] * exp(-kappa[13] * DM)[1,2] / sqrt(kappa[1]),    0.5 * sqrt(kappa[6]) * kappa[10] * exp(-kappa[13] * DM)[1,3] / sqrt(kappa[1]),    0.5 * sqrt(kappa[7]) * kappa[11] * exp(-kappa[13] * DM)[1,1] / sqrt(kappa[1]),    0.5 * sqrt(kappa[8]) * kappa[11] * exp(-kappa[13] * DM)[1,2] / sqrt(kappa[1]),    0.5 * sqrt(kappa[9]) * kappa[11] * exp(-kappa[13] * DM)[1,3] / sqrt(kappa[1])),
                       cbind(0.5 * sqrt(kappa[2]) * exp(-kappa[13] * DM)[2,1] / sqrt(kappa[1]), 0,    0,    0,    0,    0,    0,    0,    0),
                       cbind(0.5 * sqrt(kappa[3]) * exp(-kappa[13] * DM)[3,1] / sqrt(kappa[1]), 0,    0,    0,    0,    0,    0,    0,    0),
                       cbind(0.5 * sqrt(kappa[4]) * kappa[10] * exp(-kappa[13] * DM)[1,1] / sqrt(kappa[1]), 0,    0,    0,    0,    0,    0,    0,    0),
                       cbind(0.5 * sqrt(kappa[5]) * kappa[10] * exp(-kappa[13] * DM)[2,1] / sqrt(kappa[1]), 0,    0,    0,    0,    0,    0,    0,    0),
                       cbind(0.5 * sqrt(kappa[6]) * kappa[10] * exp(-kappa[13] * DM)[3,1] / sqrt(kappa[1]), 0,    0,    0,    0,    0,    0,    0,    0),
                       cbind(0.5 * sqrt(kappa[7]) * kappa[11] * exp(-kappa[13] * DM)[1,1] / sqrt(kappa[1]), 0,    0,    0,    0,    0,    0,    0,    0),
                       cbind(0.5 * sqrt(kappa[8]) * kappa[11] * exp(-kappa[13] * DM)[2,1] / sqrt(kappa[1]), 0,    0,    0,    0,    0,    0,    0,    0),
                       cbind(0.5 * sqrt(kappa[9]) * kappa[11] * exp(-kappa[13] * DM)[3,1] / sqrt(kappa[1]), 0,    0,    0,    0,    0,    0,    0,    0))

  expect_equal(result[[2]], expectation)

  ## Derivative wrt variance 2 =================================================
  expectation <- rbind(cbind(0,    0.5 * sqrt(kappa[1]) * exp(-kappa[13] * DM)[1,2] / sqrt(kappa[2]),    0,    0,    0,    0,    0,    0,    0),
                       cbind(0.5 * sqrt(kappa[1]) * exp(-kappa[13] * DM)[2,1] / sqrt(kappa[2]),    exp(-kappa[13] * DM)[2,2],    0.5 * sqrt(kappa[3]) * exp(-kappa[13] * DM)[2,3] / sqrt(kappa[2]),    0.5 * sqrt(kappa[4]) * kappa[10] * exp(-kappa[13] * DM)[2,1] / sqrt(kappa[2]),    0.5 * sqrt(kappa[5]) * kappa[10] * exp(-kappa[13] * DM)[2,2] / sqrt(kappa[2]),    0.5 * sqrt(kappa[6]) * kappa[10] * exp(-kappa[13] * DM)[2,3] / sqrt(kappa[2]),    0.5 * sqrt(kappa[7]) * kappa[11] * exp(-kappa[13] * DM)[2,1] / sqrt(kappa[2]),    0.5 * sqrt(kappa[8]) * kappa[11] * exp(-kappa[13] * DM)[2,2] / sqrt(kappa[2]),    0.5 * sqrt(kappa[9]) * kappa[11] * exp(-kappa[13] * DM)[2,3] / sqrt(kappa[2])),
                       cbind(0,    0.5 * sqrt(kappa[3]) * exp(-kappa[13] * DM)[3,2] / sqrt(kappa[2]),    0,    0,    0,    0,    0,    0,    0),
                       cbind(0,    0.5 * sqrt(kappa[4]) * kappa[10] * exp(-kappa[13] * DM)[1,2] / sqrt(kappa[2]),    0,    0,    0,    0,    0,    0,    0),
                       cbind(0,    0.5 * sqrt(kappa[5]) * kappa[10] * exp(-kappa[13] * DM)[2,2] / sqrt(kappa[2]),    0,    0,    0,    0,    0,    0,    0),
                       cbind(0,    0.5 * sqrt(kappa[6]) * kappa[10] * exp(-kappa[13] * DM)[3,2] / sqrt(kappa[2]),    0,    0,    0,    0,    0,    0,    0),
                       cbind(0,    0.5 * sqrt(kappa[7]) * kappa[11] * exp(-kappa[13] * DM)[1,2] / sqrt(kappa[2]),    0,    0,    0,    0,    0,    0,    0),
                       cbind(0,    0.5 * sqrt(kappa[8]) * kappa[11] * exp(-kappa[13] * DM)[2,2] / sqrt(kappa[2]),    0,    0,    0,    0,    0,    0,    0),
                       cbind(0,    0.5 * sqrt(kappa[9]) * kappa[11] * exp(-kappa[13] * DM)[3,2] / sqrt(kappa[2]),    0,    0,    0,    0,    0,    0,    0))

  expect_equal(result[[3]], expectation)

  ## Derivative wrt variance 3 =================================================
  expectation <- rbind(cbind(0,    0,    0.5 * sqrt(kappa[1]) * exp(-kappa[13] * DM)[1,3] / sqrt(kappa[3]),    0,    0,    0,    0,    0,    0),
                       cbind(0,    0,    0.5 * sqrt(kappa[2]) * exp(-kappa[13] * DM)[2,3] / sqrt(kappa[3]),    0,    0,    0,    0,    0,    0),
                       cbind(0.5 * sqrt(kappa[1]) * exp(-kappa[13] * DM)[3,1] / sqrt(kappa[3]),    0.5 * sqrt(kappa[2]) * exp(-kappa[13] * DM)[3,2] / sqrt(kappa[3]),    exp(-kappa[13] * DM)[3,3],    0.5 * sqrt(kappa[4]) * kappa[10] * exp(-kappa[13] * DM)[3,1] / sqrt(kappa[3]),    0.5 * sqrt(kappa[5]) * kappa[10] * exp(-kappa[13] * DM)[3,2] / sqrt(kappa[3]),    0.5 * sqrt(kappa[6]) * kappa[10] * exp(-kappa[13] * DM)[3,3] / sqrt(kappa[3]),    0.5 * sqrt(kappa[7]) * kappa[11] * exp(-kappa[13] * DM)[3,1] / sqrt(kappa[3]),    0.5 * sqrt(kappa[8]) * kappa[11] * exp(-kappa[13] * DM)[3,2] / sqrt(kappa[3]),    0.5 * sqrt(kappa[9]) * kappa[11] * exp(-kappa[13] * DM)[3,3] / sqrt(kappa[3])),
                       cbind(0,    0,    0.5 * sqrt(kappa[4]) * kappa[10] * exp(-kappa[13] * DM)[1,3] / sqrt(kappa[3]),    0,    0,    0,    0,    0,    0),
                       cbind(0,    0,    0.5 * sqrt(kappa[5]) * kappa[10] * exp(-kappa[13] * DM)[2,3] / sqrt(kappa[3]),    0,    0,    0,    0,    0,    0),
                       cbind(0,    0,    0.5 * sqrt(kappa[6]) * kappa[10] * exp(-kappa[13] * DM)[3,3] / sqrt(kappa[3]),    0,    0,    0,    0,    0,    0),
                       cbind(0,    0,    0.5 * sqrt(kappa[7]) * kappa[11] * exp(-kappa[13] * DM)[1,3] / sqrt(kappa[3]),    0,    0,    0,    0,    0,    0),
                       cbind(0,    0,    0.5 * sqrt(kappa[8]) * kappa[11] * exp(-kappa[13] * DM)[2,3] / sqrt(kappa[3]),    0,    0,    0,    0,    0,    0),
                       cbind(0,    0,    0.5 * sqrt(kappa[9]) * kappa[11] * exp(-kappa[13] * DM)[3,3] / sqrt(kappa[3]),    0,    0,    0,    0,    0,    0))

  expect_equal(result[[4]], expectation)

  ## Derivative wrt variance 4 =================================================
  expectation <- rbind(cbind(0,    0,    0,    0.5 * sqrt(kappa[1]) * kappa[10] * exp(-kappa[13] * DM)[1,1] / sqrt(kappa[4]),    0,    0,    0,    0,    0),
                       cbind(0,    0,    0,    0.5 * sqrt(kappa[2]) * kappa[10] * exp(-kappa[13] * DM)[2,1] / sqrt(kappa[4]),    0,    0,    0,    0,    0),
                       cbind(0,    0,    0,    0.5 * sqrt(kappa[3]) * kappa[10] * exp(-kappa[13] * DM)[3,1] / sqrt(kappa[4]),    0,    0,    0,    0,    0),
                       cbind(0.5 * sqrt(kappa[1]) * kappa[10] * exp(-kappa[13] * DM)[1,1] / sqrt(kappa[4]),    0.5 * sqrt(kappa[2]) * kappa[10] * exp(-kappa[13] * DM)[1,2] / sqrt(kappa[4]),    0.5 * sqrt(kappa[3]) * kappa[10] * exp(-kappa[13] * DM)[1,3] / sqrt(kappa[4]),    exp(-kappa[13] * DM)[1,1],    0.5 * sqrt(kappa[5]) * exp(-kappa[13] * DM)[1,2] / sqrt(kappa[4]),    0.5 * sqrt(kappa[6]) * exp(-kappa[13] * DM)[1,3] / sqrt(kappa[4]),    0.5 * sqrt(kappa[7]) * kappa[12] * exp(-kappa[13] * DM)[1,1] / sqrt(kappa[4]),    0.5 * sqrt(kappa[8]) * kappa[12] * exp(-kappa[13] * DM)[1,2] / sqrt(kappa[4]),    0.5 * sqrt(kappa[9]) * kappa[12] * exp(-kappa[13] * DM)[1,3] / sqrt(kappa[4])),
                       cbind(0,    0,    0,    0.5 * sqrt(kappa[5]) * exp(-kappa[13] * DM)[2,1] / sqrt(kappa[4]),    0,    0,    0,    0,    0),
                       cbind(0,    0,    0,    0.5 * sqrt(kappa[6]) * exp(-kappa[13] * DM)[3,1] / sqrt(kappa[4]),    0,    0,    0,    0,    0),
                       cbind(0,    0,    0,    0.5 * sqrt(kappa[7]) * kappa[12] * exp(-kappa[13] * DM)[1,1] / sqrt(kappa[4]),    0,    0,    0,    0,    0),
                       cbind(0,    0,    0,    0.5 * sqrt(kappa[8]) * kappa[12] * exp(-kappa[13] * DM)[2,1] / sqrt(kappa[4]),    0,    0,    0,    0,    0),
                       cbind(0,    0,    0,    0.5 * sqrt(kappa[9]) * kappa[12] * exp(-kappa[13] * DM)[3,1] / sqrt(kappa[4]),    0,    0,    0,    0,    0))

  expect_equal(result[[5]], expectation)

  ## Derivative wrt variance 5 =================================================
  expectation <- rbind(cbind(0,    0,    0,    0,    0.5 * sqrt(kappa[1]) * kappa[10] * exp(-kappa[13] * DM)[1,2] / sqrt(kappa[5]),    0,    0,    0,    0),
                       cbind(0,    0,    0,    0,    0.5 * sqrt(kappa[2]) * kappa[10] * exp(-kappa[13] * DM)[2,2] / sqrt(kappa[5]),    0,    0,    0,    0),
                       cbind(0,    0,    0,    0,    0.5 * sqrt(kappa[3]) * kappa[10] * exp(-kappa[13] * DM)[3,2] / sqrt(kappa[5]),    0,    0,    0,    0),
                       cbind(0,    0,    0,    0,    0.5 * sqrt(kappa[4]) * exp(-kappa[13] * DM)[1,2] / sqrt(kappa[5]),    0,    0,    0,    0),
                       cbind(0.5 * sqrt(kappa[1]) * kappa[10] * exp(-kappa[13] * DM)[2,1] / sqrt(kappa[5]),    0.5 * sqrt(kappa[2]) * kappa[10] * exp(-kappa[13] * DM)[2,2] / sqrt(kappa[5]),    0.5 * sqrt(kappa[3]) * kappa[10] * exp(-kappa[13] * DM)[2,3] / sqrt(kappa[5]),    0.5 * sqrt(kappa[4]) * exp(-kappa[13] * DM)[2,1] / sqrt(kappa[5]),    exp(-kappa[13] * DM)[2,2],    0.5 * sqrt(kappa[6]) * exp(-kappa[13] * DM)[2,3] / sqrt(kappa[5]),    0.5 * sqrt(kappa[7]) * kappa[12] * exp(-kappa[13] * DM)[2,1] / sqrt(kappa[5]),    0.5 * sqrt(kappa[8]) * kappa[12] * exp(-kappa[13] * DM)[2,2] / sqrt(kappa[5]),    0.5 * sqrt(kappa[9]) * kappa[12] * exp(-kappa[13] * DM)[2,3] / sqrt(kappa[5])),
                       cbind(0,    0,    0,    0,    0.5 * sqrt(kappa[6]) * exp(-kappa[13] * DM)[3,2] / sqrt(kappa[5]),    0,    0,    0,    0),
                       cbind(0,    0,    0,    0,    0.5 * sqrt(kappa[7]) * kappa[12] * exp(-kappa[13] * DM)[1,2] / sqrt(kappa[5]),    0,    0,    0,    0),
                       cbind(0,    0,    0,    0,    0.5 * sqrt(kappa[8]) * kappa[12] * exp(-kappa[13] * DM)[2,2] / sqrt(kappa[5]),    0,    0,    0,    0),
                       cbind(0,    0,    0,    0,    0.5 * sqrt(kappa[9]) * kappa[12] * exp(-kappa[13] * DM)[3,2] / sqrt(kappa[5]),    0,    0,    0,    0))

  expect_equal(result[[6]], expectation)

  ## Derivative wrt variance 6 =================================================
  expectation <- rbind(cbind(0,    0,    0,    0,    0,    0.5 * sqrt(kappa[1]) * kappa[10] * exp(-kappa[13] * DM)[1,3] / sqrt(kappa[6]),    0,    0,    0),
                       cbind(0,    0,    0,    0,    0,    0.5 * sqrt(kappa[2]) * kappa[10] * exp(-kappa[13] * DM)[2,3] / sqrt(kappa[6]),    0,    0,    0),
                       cbind(0,    0,    0,    0,    0,    0.5 * sqrt(kappa[3]) * kappa[10] * exp(-kappa[13] * DM)[3,3] / sqrt(kappa[6]),    0,    0,    0),
                       cbind(0,    0,    0,    0,    0,    0.5 * sqrt(kappa[4]) * exp(-kappa[13] * DM)[1,3] / sqrt(kappa[6]),    0,    0,    0),
                       cbind(0,    0,    0,    0,    0,    0.5 * sqrt(kappa[5]) * exp(-kappa[13] * DM)[2,3] / sqrt(kappa[6]),    0,    0,    0),
                       cbind(0.5 * sqrt(kappa[1]) * kappa[10] * exp(-kappa[13] * DM)[3,1] / sqrt(kappa[6]),    0.5 * sqrt(kappa[2]) * kappa[10] * exp(-kappa[13] * DM)[3,2] / sqrt(kappa[6]),    0.5 * sqrt(kappa[3]) * kappa[10] * exp(-kappa[13] * DM)[3,3] / sqrt(kappa[6]),    0.5 * sqrt(kappa[4]) * exp(-kappa[13] * DM)[3,1] / sqrt(kappa[6]),    0.5 * sqrt(kappa[5]) * exp(-kappa[13] * DM)[3,2] / sqrt(kappa[6]),    exp(-kappa[13] * DM)[3,3],    0.5 * sqrt(kappa[7]) * kappa[12] * exp(-kappa[13] * DM)[3,1] / sqrt(kappa[6]),    0.5 * sqrt(kappa[8]) * kappa[12] * exp(-kappa[13] * DM)[3,2] / sqrt(kappa[6]),    0.5 * sqrt(kappa[9]) * kappa[12] * exp(-kappa[13] * DM)[3,3] / sqrt(kappa[6])),
                       cbind(0,    0,    0,    0,    0,    0.5 * sqrt(kappa[7]) * kappa[12] * exp(-kappa[13] * DM)[1,3] / sqrt(kappa[6]),    0,    0,    0),
                       cbind(0,    0,    0,    0,    0,    0.5 * sqrt(kappa[8]) * kappa[12] * exp(-kappa[13] * DM)[2,3] / sqrt(kappa[6]),    0,    0,    0),
                       cbind(0,    0,    0,    0,    0,    0.5 * sqrt(kappa[9]) * kappa[12] * exp(-kappa[13] * DM)[3,3] / sqrt(kappa[6]),    0,    0,    0))

  expect_equal(result[[7]], expectation)

  ## Derivative wrt variance 7 =================================================
  expectation <- rbind(cbind(0,    0,    0,    0,    0,    0,    0.5 * sqrt(kappa[1]) * kappa[11] * exp(-kappa[13] * DM)[1,1] / sqrt(kappa[7]),    0,    0),
                       cbind(0,    0,    0,    0,    0,    0,    0.5 * sqrt(kappa[2]) * kappa[11] * exp(-kappa[13] * DM)[2,1] / sqrt(kappa[7]),    0,    0),
                       cbind(0,    0,    0,    0,    0,    0,    0.5 * sqrt(kappa[3]) * kappa[11] * exp(-kappa[13] * DM)[3,1] / sqrt(kappa[7]),    0,    0),
                       cbind(0,    0,    0,    0,    0,    0,    0.5 * sqrt(kappa[4]) * kappa[12] * exp(-kappa[13] * DM)[1,1] / sqrt(kappa[7]),    0,    0),
                       cbind(0,    0,    0,    0,    0,    0,    0.5 * sqrt(kappa[5]) * kappa[12] * exp(-kappa[13] * DM)[2,1] / sqrt(kappa[7]),    0,    0),
                       cbind(0,    0,    0,    0,    0,    0,    0.5 * sqrt(kappa[6]) * kappa[12] * exp(-kappa[13] * DM)[3,1] / sqrt(kappa[7]),    0,    0),
                       cbind(0.5 * sqrt(kappa[1]) * kappa[11] * exp(-kappa[13] * DM)[1,1] / sqrt(kappa[7]),    0.5 * sqrt(kappa[2]) * kappa[11] * exp(-kappa[13] * DM)[1,2] / sqrt(kappa[7]),    0.5 * sqrt(kappa[3]) * kappa[11] * exp(-kappa[13] * DM)[1,3] / sqrt(kappa[7]),    0.5 * sqrt(kappa[4]) * kappa[12] * exp(-kappa[13] * DM)[1,1] / sqrt(kappa[7]),    0.5 * sqrt(kappa[5]) * kappa[12] * exp(-kappa[13] * DM)[1,2] / sqrt(kappa[7]),    0.5 * sqrt(kappa[6]) * kappa[12] * exp(-kappa[13] * DM)[1,3] / sqrt(kappa[7]),    exp(-kappa[13] * DM)[1,1],    0.5 * sqrt(kappa[8]) * exp(-kappa[13] * DM)[1,2] / sqrt(kappa[7]),    0.5 * sqrt(kappa[9]) * exp(-kappa[13] * DM)[1,3] / sqrt(kappa[7])),
                       cbind(0,    0,    0,    0,    0,    0,    0.5 * sqrt(kappa[8]) * exp(-kappa[13] * DM)[2,1] / sqrt(kappa[7]),    0,    0),
                       cbind(0,    0,    0,    0,    0,    0,    0.5 * sqrt(kappa[9]) * exp(-kappa[13] * DM)[3,1] / sqrt(kappa[7]),    0,    0))

  expect_equal(result[[8]], expectation)

  ## Derivative wrt variance 8 =================================================
  expectation <- rbind(cbind(0,    0,    0,    0,    0,    0,    0,    0.5 * sqrt(kappa[1]) * kappa[11] * exp(-kappa[13] * DM)[1,2] / sqrt(kappa[8]),    0),
                       cbind(0,    0,    0,    0,    0,    0,    0,    0.5 * sqrt(kappa[2]) * kappa[11] * exp(-kappa[13] * DM)[2,2] / sqrt(kappa[8]),    0),
                       cbind(0,    0,    0,    0,    0,    0,    0,    0.5 * sqrt(kappa[3]) * kappa[11] * exp(-kappa[13] * DM)[3,2] / sqrt(kappa[8]),    0),
                       cbind(0,    0,    0,    0,    0,    0,    0,    0.5 * sqrt(kappa[4]) * kappa[12] * exp(-kappa[13] * DM)[1,2] / sqrt(kappa[8]),    0),
                       cbind(0,    0,    0,    0,    0,    0,    0,    0.5 * sqrt(kappa[5]) * kappa[12] * exp(-kappa[13] * DM)[2,2] / sqrt(kappa[8]),    0),
                       cbind(0,    0,    0,    0,    0,    0,    0,    0.5 * sqrt(kappa[6]) * kappa[12] * exp(-kappa[13] * DM)[3,2] / sqrt(kappa[8]),    0),
                       cbind(0,    0,    0,    0,    0,    0,    0,    0.5 * sqrt(kappa[7]) * exp(-kappa[13] * DM)[1,2] / sqrt(kappa[8]),    0),
                       cbind(0.5 * sqrt(kappa[1]) * kappa[11] * exp(-kappa[13] * DM)[2,1] / sqrt(kappa[8]),    0.5 * sqrt(kappa[2]) * kappa[11] * exp(-kappa[13] * DM)[2,2] / sqrt(kappa[8]),    0.5 * sqrt(kappa[3]) * kappa[11] * exp(-kappa[13] * DM)[2,3] / sqrt(kappa[8]),    0.5 * sqrt(kappa[4]) * kappa[12] * exp(-kappa[13] * DM)[2,1] / sqrt(kappa[8]),    0.5 * sqrt(kappa[5]) * kappa[12] * exp(-kappa[13] * DM)[2,2] / sqrt(kappa[8]),    0.5 * sqrt(kappa[6]) * kappa[12] * exp(-kappa[13] * DM)[2,3] / sqrt(kappa[8]),    0.5 * sqrt(kappa[7]) * exp(-kappa[13] * DM)[2,1] / sqrt(kappa[8]),    exp(-kappa[13] * DM)[2,2],    0.5 * sqrt(kappa[9]) * exp(-kappa[13] * DM)[2,3] / sqrt(kappa[8])),
                       cbind(0,    0,    0,    0,    0,    0,    0,    0.5 * sqrt(kappa[9]) * exp(-kappa[13] * DM)[3,2] / sqrt(kappa[8]),    0))

  expect_equal(result[[9]], expectation)

  ## Derivative wrt variance 9 =================================================
  expectation <- rbind(cbind(0,    0,    0,    0,    0,    0,    0,    0,    0.5 * sqrt(kappa[1]) * kappa[11] * exp(-kappa[13] * DM)[1,3] / sqrt(kappa[9])),
                       cbind(0,    0,    0,    0,    0,    0,    0,    0,    0.5 * sqrt(kappa[2]) * kappa[11] * exp(-kappa[13] * DM)[2,3] / sqrt(kappa[9])),
                       cbind(0,    0,    0,    0,    0,    0,    0,    0,    0.5 * sqrt(kappa[3]) * kappa[11] * exp(-kappa[13] * DM)[3,3] / sqrt(kappa[9])),
                       cbind(0,    0,    0,    0,    0,    0,    0,    0,    0.5 * sqrt(kappa[4]) * kappa[12] * exp(-kappa[13] * DM)[1,3] / sqrt(kappa[9])),
                       cbind(0,    0,    0,    0,    0,    0,    0,    0,    0.5 * sqrt(kappa[5]) * kappa[12] * exp(-kappa[13] * DM)[2,3] / sqrt(kappa[9])),
                       cbind(0,    0,    0,    0,    0,    0,    0,    0,    0.5 * sqrt(kappa[6]) * kappa[12] * exp(-kappa[13] * DM)[3,3] / sqrt(kappa[9])),
                       cbind(0,    0,    0,    0,    0,    0,    0,    0,    0.5 * sqrt(kappa[7]) * exp(-kappa[13] * DM)[1,3] / sqrt(kappa[9])),
                       cbind(0,    0,    0,    0,    0,    0,    0,    0,    0.5 * sqrt(kappa[8]) * exp(-kappa[13] * DM)[2,3] / sqrt(kappa[9])),
                       cbind(0.5 * sqrt(kappa[1]) * kappa[11] * exp(-kappa[13] * DM)[3,1] / sqrt(kappa[9]),    0.5 * sqrt(kappa[2]) * kappa[11] * exp(-kappa[13] * DM)[3,2] / sqrt(kappa[9]),    0.5 * sqrt(kappa[3]) * kappa[11] * exp(-kappa[13] * DM)[3,3] / sqrt(kappa[9]),    0.5 * sqrt(kappa[4]) * kappa[12] * exp(-kappa[13] * DM)[3,1] / sqrt(kappa[9]),    0.5 * sqrt(kappa[5]) * kappa[12] * exp(-kappa[13] * DM)[3,2] / sqrt(kappa[9]),    0.5 * sqrt(kappa[6]) * kappa[12] * exp(-kappa[13] * DM)[3,3] / sqrt(kappa[9]),    0.5 * sqrt(kappa[7]) * exp(-kappa[13] * DM)[3,1] / sqrt(kappa[9]),    0.5 * sqrt(kappa[8]) * exp(-kappa[13] * DM)[3,2] / sqrt(kappa[9]),    exp(-kappa[13] * DM)[3,3]))

  expect_equal(result[[10]], expectation)

  ## Derivative wrt correlation 1-2 ============================================
  expectation <- rbind(cbind(0,    0,    0,    sqrt(kappa[1] * kappa[4]) * exp(-kappa[13] * DM)[1,1],    sqrt(kappa[1] * kappa[5]) * exp(-kappa[13] * DM)[1,2],    sqrt(kappa[1] * kappa[6]) * exp(-kappa[13] * DM)[1,3],    0,    0,    0),
                       cbind(0,    0,    0,    sqrt(kappa[2] * kappa[4]) * exp(-kappa[13] * DM)[2,1],    sqrt(kappa[2] * kappa[5]) * exp(-kappa[13] * DM)[2,2],    sqrt(kappa[2] * kappa[6]) * exp(-kappa[13] * DM)[2,3],    0,    0,    0),
                       cbind(0,    0,    0,    sqrt(kappa[3] * kappa[4]) * exp(-kappa[13] * DM)[3,1],    sqrt(kappa[3] * kappa[5]) * exp(-kappa[13] * DM)[3,2],    sqrt(kappa[3] * kappa[6]) * exp(-kappa[13] * DM)[3,3],    0,    0,    0),
                       cbind(sqrt(kappa[4] * kappa[1]) * exp(-kappa[13] * DM)[1,1],    sqrt(kappa[4] * kappa[2]) * exp(-kappa[13] * DM)[1,2],    sqrt(kappa[4] * kappa[3]) * exp(-kappa[13] * DM)[1,3],    0,    0,    0,    0,    0,    0),
                       cbind(sqrt(kappa[5] * kappa[1]) * exp(-kappa[13] * DM)[2,1],    sqrt(kappa[5] * kappa[2]) * exp(-kappa[13] * DM)[2,2],    sqrt(kappa[5] * kappa[3]) * exp(-kappa[13] * DM)[2,3],    0,    0,    0,    0,    0,    0),
                       cbind(sqrt(kappa[6] * kappa[1]) * exp(-kappa[13] * DM)[3,1],    sqrt(kappa[6] * kappa[2]) * exp(-kappa[13] * DM)[3,2],    sqrt(kappa[6] * kappa[3]) * exp(-kappa[13] * DM)[3,3],    0,    0,    0,    0,    0,    0),
                       cbind(0,    0,    0,    0,    0,    0,    0,    0,    0),
                       cbind(0,    0,    0,    0,    0,    0,    0,    0,    0),
                       cbind(0,    0,    0,    0,    0,    0,    0,    0,    0))

  expect_equal(result[[11]], expectation)

  ## Derivative wrt correlation 1-3 ============================================
  expectation <- rbind(cbind(0,    0,    0,    0,    0,    0,    sqrt(kappa[1] * kappa[7]) * exp(-kappa[13] * DM)[1,1],    sqrt(kappa[1] * kappa[8]) * exp(-kappa[13] * DM)[1,2],    sqrt(kappa[1] * kappa[9]) * exp(-kappa[13] * DM)[1,3]),
                       cbind(0,    0,    0,    0,    0,    0,    sqrt(kappa[2] * kappa[7]) * exp(-kappa[13] * DM)[2,1],    sqrt(kappa[2] * kappa[8]) * exp(-kappa[13] * DM)[2,2],    sqrt(kappa[2] * kappa[9]) * exp(-kappa[13] * DM)[2,3]),
                       cbind(0,    0,    0,    0,    0,    0,    sqrt(kappa[3] * kappa[7]) * exp(-kappa[13] * DM)[3,1],    sqrt(kappa[3] * kappa[8]) * exp(-kappa[13] * DM)[3,2],    sqrt(kappa[3] * kappa[9]) * exp(-kappa[13] * DM)[3,3]),
                       cbind(0,    0,    0,    0,    0,    0,    0,    0,    0),
                       cbind(0,    0,    0,    0,    0,    0,    0,    0,    0),
                       cbind(0,    0,    0,    0,    0,    0,    0,    0,    0),
                       cbind(sqrt(kappa[7] * kappa[1]) * exp(-kappa[13] * DM)[1,1],    sqrt(kappa[7] * kappa[2]) * exp(-kappa[13] * DM)[1,2],    sqrt(kappa[7] * kappa[3]) * exp(-kappa[13] * DM)[1,3],    0,    0,    0,    0,    0,    0),
                       cbind(sqrt(kappa[8] * kappa[1]) * exp(-kappa[13] * DM)[2,1],    sqrt(kappa[8] * kappa[2]) * exp(-kappa[13] * DM)[2,2],    sqrt(kappa[8] * kappa[3]) * exp(-kappa[13] * DM)[2,3],    0,    0,    0,    0,    0,    0),
                       cbind(sqrt(kappa[9] * kappa[1]) * exp(-kappa[13] * DM)[3,1],    sqrt(kappa[9] * kappa[2]) * exp(-kappa[13] * DM)[3,2],    sqrt(kappa[9] * kappa[3]) * exp(-kappa[13] * DM)[3,3],    0,    0,    0,    0,    0,    0))

  expect_equal(result[[12]], expectation)

  ## Derivative wrt correlation 2-3 ============================================
  expectation <- rbind(cbind(0,    0,    0,    0,    0,    0,    0,    0,    0),
                       cbind(0,    0,    0,    0,    0,    0,    0,    0,    0),
                       cbind(0,    0,    0,    0,    0,    0,    0,    0,    0),
                       cbind(0,    0,    0,    0,    0,    0,    sqrt(kappa[4] * kappa[7]) * exp(-kappa[13] * DM)[1,1],    sqrt(kappa[4] * kappa[8]) * exp(-kappa[13] * DM)[1,2],    sqrt(kappa[4] * kappa[9]) * exp(-kappa[13] * DM)[1,3]),
                       cbind(0,    0,    0,    0,    0,    0,    sqrt(kappa[5] * kappa[7]) * exp(-kappa[13] * DM)[2,1],    sqrt(kappa[5] * kappa[8]) * exp(-kappa[13] * DM)[2,2],    sqrt(kappa[5] * kappa[9]) * exp(-kappa[13] * DM)[2,3]),
                       cbind(0,    0,    0,    0,    0,    0,    sqrt(kappa[6] * kappa[7]) * exp(-kappa[13] * DM)[3,1],    sqrt(kappa[6] * kappa[8]) * exp(-kappa[13] * DM)[3,2],    sqrt(kappa[6] * kappa[9]) * exp(-kappa[13] * DM)[3,3]),
                       cbind(0,    0,    0,    sqrt(kappa[7] * kappa[4]) * exp(-kappa[13] * DM)[1,1],    sqrt(kappa[7] * kappa[5]) * exp(-kappa[13] * DM)[1,2],    sqrt(kappa[7] * kappa[6]) * exp(-kappa[13] * DM)[1,3],    0,    0,    0),
                       cbind(0,    0,    0,    sqrt(kappa[8] * kappa[4]) * exp(-kappa[13] * DM)[2,1],    sqrt(kappa[8] * kappa[5]) * exp(-kappa[13] * DM)[2,2],    sqrt(kappa[8] * kappa[6]) * exp(-kappa[13] * DM)[2,3],    0,    0,    0),
                       cbind(0,    0,    0,    sqrt(kappa[9] * kappa[4]) * exp(-kappa[13] * DM)[3,1],    sqrt(kappa[9] * kappa[5]) * exp(-kappa[13] * DM)[3,2],    sqrt(kappa[9] * kappa[6]) * exp(-kappa[13] * DM)[3,3],    0,    0,    0))

  expect_equal(result[[13]], expectation)

  ## Derivative wrt the bandwidth ==============================================
  expectation <- rbind(cbind(-DM[1,1] * sqrt(kappa[1] * kappa[1]) * exp(-kappa[13] * DM)[1,1],    -DM[1,2] * sqrt(kappa[1] * kappa[2]) * exp(-kappa[13] * DM)[1,2],    -DM[1,3] * sqrt(kappa[1] * kappa[3]) * exp(-kappa[13] * DM)[1,3],    -DM[1,1] * sqrt(kappa[1] * kappa[4]) * kappa[10] * exp(-kappa[13] * DM)[1,1],    -DM[1,2] * sqrt(kappa[1] * kappa[5]) * kappa[10] * exp(-kappa[13] * DM)[1,2],    -DM[1,3] * sqrt(kappa[1] * kappa[6]) * kappa[10] * exp(-kappa[13] * DM)[1,3],    -DM[1,1] * sqrt(kappa[1] * kappa[7]) * kappa[11] * exp(-kappa[13] * DM)[1,1],    -DM[1,2] * sqrt(kappa[1] * kappa[8]) * kappa[11] * exp(-kappa[13] * DM)[1,2],    -DM[1,3] * sqrt(kappa[1] * kappa[9]) * kappa[11] * exp(-kappa[13] * DM)[1,3]),
                       cbind(-DM[2,1] * sqrt(kappa[2] * kappa[1]) * exp(-kappa[13] * DM)[2,1],    -DM[2,2] * sqrt(kappa[2] * kappa[2]) * exp(-kappa[13] * DM)[2,2],    -DM[2,3] * sqrt(kappa[2] * kappa[3]) * exp(-kappa[13] * DM)[2,3],    -DM[2,1] * sqrt(kappa[2] * kappa[4]) * kappa[10] * exp(-kappa[13] * DM)[2,1],    -DM[2,2] * sqrt(kappa[2] * kappa[5]) * kappa[10] * exp(-kappa[13] * DM)[2,2],    -DM[2,3] * sqrt(kappa[2] * kappa[6]) * kappa[10] * exp(-kappa[13] * DM)[2,3],    -DM[2,1] * sqrt(kappa[2] * kappa[7]) * kappa[11] * exp(-kappa[13] * DM)[2,1],    -DM[2,2] * sqrt(kappa[2] * kappa[8]) * kappa[11] * exp(-kappa[13] * DM)[2,2],    -DM[2,3] * sqrt(kappa[2] * kappa[9]) * kappa[11] * exp(-kappa[13] * DM)[2,3]),
                       cbind(-DM[3,1] * sqrt(kappa[3] * kappa[1]) * exp(-kappa[13] * DM)[3,1],    -DM[3,2] * sqrt(kappa[3] * kappa[2]) * exp(-kappa[13] * DM)[3,2],    -DM[3,3] * sqrt(kappa[3] * kappa[3]) * exp(-kappa[13] * DM)[3,3],    -DM[3,1] * sqrt(kappa[3] * kappa[4]) * kappa[10] * exp(-kappa[13] * DM)[3,1],    -DM[3,2] * sqrt(kappa[3] * kappa[5]) * kappa[10] * exp(-kappa[13] * DM)[3,2],    -DM[3,3] * sqrt(kappa[3] * kappa[6]) * kappa[10] * exp(-kappa[13] * DM)[3,3],    -DM[3,1] * sqrt(kappa[3] * kappa[7]) * kappa[11] * exp(-kappa[13] * DM)[3,1],    -DM[3,2] * sqrt(kappa[3] * kappa[8]) * kappa[11] * exp(-kappa[13] * DM)[3,2],    -DM[3,3] * sqrt(kappa[3] * kappa[9]) * kappa[11] * exp(-kappa[13] * DM)[3,3]),

                       cbind(-DM[1,1] * sqrt(kappa[4] * kappa[1]) * kappa[10] * exp(-kappa[13] * DM)[1,1],    -DM[1,2] * sqrt(kappa[4] * kappa[2]) * kappa[10] * exp(-kappa[13] * DM)[1,2],    -DM[1,3] * sqrt(kappa[4] * kappa[3]) * kappa[10] * exp(-kappa[13] * DM)[1,3],    -DM[1,1] * sqrt(kappa[4] * kappa[4]) * exp(-kappa[13] * DM)[1,1],    -DM[1,2] * sqrt(kappa[4] * kappa[5]) * exp(-kappa[13] * DM)[1,2],    -DM[1,3] * sqrt(kappa[4] * kappa[6]) * exp(-kappa[13] * DM)[1,3],    -DM[1,1] * sqrt(kappa[4] * kappa[7]) * kappa[12] * exp(-kappa[13] * DM)[1,1],    -DM[1,2] * sqrt(kappa[4] * kappa[8]) * kappa[12] * exp(-kappa[13] * DM)[1,2],    -DM[1,3] * sqrt(kappa[4] * kappa[9]) * kappa[12] * exp(-kappa[13] * DM)[1,3]),
                       cbind(-DM[2,1] * sqrt(kappa[5] * kappa[1]) * kappa[10] * exp(-kappa[13] * DM)[2,1],    -DM[2,2] * sqrt(kappa[5] * kappa[2]) * kappa[10] * exp(-kappa[13] * DM)[2,2],    -DM[2,3] * sqrt(kappa[5] * kappa[3]) * kappa[10] * exp(-kappa[13] * DM)[2,3],    -DM[2,1] * sqrt(kappa[5] * kappa[4]) * exp(-kappa[13] * DM)[2,1],    -DM[2,2] * sqrt(kappa[5] * kappa[5]) * exp(-kappa[13] * DM)[2,2],    -DM[2,3] * sqrt(kappa[5] * kappa[6]) * exp(-kappa[13] * DM)[2,3],    -DM[2,1] * sqrt(kappa[5] * kappa[7]) * kappa[12] * exp(-kappa[13] * DM)[2,1],    -DM[2,2] * sqrt(kappa[5] * kappa[8]) * kappa[12] * exp(-kappa[13] * DM)[2,2],    -DM[2,3] * sqrt(kappa[5] * kappa[9]) * kappa[12] * exp(-kappa[13] * DM)[2,3]),
                       cbind(-DM[3,1] * sqrt(kappa[6] * kappa[1]) * kappa[10] * exp(-kappa[13] * DM)[3,1],    -DM[3,2] * sqrt(kappa[6] * kappa[2]) * kappa[10] * exp(-kappa[13] * DM)[3,2],    -DM[3,3] * sqrt(kappa[6] * kappa[3]) * kappa[10] * exp(-kappa[13] * DM)[3,3],    -DM[3,1] * sqrt(kappa[6] * kappa[4]) * exp(-kappa[13] * DM)[3,1],    -DM[3,2] * sqrt(kappa[6] * kappa[5]) * exp(-kappa[13] * DM)[3,2],    -DM[3,3] * sqrt(kappa[6] * kappa[6]) * exp(-kappa[13] * DM)[3,3],    -DM[3,1] * sqrt(kappa[6] * kappa[7]) * kappa[12] * exp(-kappa[13] * DM)[3,1],    -DM[3,2] * sqrt(kappa[6] * kappa[8]) * kappa[12] * exp(-kappa[13] * DM)[3,2],    -DM[3,3] * sqrt(kappa[6] * kappa[9]) * kappa[12] * exp(-kappa[13] * DM)[3,3]),

                       cbind(-DM[1,1] * sqrt(kappa[7] * kappa[1]) * kappa[11] * exp(-kappa[13] * DM)[1,1],    -DM[1,2] * sqrt(kappa[7] * kappa[2]) * kappa[11] * exp(-kappa[13] * DM)[1,2],    -DM[1,3] * sqrt(kappa[7] * kappa[3]) * kappa[11] * exp(-kappa[13] * DM)[1,3],    -DM[1,1] * sqrt(kappa[7] * kappa[4]) * kappa[12] * exp(-kappa[13] * DM)[1,1],    -DM[1,2] * sqrt(kappa[7] * kappa[5]) * kappa[12] * exp(-kappa[13] * DM)[1,2],    -DM[1,3] * sqrt(kappa[7] * kappa[6]) * kappa[12] * exp(-kappa[13] * DM)[1,3],    -DM[1,1] * sqrt(kappa[7] * kappa[7]) * exp(-kappa[13] * DM)[1,1],    -DM[1,2] * sqrt(kappa[7] * kappa[8]) * exp(-kappa[13] * DM)[1,2],    -DM[1,3] * sqrt(kappa[7] * kappa[9]) * exp(-kappa[13] * DM)[1,3]),
                       cbind(-DM[2,1] * sqrt(kappa[8] * kappa[1]) * kappa[11] * exp(-kappa[13] * DM)[2,1],    -DM[2,2] * sqrt(kappa[8] * kappa[2]) * kappa[11] * exp(-kappa[13] * DM)[2,2],    -DM[2,3] * sqrt(kappa[8] * kappa[3]) * kappa[11] * exp(-kappa[13] * DM)[2,3],    -DM[2,1] * sqrt(kappa[8] * kappa[4]) * kappa[12] * exp(-kappa[13] * DM)[2,1],    -DM[2,2] * sqrt(kappa[8] * kappa[5]) * kappa[12] * exp(-kappa[13] * DM)[2,2],    -DM[2,3] * sqrt(kappa[8] * kappa[6]) * kappa[12] * exp(-kappa[13] * DM)[2,3],    -DM[2,1] * sqrt(kappa[8] * kappa[7]) * exp(-kappa[13] * DM)[2,1],    -DM[2,2] * sqrt(kappa[8] * kappa[8]) * exp(-kappa[13] * DM)[2,2],    -DM[2,3] * sqrt(kappa[8] * kappa[9]) * exp(-kappa[13] * DM)[2,3]),
                       cbind(-DM[3,1] * sqrt(kappa[9] * kappa[1]) * kappa[11] * exp(-kappa[13] * DM)[3,1],    -DM[3,2] * sqrt(kappa[9] * kappa[2]) * kappa[11] * exp(-kappa[13] * DM)[3,2],    -DM[3,3] * sqrt(kappa[9] * kappa[3]) * kappa[11] * exp(-kappa[13] * DM)[3,3],    -DM[3,1] * sqrt(kappa[9] * kappa[4]) * kappa[12] * exp(-kappa[13] * DM)[3,1],    -DM[3,2] * sqrt(kappa[9] * kappa[5]) * kappa[12] * exp(-kappa[13] * DM)[3,2],    -DM[3,3] * sqrt(kappa[9] * kappa[6]) * kappa[12] * exp(-kappa[13] * DM)[3,3],    -DM[3,1] * sqrt(kappa[9] * kappa[7]) * exp(-kappa[13] * DM)[3,1],    -DM[3,2] * sqrt(kappa[9] * kappa[8]) * exp(-kappa[13] * DM)[3,2],    -DM[3,3] * sqrt(kappa[9] * kappa[9]) * exp(-kappa[13] * DM)[3,3]))

  expect_equal(result[[14]], expectation)

  # Testing the GxE case =======================================================
  LK <- matrix(c(1.0, 0.8, 0.7,
                 0.8, 1.0, 0.9,
                 0.7, 0.9, 1.0),
               nrow = 3, byrow = TRUE)
  DM <- log(LK) / -0.15

  order <- 3
  kappa <- c(2, 3, 4, 0.1)
  vf <- mvgk(D = DM)
  result <- vf(order = order, kappa = kappa)

  ## Covariance matrix =========================================================
  expectation <- rbind(cbind(sqrt(kappa[1] * kappa[1]) * exp(-kappa[4] * DM)[1,1],    sqrt(kappa[1] * kappa[2]) * exp(-kappa[4] * DM)[1,2],    sqrt(kappa[1] * kappa[3]) * exp(-kappa[4] * DM)[1,3]),
                       cbind(sqrt(kappa[2] * kappa[1]) * exp(-kappa[4] * DM)[2,1],    sqrt(kappa[2] * kappa[2]) * exp(-kappa[4] * DM)[2,2],    sqrt(kappa[2] * kappa[3]) * exp(-kappa[4] * DM)[2,3]),
                       cbind(sqrt(kappa[3] * kappa[1]) * exp(-kappa[4] * DM)[3,1],    sqrt(kappa[3] * kappa[2]) * exp(-kappa[4] * DM)[3,2],    sqrt(kappa[3] * kappa[3]) * exp(-kappa[4] * DM)[3,3]))

  expect_equal(result[[1]], expectation)

  ## Derivative wrt variance 1 =================================================
  expectation <- rbind(cbind(exp(-kappa[4] * DM)[1,1],    0.5 * sqrt(kappa[2]) * exp(-kappa[4] * DM)[1,2] / sqrt(kappa[1]),    0.5 * sqrt(kappa[3]) * exp(-kappa[4] * DM)[1,3] / sqrt(kappa[1])),
                       cbind(0.5 * sqrt(kappa[2]) * exp(-kappa[4] * DM)[2,1] / sqrt(kappa[1]),    0,    0),
                       cbind(0.5 * sqrt(kappa[3]) * exp(-kappa[4] * DM)[3,1] / sqrt(kappa[1]),    0,    0))

  expect_equal(result[[2]], expectation)

  ## Derivative wrt variance 2 =================================================
  expectation <- rbind(cbind(0,    0.5 * sqrt(kappa[1]) * exp(-kappa[4] * DM)[1,2] / sqrt(kappa[2]),    0),
                       cbind(0.5 * sqrt(kappa[1]) * exp(-kappa[4] * DM)[2,1] / sqrt(kappa[2]),    exp(-kappa[4] * DM)[2,2],    0.5 * sqrt(kappa[3]) * exp(-kappa[4] * DM)[2,3] / sqrt(kappa[2])),
                       cbind(0,    0.5 * sqrt(kappa[3]) * exp(-kappa[4] * DM)[3,2] / sqrt(kappa[2]),    0))

  expect_equal(result[[3]], expectation)

  ## Derivative wrt variance 3 =================================================
  expectation <- rbind(cbind(0,    0,    0.5 * sqrt(kappa[1]) * exp(-kappa[4] * DM)[1,3] / sqrt(kappa[3])),
                       cbind(0,    0,    0.5 * sqrt(kappa[2]) * exp(-kappa[4] * DM)[2,3] / sqrt(kappa[3])),
                       cbind(0.5 * sqrt(kappa[1]) * exp(-kappa[4] * DM)[3,1] / sqrt(kappa[3]),    0.5 * sqrt(kappa[2]) * exp(-kappa[4] * DM)[3,2] / sqrt(kappa[3]),    exp(-kappa[4] * DM)[3,3]))

  expect_equal(result[[4]], expectation)

  ## Derivative wrt the bandwidth ==============================================
  expectation <- rbind(cbind(-DM[1,1] * sqrt(kappa[1] * kappa[1]) * exp(-kappa[4] * DM)[1,1],    -DM[1,2] * sqrt(kappa[1] * kappa[2]) * exp(-kappa[4] * DM)[1,2],    -DM[1,3] * sqrt(kappa[1] * kappa[3]) * exp(-kappa[4] * DM)[1,3]),
                       cbind(-DM[2,1] * sqrt(kappa[2] * kappa[1]) * exp(-kappa[4] * DM)[2,1],    -DM[2,2] * sqrt(kappa[2] * kappa[2]) * exp(-kappa[4] * DM)[2,2],    -DM[2,3] * sqrt(kappa[2] * kappa[3]) * exp(-kappa[4] * DM)[2,3]),
                       cbind(-DM[3,1] * sqrt(kappa[3] * kappa[1]) * exp(-kappa[4] * DM)[3,1],    -DM[3,2] * sqrt(kappa[3] * kappa[2]) * exp(-kappa[4] * DM)[3,2],    -DM[3,3] * sqrt(kappa[3] * kappa[3]) * exp(-kappa[4] * DM)[3,3]))

  expect_equal(result[[5]], expectation)
})


