test_that("agree_test detects a strong gradient and returns finite output", {
  d <- simulate_agree(n_schools = 30, n_per_school = 20, gamma = 0.8,
                      reliability = 1, n_pv = 1, seed = 6)
  x <- ca_data(d, "first", "second", "item", "PV1", cluster = "school", id = "student")
  at <- agree_test(x, adjust = c("none", "score"), unit = 1)
  expect_s3_class(at, "agree_test")
  expect_equal(nrow(at$results), 2L)
  expect_lt(at$results$trend_p[1], 0.001)
  expect_gt(at$results$slope_sd[1], 0)
  at1 <- agree_test(x, df = 1)
  expect_equal(at1$results$omnibus_p, at1$results$trend_p)
})

test_that("severity_test detects a harsh second rater", {
  d <- simulate_agree(n_schools = 30, n_per_school = 20, severity = 2,
                      alpha = c(zero = -2.5, partial = -1.5, full = -2.5),
                      reliability = 1, n_pv = 1, seed = 7)
  x <- ca_data(d, "first", "second", "item", "PV1", cluster = "school", id = "student")
  st <- severity_test(x)
  expect_gt(st$prop_second_lower, 0.5)
  expect_lt(st$joint["p"], 0.01)
})

test_that("rater_model scores match numerical derivatives and recover gamma", {
  d <- simulate_agree(n_schools = 40, n_per_school = 25, gamma = 0.5,
                      alpha = c(zero = -4, partial = -2, full = -3),
                      reliability = 1, n_pv = 1, seed = 8)
  x <- ca_data(d, "first", "second", "item", "PV1", cluster = "school", id = "student")
  S <- condagree:::.rm_setup(x)
  z <- as.numeric(scale(x$ability[, 1]))
  set.seed(1)
  p <- condagree:::.rm_start(S, z) + rnorm(S$P, 0, 0.05)
  an <- colSums(condagree:::.rm_core(S, p, z)$scores)
  nu <- vapply(seq_along(p), function(k) {
    h <- 1e-5
    e <- replace(numeric(length(p)), k, h)
    (sum(condagree:::.rm_core(S, p + e, z, FALSE)$ll) -
        sum(condagree:::.rm_core(S, p - e, z, FALSE)$ll)) / (2 * h)
  }, numeric(1))
  expect_lt(max(abs(an - nu)), 1e-4)
  rm <- rater_model(x)
  expect_true(all(rm$converged))
  g <- rm$coefficients[rm$coefficients$parameter == "gamma", ]
  expect_true(g$lower < 0.5 && g$upper > 0.5)
  expect_equal(nrow(rm$fit_check), 5L)
  expect_invisible(plot(rm))
})

test_that("rater_model handles dichotomous items only", {
  d <- simulate_agree(n_schools = 20, n_per_school = 20, max_scores = rep(1, 6),
                      alpha = c(zero = -3, partial = -3, full = -3),
                      reliability = 1, n_pv = 1, seed = 9)
  x <- ca_data(d, "first", "second", "item", "PV1", cluster = "school", id = "student")
  rm <- rater_model(x)
  expect_false("severity" %in% rm$coefficients$parameter[!is.na(rm$coefficients$estimate)])
})

test_that("item-specific error intercepts: exact scores and estimation", {
  d <- simulate_agree(n_schools = 30, n_per_school = 20, gamma = 0.3,
                      alpha = c(zero = -4, partial = -2, full = -3),
                      reliability = 1, n_pv = 1, seed = 10)
  x <- ca_data(d, "first", "second", "item", "PV1", cluster = "school", id = "student")
  for (ab in c("maxscore", "item")) {
    S <- condagree:::.rm_setup(x, ab)
    z <- as.numeric(scale(x$ability[, 1]))
    set.seed(2)
    p <- condagree:::.rm_start(S, z) + rnorm(S$P, 0, 0.05)
    an <- colSums(condagree:::.rm_core(S, p, z)$scores)
    nu <- vapply(seq_along(p), function(k) {
      h <- 1e-5
      e <- replace(numeric(length(p)), k, h)
      (sum(condagree:::.rm_core(S, p + e, z, FALSE)$ll) -
          sum(condagree:::.rm_core(S, p - e, z, FALSE)$ll)) / (2 * h)
    }, numeric(1))
    expect_lt(max(abs(an - nu)), 1e-4)
  }
  rm <- rater_model(x, alpha_by = "item")
  expect_true(all(rm$converged))
  expect_true("gamma" %in% rm$coefficients$parameter)
  expect_gt(nrow(rm$alpha), 3L)
})
