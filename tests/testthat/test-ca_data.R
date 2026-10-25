test_that("ca_data validates input and stores structure", {
  d <- simulate_agree(n_schools = 5, n_per_school = 6, seed = 1)
  x <- ca_data(d, "first", "second", "item", c("PV1", "PV2"),
               cluster = "school", id = "student")
  expect_s3_class(x, "ca_data")
  expect_equal(x$n, nrow(d))
  expect_equal(x$M, 2L)
  expect_equal(nlevels(x$cluster), 5L)
  expect_error(ca_data(d, "first", "nope", "item", "PV1"), "not found")
  d$first[1] <- -1
  expect_error(ca_data(d, "first", "second", "item", "PV1"), "non-negative")
})

test_that("incomplete rows are dropped", {
  d <- simulate_agree(n_schools = 4, n_per_school = 5, seed = 2)
  d$second[1:3] <- NA
  x <- ca_data(d, "first", "second", "item", "theta")
  expect_equal(x$n_dropped, 3L)
})

test_that("Rubin and D1 helpers behave as expected", {
  r <- condagree:::.rubin(matrix(c(1, 1, 1), ncol = 1), matrix(c(4, 4, 4), ncol = 1))
  expect_equal(unname(r$est), 1)
  expect_equal(unname(r$se), 2)
  d1 <- condagree:::.d1(list(c(1, 0)), list(diag(2)))
  expect_equal(unname(d1["p"]), stats::pchisq(1, 2, lower.tail = FALSE))
  d2 <- condagree:::.d1(list(c(1, 0)), list(diag(2)), ncom = 10)
  expect_equal(unname(d2["p"]), stats::pf(0.5, 2, 10, lower.tail = FALSE))
})
