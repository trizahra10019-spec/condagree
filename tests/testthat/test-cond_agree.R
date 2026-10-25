test_that("kappa equals Cohen's kappa for a single item without conditioning", {
  set.seed(3)
  n <- 400
  s1 <- sample(0:2, n, replace = TRUE, prob = c(.5, .3, .2))
  s2 <- ifelse(runif(n) < .85, s1, sample(0:2, n, replace = TRUE))
  d <- data.frame(item = "A", first = s1, second = s2, ab = rnorm(n), id = 1:n)
  x <- ca_data(d, "first", "second", "item", "ab", id = "id")
  tab <- table(factor(s1, 0:2), factor(s2, 0:2)) / n
  po <- sum(diag(tab)); pe <- sum(rowSums(tab) * colSums(tab))
  ca <- cond_agree(x, B = 0)
  expect_equal(unname(ca$global["kappa"]), (po - pe) / (1 - pe), tolerance = 1e-12)
  expect_equal(unname(ca$global["agree"]), po, tolerance = 1e-12)
  # with very large bandwidth, conditional values equal global values
  cb <- cond_agree(x, bandwidth = 1e6, B = 0)
  expect_equal(cb$curve$kappa, rep(unname(ca$global["kappa"]), nrow(cb$curve)),
               tolerance = 1e-8)
})

test_that("quadratic weights reproduce weighted kappa", {
  set.seed(4)
  n <- 300
  s1 <- sample(0:3, n, replace = TRUE)
  s2 <- pmin(3, pmax(0, s1 + sample(-1:1, n, replace = TRUE, prob = c(.1, .8, .1))))
  d <- data.frame(item = "A", first = s1, second = s2, ab = rnorm(n), id = 1:n)
  x <- ca_data(d, "first", "second", "item", "ab", id = "id")
  t <- table(factor(s1, 0:3), factor(s2, 0:3)) / n
  W <- 1 - outer(0:3, 0:3, function(a, b) (a - b)^2) / 9
  po <- sum(W * t); pe <- sum(W * outer(rowSums(t), colSums(t)))
  ca <- cond_agree(x, B = 0, weights = "quadratic")
  expect_equal(unname(ca$global["kappa"]), (po - pe) / (1 - pe), tolerance = 1e-12)
})

test_that("bands, bootstrap and plausible-value pooling run", {
  d <- simulate_agree(n_schools = 12, n_per_school = 10, gamma = 0.5, seed = 5)
  x <- ca_data(d, "first", "second", "item", paste0("PV", 1:3),
               cluster = "school", id = "student")
  ca <- cond_agree(x, bands = 3, B = 20, seed = 1)
  expect_s3_class(ca, "cond_agree")
  expect_equal(nrow(ca$curve), 3L)
  expect_true(all(ca$curve$kappa_lo <= ca$curve$kappa_hi))
  expect_true(all(is.finite(ca$contrast$se)))
  expect_invisible(plot(ca))
})
