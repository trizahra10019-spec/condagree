#' Test whether rater agreement is uniform across ability
#'
#' Fits a fixed-effects logistic regression of disagreement
#' (\eqn{x_1 \neq x_2}) on the standardized conditioning variable, with
#' stratum intercepts that absorb item differences, and tests the ability
#' terms with cluster-robust (sandwich) covariance matrices.
#'
#' With `adjust = "none"` (the default) the strata are items, so the
#' ability coefficient captures the total ability gradient in
#' disagreement, and the test is a valid test of the null hypothesis that
#' agreement is uniform across ability. With `adjust = "score"` the strata
#' are item-by-score categories of the first (operational) rater. This
#' naive adjustment is provided for comparison only: because the first
#' rater's score is itself error-prone, conditioning on it does not remove
#' the score-opportunity effect, and the resulting test is anti-conservative
#' (see the simulation study in the package documentation). Use
#' [rater_model()] to separate score opportunity from a residual ability
#' effect.
#'
#' For each specification, a 1-df trend test (linear term) and an omnibus
#' test of uniformity (natural cubic spline with `df` degrees of freedom)
#' are reported. Strata in which all responses agree (or all disagree)
#' carry no information about the ability terms and are dropped. With
#' several plausible values, trend estimates are pooled with Rubin's rules
#' and omnibus tests with the D1 statistic of Li, Raghunathan and Rubin
#' (1991). Because cluster-robust covariance matrices are imprecise with a
#' limited number of clusters \eqn{G}, the trend test uses a \eqn{t}
#' reference distribution with \eqn{G - 1} complete-data degrees of freedom
#' (with the Barnard-Rubin adjustment for plausible values) and the omnibus
#' test an \eqn{F} reference distribution whose denominator degrees of
#' freedom are at most \eqn{G - 1}. Spline knots are fixed at quantiles of
#' the pooled conditioning
#' variable so that the basis is identical across plausible values.
#'
#' @param x A `ca_data` object.
#' @param adjust Character vector with one or both of `"none"` and
#'   `"score"`.
#' @param df Degrees of freedom of the spline for the omnibus test (1 gives
#'   a linear omnibus test identical to the trend test).
#' @param unit Optional scaling constant for reporting the trend on the
#'   original metric (e.g., 100 for "per 100 scale points"). The slope per
#'   pooled standard deviation is always reported.
#'
#' @return An object of class `agree_test` with a data frame `results`.
#'
#' @references
#' Li, K. H., Raghunathan, T. E., & Rubin, D. B. (1991). Large-sample
#' significance levels from multiply imputed data using moment-based
#' statistics and an F reference distribution. \emph{Journal of the
#' American Statistical Association, 86}(416), 1065--1073.
#' \doi{10.1080/01621459.1991.10475152}
#'
#' Liang, K.-Y., & Zeger, S. L. (1986). Longitudinal data analysis using
#' generalized linear models. \emph{Biometrika, 73}(1), 13--22.
#' \doi{10.1093/biomet/73.1.13}
#'
#' @examples
#' d <- simulate_agree(n_schools = 30, n_per_school = 10, gamma = 0.4, seed = 3)
#' x <- ca_data(d, "first", "second", "item", paste0("PV", 1:5),
#'              cluster = "school", id = "student")
#' agree_test(x)
#' @export
agree_test <- function(x, adjust = "none", df = 3, unit = NULL) {
  if (!inherits(x, "ca_data")) stop("'x' must be a 'ca_data' object.", call. = FALSE)
  adjust <- match.arg(adjust, c("none", "score"), several.ok = TRUE)
  sc <- .ab_scale(x)
  y <- as.numeric(x$first != x$second)
  u <- !duplicated(x$id)
  pooled_z <- (as.vector(x$ability[u, , drop = FALSE]) - sc["mean"]) / sc["sd"]
  if (df > 1) {
    kn <- stats::quantile(pooled_z, seq(0, 1, length.out = df + 1)[-c(1, df + 1)],
                          names = FALSE)
    bk <- range(pooled_z)
  }
  basis <- function(z) {
    if (df == 1) return(matrix(z, ncol = 1))
    unclass(splines::ns(z, knots = kn, Boundary.knots = bk))[, , drop = FALSE]
  }
  res <- list()
  slopes <- c()
  for (adj in adjust) {
    strat <- if (adj == "none") x$item else interaction(x$item, x$first, drop = TRUE)
    ev <- tapply(y, strat, sum)
    nn <- tapply(y, strat, length)
    inf <- names(ev)[ev > 0 & ev < nn]
    keep <- as.character(strat) %in% inf
    if (sum(keep) == 0 || length(unique(y[keep])) < 2) {
      stop("No informative strata (no disagreements).", call. = FALSE)
    }
    st <- droplevels(factor(strat[keep]))
    yk <- y[keep]
    clk <- droplevels(x$cluster[keep])
    D <- stats::model.matrix(~ st - 1)
    trQ <- list(); trU <- list(); omQ <- list(); omU <- list()
    for (m in seq_len(x$M)) {
      z <- (x$ability[keep, m] - sc["mean"]) / sc["sd"]
      f1 <- .fe_logit(D, matrix(z, ncol = 1), yk, clk)
      trQ[[m]] <- f1$b
      trU[[m]] <- f1$V
      if (df > 1) {
        fo <- .fe_logit(D, basis(z), yk, clk)
        omQ[[m]] <- fo$b
        omU[[m]] <- fo$V
      }
    }
    ncom <- nlevels(clk) - 1
    tr <- .rubin(matrix(unlist(trQ), ncol = 1),
                 matrix(vapply(trU, function(v) v[1, 1], numeric(1)), ncol = 1),
                 ncom = ncom)
    tstat <- tr$est / tr$se
    ptr <- 2 * stats::pt(-abs(tstat), pmax(tr$df, 1))
    om <- if (df > 1) .d1(omQ, omU, ncom = ncom) else c(stat = tstat^2, df1 = 1, df2 = tr$df, p = ptr)
    slopes[adj] <- tr$est
    res[[adj]] <- data.frame(
      adjust = adj, strata = nlevels(st), responses = length(yk),
      disagreements = sum(yk),
      slope_sd = tr$est, se_sd = tr$se, or_sd = exp(tr$est),
      slope_unit = if (is.null(unit)) NA_real_ else tr$est * unit / sc["sd"],
      or_unit = if (is.null(unit)) NA_real_ else exp(tr$est * unit / sc["sd"]),
      trend_t = tstat, trend_df = tr$df, trend_p = ptr,
      omnibus_F = om["stat"], omnibus_df1 = om["df1"], omnibus_df2 = om["df2"],
      omnibus_p = om["p"], row.names = NULL)
  }
  out <- list(results = do.call(rbind, res), df = df, unit = unit, M = x$M,
              scale = sc,
              share_score = if (all(c("none", "score") %in% adjust)) {
                unname(1 - slopes["score"] / slopes["none"])
              } else NA_real_)
  class(out) <- "agree_test"
  out
}

#' @export
print.agree_test <- function(x, digits = 3, ...) {
  cat(sprintf("Uniformity of rater agreement across ability (M = %d)\n", x$M))
  cat("Trend: 1-df logit slope per pooled SD (cluster-robust);",
      sprintf("omnibus: spline with %d df.\n", x$df))
  r <- x$results
  show <- r[, c("adjust", "strata", "responses", "disagreements", "slope_sd",
                "se_sd", "or_sd", "trend_p", "omnibus_F", "omnibus_p")]
  print(format(show, digits = digits), row.names = FALSE)
  invisible(x)
}

#' Test for directional rater severity across ability
#'
#' Among discordant responses, models the probability that the second rater
#' awarded the lower score as a logistic function of standardized ability.
#' The intercept tests symmetric disagreement at average ability
#' (\eqn{P = .5}); the slope tests whether the direction of disagreement
#' changes with ability; a joint 2-df test (D1 when plausible values are
#' supplied) tests for any systematic severity. Reference distributions use
#' \eqn{G - 1} complete-data degrees of freedom, as in [agree_test()].
#'
#' @param x A `ca_data` object.
#' @return An object of class `severity_test`.
#' @examples
#' d <- simulate_agree(n_schools = 30, n_per_school = 10, severity = 0.8, seed = 4)
#' x <- ca_data(d, "first", "second", "item", "theta", cluster = "school",
#'              id = "student")
#' severity_test(x)
#' @export
severity_test <- function(x) {
  if (!inherits(x, "ca_data")) stop("'x' must be a 'ca_data' object.", call. = FALSE)
  sc <- .ab_scale(x)
  dis <- x$first != x$second
  if (sum(dis) < 10) stop("Fewer than 10 discordant responses.", call. = FALSE)
  y <- as.numeric(x$second[dis] < x$first[dis])
  clk <- droplevels(x$cluster[dis])
  Q <- list(); U <- list()
  for (m in seq_len(x$M)) {
    z <- (x$ability[dis, m] - sc["mean"]) / sc["sd"]
    f <- .fe_logit(matrix(1, length(y), 1), matrix(z, ncol = 1), y, clk,
                   return_all = TRUE)
    Q[[m]] <- f$b_all
    U[[m]] <- f$V_all
  }
  ncom <- nlevels(clk) - 1
  r0 <- .rubin(do.call(rbind, lapply(Q, function(q) q)),
               do.call(rbind, lapply(U, function(v) diag(v))), ncom = ncom)
  tt <- r0$est / r0$se
  p <- 2 * stats::pt(-abs(tt), pmax(r0$df, 1))
  joint <- .d1(Q, U, ncom = ncom)
  tab <- data.frame(term = c("intercept", "ability_sd"), estimate = r0$est,
                    se = r0$se, t = tt, df = r0$df, p = p, row.names = NULL)
  out <- list(coefficients = tab, joint = joint, n_discordant = length(y),
              prop_second_lower = mean(y), M = x$M)
  class(out) <- "severity_test"
  out
}

#' @export
print.severity_test <- function(x, digits = 3, ...) {
  cat(sprintf("Directional severity among %d discordant responses (M = %d)\n",
              x$n_discordant, x$M))
  cat(sprintf("Proportion in which the second rater scored lower: %.3f\n",
              x$prop_second_lower))
  print(format(x$coefficients, digits = digits), row.names = FALSE)
  cat(sprintf("Joint test of no severity: F(%g, %.1f) = %.3f, p = %.4f\n",
              x$joint["df1"], x$joint["df2"], x$joint["stat"], x$joint["p"]))
  invisible(x)
}

## ---- internal --------------------------------------------------------

# Logistic regression with stratum dummies D and covariates Z; returns the
# Z coefficients and their cluster-robust (CR1) covariance.
.fe_logit <- function(D, Z, y, cluster, return_all = FALSE) {
  X <- cbind(D, Z)
  fit <- suppressWarnings(stats::glm.fit(X, y, family = stats::binomial(),
                                         control = list(maxit = 50)))
  b <- fit$coefficients
  ok <- !is.na(b)
  X <- X[, ok, drop = FALSE]
  b <- b[ok]
  mu <- fit$fitted.values
  w <- mu * (1 - mu)
  bread <- solve(crossprod(X * sqrt(w)))
  sc <- rowsum(X * (y - mu), cluster)
  G <- nrow(sc)
  meat <- crossprod(sc) * G / (G - 1)
  V <- bread %*% meat %*% bread
  zi <- (ncol(X) - ncol(Z) + 1):ncol(X)
  out <- list(b = b[zi], V = V[zi, zi, drop = FALSE])
  if (return_all) {
    out$b_all <- b
    out$V_all <- V
  }
  out
}
