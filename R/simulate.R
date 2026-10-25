#' Simulate double-scored constructed responses
#'
#' Generates two-rater scores for polytomous constructed-response items in
#' a clustered (students-within-schools) design. True scores follow a
#' generalized partial credit model with unit discrimination. Each rater
#' independently misassigns a response to an adjacent category with
#' probability \eqn{\mathrm{logit}^{-1}(\alpha_c + \gamma \theta)}, where
#' \eqn{\alpha_c} depends on whether the true category is the lowest
#' (`zero`), an intermediate (`partial`), or the highest (`full`) score.
#' Unequal \eqn{\alpha_c} create a score-opportunity effect; \eqn{\gamma}
#' creates a residual ability effect. Optionally the second rater
#' errs downward more often (`severity`, on the logit scale). Plausible
#' values are drawn from the normal posterior of ability given an error-prone
#' observed score with the stated reliability.
#'
#' @param n_schools,n_per_school Number of clusters and examinees per
#'   cluster.
#' @param max_scores Integer vector of item maxima (one element per item).
#' @param icc Intraclass correlation of ability.
#' @param gamma Residual ability effect on each rater's error logit (per SD).
#' @param alpha Named vector with error-logit intercepts `zero`, `partial`,
#'   and `full`.
#' @param alpha_item Optional numeric vector (one element per item) added to
#'   the error intercepts of every category of that item, to generate
#'   item-specific scoring difficulty.
#' @param severity Logit of the probability that an error by the second
#'   rater is downward minus the logit for the first rater (0 = none).
#' @param reliability Reliability of the observed ability measure used to
#'   generate plausible values (1 gives plausible values equal to the true
#'   ability).
#' @param n_pv Number of plausible values.
#' @param seed Optional random seed.
#'
#' @return A long data frame with columns `school`, `student`, `item`,
#'   `max_score`, `true`, `first`, `second`, `theta`, and `PV1`, ...
#' @examples
#' d <- simulate_agree(n_schools = 10, n_per_school = 10, seed = 1)
#' head(d)
#' @export
simulate_agree <- function(n_schools = 50, n_per_school = 20,
                           max_scores = c(1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2),
                           icc = 0.2, gamma = 0,
                           alpha = c(zero = -3.2, partial = -3.2, full = -3.2),
                           alpha_item = NULL, severity = 0, reliability = 0.85,
                           n_pv = 5,
                           seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (!all(c("zero", "partial", "full") %in% names(alpha))) {
    stop("'alpha' must have elements 'zero', 'partial', 'full'.", call. = FALSE)
  }
  N <- n_schools * n_per_school
  J <- length(max_scores)
  school <- rep(seq_len(n_schools), each = n_per_school)
  theta <- sqrt(icc) * stats::rnorm(n_schools)[school] +
    sqrt(1 - icc) * stats::rnorm(N)
  loc <- if (J > 1) seq(-1, 1, length.out = J) else 0
  if (is.null(alpha_item)) alpha_item <- rep(0, J)
  if (length(alpha_item) != J) stop("'alpha_item' must have one element per item.", call. = FALSE)

  rows <- vector("list", J)
  for (j in seq_len(J)) {
    m <- max_scores[j]
    steps <- if (m == 1) loc[j] else loc[j] + seq(-0.6, 0.6, length.out = m)
    # GPCM category probabilities
    lin <- cbind(0, sapply(seq_len(m), function(k) k * theta - sum(steps[seq_len(k)])))
    pr <- exp(lin - apply(lin, 1, max))
    pr <- pr / rowSums(pr)
    cu <- t(apply(pr, 1, cumsum))
    tr <- rowSums(stats::runif(N) > cu)
    pos <- ifelse(tr == 0, "zero", ifelse(tr == m, "full", "partial"))
    rate <- function(down_shift) {
      err <- stats::runif(N) < stats::plogis(alpha[pos] + alpha_item[j] + gamma * theta)
      down <- stats::runif(N) < stats::plogis(down_shift)
      dir <- ifelse(tr == 0, 1, ifelse(tr == m, -1, ifelse(down, -1, 1)))
      ifelse(err, tr + dir, tr)
    }
    rows[[j]] <- data.frame(school = school, student = seq_len(N),
                            item = sprintf("I%02d", j), max_score = m,
                            true = tr, first = rate(0), second = rate(severity),
                            theta = theta)
  }
  d <- do.call(rbind, rows)
  s2e <- (1 - reliability) / reliability
  if (reliability < 1) {
    obs <- theta + stats::rnorm(N, 0, sqrt(s2e))
    pm <- obs / (1 + s2e)
    psd <- sqrt(s2e / (1 + s2e))
  }
  for (k in seq_len(n_pv)) {
    pv <- if (reliability < 1) pm + stats::rnorm(N, 0, psd) else theta
    d[[paste0("PV", k)]] <- pv[d$student]
  }
  d
}
