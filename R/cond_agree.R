#' Ability-conditioned agreement function
#'
#' Estimates exact agreement and item-stratified Cohen's kappa between two
#' raters as a function of the conditioning variable (ability), either with
#' a Gaussian kernel evaluated on a grid or within quantile bands.
#'
#' At ability value \eqn{\theta}, each response receives kernel weight
#' \eqn{w_r(\theta) = \phi((a_r - \theta)/h)}. Conditional exact agreement is
#' \eqn{p_o(\theta) = \sum_r w_r 1(x_{1r} = x_{2r}) / \sum_r w_r}. Chance
#' agreement is computed within items from the weighted marginal score
#' distributions of each rater and pooled over items in proportion to their
#' weight, \eqn{p_e(\theta) = \sum_i W_i(\theta) \sum_k p_{1ik}(\theta)
#' p_{2ik}(\theta) / \sum_i W_i(\theta)}, so that the kappa is not inflated
#' by differences in score distributions across items. With
#' `weights = "quadratic"`, disagreements are weighted by
#' \eqn{(x_1 - x_2)^2 / m_i^2}, where \eqn{m_i} is the item maximum.
#'
#' Uncertainty is quantified with a cluster bootstrap (resampling clusters,
#' typically schools, with replacement). When several plausible values are
#' supplied, the curve is computed for each, the bootstrap variance serves as
#' the within-imputation variance, and the results are pooled with Rubin's
#' rules.
#'
#' @param x A `ca_data` object.
#' @param at Ability values at which to evaluate the kernel estimator. By
#'   default, the 5th, 10th, ..., 95th percentiles of the pooled
#'   student-level conditioning variable.
#' @param bands If an integer (e.g., 5), agreement is computed within that
#'   many quantile bands instead of with a kernel; `at` and `bandwidth` are
#'   then ignored.
#' @param bandwidth Kernel bandwidth on the scale of the conditioning
#'   variable. Default: Silverman's rule of thumb, \eqn{1.06\,s\,n^{-1/5}},
#'   with \eqn{s} and \eqn{n} the pooled SD and number of examinees.
#' @param contrast Two probabilities defining the low and high ability
#'   points whose difference summarizes the gradient (default 10th vs. 90th
#'   percentile). Ignored when `bands` is used, in which case the
#'   contrast is highest minus lowest band.
#' @param contrast_at Optional pair of ability values (on the scale of the
#'   conditioning variable) overriding `contrast`.
#' @param weights `"none"` (default) for exact agreement and Cohen's kappa,
#'   or `"quadratic"` for quadratic-weighted agreement and kappa.
#' @param B Number of bootstrap replications (0 for none).
#' @param level Confidence level.
#' @param seed Optional random seed.
#'
#' @return An object of class `cond_agree` with components `curve` (a data
#'   frame with estimates, standard errors and confidence limits of exact
#'   agreement and kappa at each point), `contrast` (high-minus-low
#'   differences), and settings.
#'
#' @references
#' Cohen, J. (1960). A coefficient of agreement for nominal scales.
#' \emph{Educational and Psychological Measurement, 20}(1), 37--46.
#' \doi{10.1177/001316446002000104}
#'
#' Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse in
#' surveys}. Wiley. \doi{10.1002/9780470316696}
#'
#' @examples
#' d <- simulate_agree(n_schools = 30, n_per_school = 10, gamma = 0.4, seed = 2)
#' x <- ca_data(d, "first", "second", "item", paste0("PV", 1:5),
#'              cluster = "school", id = "student")
#' ca <- cond_agree(x, bands = 4, B = 50, seed = 1)
#' ca
#' @export
cond_agree <- function(x, at = NULL, bands = NULL, bandwidth = NULL,
                       contrast = c(0.10, 0.90), contrast_at = NULL,
                       weights = c("none", "quadratic"), B = 200,
                       level = 0.95, seed = NULL) {
  if (!inherits(x, "ca_data")) stop("'x' must be a 'ca_data' object.", call. = FALSE)
  weights <- match.arg(weights)
  if (!is.null(seed)) set.seed(seed)
  u <- !duplicated(x$id)
  pooled <- as.vector(x$ability[u, , drop = FALSE])
  n_ex <- sum(u)

  agree <- as.numeric(x$first == x$second)
  qw <- if (weights == "quadratic") {
    1 - (x$first - x$second)^2 / pmax(x$max_score[as.character(x$item)], 1)^2
  } else {
    agree
  }
  cats <- .cat_index(x)

  if (!is.null(bands)) {
    bands <- as.integer(bands)
    if (bands < 2) stop("'bands' must be at least 2.", call. = FALSE)
    brk <- stats::quantile(pooled, seq(0, 1, length.out = bands + 1))
    brk[1] <- -Inf
    brk[length(brk)] <- Inf
    pts <- seq_len(bands)
    lab <- paste0("B", pts)
    mkW <- function(a) {
      g <- cut(a, brk, labels = FALSE)
      W <- matrix(0, length(a), bands)
      W[cbind(seq_along(a), g)] <- 1
      W
    }
    lo_i <- 1L
    hi_i <- bands
    at_out <- vapply(pts, function(j) {
      stats::median(pooled[pooled > brk[j] & pooled <= brk[j + 1]])
    }, numeric(1))
    h <- NA_real_
  } else {
    h <- if (is.null(bandwidth)) 1.06 * stats::sd(pooled) * n_ex^(-1 / 5) else bandwidth
    if (is.null(at)) at <- stats::quantile(pooled, seq(0.05, 0.95, by = 0.05), names = FALSE)
    cpts <- if (is.null(contrast_at)) {
      stats::quantile(pooled, contrast, names = FALSE)
    } else {
      as.numeric(contrast_at)[1:2]
    }
    at_out <- c(at, cpts)
    lab <- c(rep("grid", length(at)), "contrast_low", "contrast_high")
    lo_i <- length(at) + 1L
    hi_i <- length(at) + 2L
    mkW <- function(a) {
      exp(-0.5 * ((a - rep(at_out, each = length(a))) / h)^2)
    }
  }

  G <- length(at_out)
  stat_fun <- function(W) {
    s <- .agree_stats(W, agree, qw, cats, weighted = weights == "quadratic")
    c(s$po, s$kappa, s$po[hi_i] - s$po[lo_i], s$kappa[hi_i] - s$kappa[lo_i])
  }
  M <- x$M
  cl <- as.integer(x$cluster)
  nc <- nlevels(x$cluster)
  # same bootstrap draws for every plausible value
  fw <- if (B > 0) {
    replicate(B, tabulate(sample.int(nc, nc, replace = TRUE), nc)[cl])
  } else {
    NULL
  }
  est <- matrix(NA_real_, M, 2 * G + 2)
  vv <- matrix(NA_real_, M, 2 * G + 2)
  for (m in seq_len(M)) {
    W <- matrix(mkW(x$ability[, m]), nrow = x$n)
    est[m, ] <- stat_fun(W)
    if (B > 0) {
      bs <- vapply(seq_len(B), function(b) stat_fun(W * fw[, b]),
                   numeric(2 * G + 2))
      vv[m, ] <- apply(bs, 1, stats::var, na.rm = TRUE)
    } else {
      vv[m, ] <- 0
    }
  }
  pool <- .rubin(est, vv)
  z <- ifelse(is.finite(pool$df),
              stats::qt(1 - (1 - level) / 2, pmax(pool$df, 1)),
              stats::qnorm(1 - (1 - level) / 2))
  lo <- pool$est - z * pool$se
  hi <- pool$est + z * pool$se
  ix_po <- seq_len(G)
  ix_k <- G + seq_len(G)
  keep <- if (is.null(bands)) which(lab == "grid") else seq_len(G)
  curve <- data.frame(
    point = if (is.null(bands)) seq_along(keep) else lab,
    ability = at_out[keep],
    agree = pool$est[ix_po][keep], agree_se = pool$se[ix_po][keep],
    agree_lo = lo[ix_po][keep], agree_hi = hi[ix_po][keep],
    kappa = pool$est[ix_k][keep], kappa_se = pool$se[ix_k][keep],
    kappa_lo = lo[ix_k][keep], kappa_hi = hi[ix_k][keep],
    stringsAsFactors = FALSE)
  if (!is.null(bands)) {
    nb <- table(factor(cut(x$ability[, 1], brk, labels = FALSE), levels = pts))
    curve$n_responses_pv1 <- as.integer(nb)
  }
  ic <- 2 * G + 1:2
  cst <- data.frame(statistic = c("agree", "kappa"),
                    low = at_out[lo_i], high = at_out[hi_i],
                    difference = pool$est[ic], se = pool$se[ic],
                    lower = lo[ic], upper = hi[ic],
                    p = 2 * stats::pt(-abs(pool$est[ic] / pool$se[ic]),
                                      pmax(pool$df[ic], 1)))
  if (B == 0) cst$p <- NA_real_
  out <- list(curve = curve, contrast = cst, bandwidth = h,
              weights = weights, B = B, M = M, level = level,
              method = if (is.null(bands)) "kernel" else "bands",
              global = c(agree = mean(agree),
                         kappa = .agree_stats(matrix(1, x$n, 1), agree, qw, cats,
                                              weighted = weights == "quadratic")$kappa))
  class(out) <- "cond_agree"
  out
}

#' @export
print.cond_agree <- function(x, digits = 3, ...) {
  cat(sprintf("Ability-conditioned agreement (%s; %s weights; B = %d; M = %d)\n",
              x$method, x$weights, x$B, x$M))
  cat(sprintf("Global: agreement = %.4f, kappa = %.4f\n",
              x$global["agree"], x$global["kappa"]))
  print(format(x$curve, digits = digits), row.names = FALSE)
  cat("\nHigh-minus-low contrast:\n")
  print(format(x$contrast, digits = digits), row.names = FALSE)
  invisible(x)
}

#' Plot an ability-conditioned agreement function
#'
#' @param x A `cond_agree` object.
#' @param which `"kappa"` (default) or `"agree"`.
#' @param xlab,ylab Axis labels.
#' @param ... Further arguments passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @examples
#' d <- simulate_agree(n_schools = 30, n_per_school = 10, gamma = 0.4, seed = 2)
#' x <- ca_data(d, "first", "second", "item", "theta", cluster = "school",
#'              id = "student")
#' plot(cond_agree(x, B = 30, seed = 1))
#' @export
plot.cond_agree <- function(x, which = c("kappa", "agree"), xlab = "Ability",
                            ylab = NULL, ...) {
  which <- match.arg(which)
  cv <- x$curve
  y <- cv[[which]]
  lo <- cv[[paste0(which, "_lo")]]
  hi <- cv[[paste0(which, "_hi")]]
  if (is.null(ylab)) ylab <- if (which == "kappa") "Conditional kappa" else "Conditional exact agreement"
  graphics::plot(cv$ability, y, type = "n", ylim = range(c(lo, hi, y), na.rm = TRUE),
                 xlab = xlab, ylab = ylab, las = 1, ...)
  if (x$method == "kernel") {
    graphics::polygon(c(cv$ability, rev(cv$ability)), c(lo, rev(hi)),
                      col = grDevices::grey(0.85), border = NA)
    graphics::lines(cv$ability, y, lwd = 2)
  } else {
    graphics::arrows(cv$ability, lo, cv$ability, hi, angle = 90, code = 3,
                     length = 0.04)
    graphics::points(cv$ability, y, pch = 19)
  }
  graphics::abline(h = x$global[which], lty = 2)
  invisible(x)
}

## ---- internal --------------------------------------------------------

# Index structure for item-by-category indicators of each rater.
.cat_index <- function(x) {
  itc <- as.character(x$item)
  lev1 <- paste(itc, x$first, sep = "\r")
  lev2 <- paste(itc, x$second, sep = "\r")
  levs <- sort(unique(c(lev1, lev2)))
  f1 <- match(lev1, levs)
  f2 <- match(lev2, levs)
  cat_item <- sub("\r.*$", "", levs)
  cat_score <- as.numeric(sub("^.*\r", "", levs))
  list(f1 = f1, f2 = f2, ncat = length(levs),
       cat_item = factor(cat_item, levels = levels(x$item)),
       cat_score = cat_score,
       max_score = x$max_score)
}

# Weighted agreement statistics for each column of the weight matrix W.
.agree_stats <- function(W, agree, qw, cats, weighted = FALSE) {
  sw <- colSums(W)
  po <- drop(crossprod(if (weighted) qw else agree, W)) / sw
  P1 <- rowsum(W, cats$f1, reorder = TRUE)
  P2 <- rowsum(W, cats$f2, reorder = TRUE)
  full1 <- matrix(0, cats$ncat, ncol(W))
  full2 <- full1
  full1[sort(unique(cats$f1)), ] <- P1
  full2[sort(unique(cats$f2)), ] <- P2
  Wi <- rowsum(full1, cats$cat_item, reorder = TRUE)
  if (!weighted) {
    num <- rowsum(full1 * full2, cats$cat_item, reorder = TRUE)
  } else {
    num <- matrix(0, nrow(Wi), ncol(W))
    for (i in seq_len(nrow(Wi))) {
      idx <- which(as.integer(cats$cat_item) == i)
      if (!length(idx)) next
      sc <- cats$cat_score[idx]
      mi <- max(cats$max_score[i], 1)
      Om <- 1 - outer(sc, sc, function(a, b) (a - b)^2) / mi^2
      num[i, ] <- colSums(full1[idx, , drop = FALSE] *
                            (Om %*% full2[idx, , drop = FALSE]))
    }
  }
  pe <- colSums(num / ifelse(Wi > 0, Wi, 1)) / colSums(Wi)
  kappa <- (po - pe) / (1 - pe)
  list(po = po, pe = pe, kappa = kappa)
}
