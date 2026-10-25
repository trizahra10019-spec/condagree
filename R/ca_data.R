#' Prepare double-scored responses for ability-conditioned agreement analysis
#'
#' Validates a long-format data set in which each row is one constructed
#' response scored independently by two raters, and stores it in the
#' structure used by all other functions of the package.
#'
#' @param data A data frame in long format (one row per double-scored
#'   response).
#' @param first,second Names of the columns holding the integer scores
#'   (0, 1, ..., m) awarded by the first (operational) and the second
#'   (reliability) rater.
#' @param item Name of the item identifier column.
#' @param ability Name(s) of the column(s) holding the conditioning
#'   variable. Supplying several columns (e.g., five plausible values)
#'   triggers pooling with Rubin's rules in the analysis functions.
#' @param cluster Optional name of the cluster column (e.g., school) used
#'   for the cluster bootstrap and cluster-robust standard errors. If
#'   `NULL`, `id` is used; if both are `NULL`, each row is its own cluster.
#' @param id Optional name of the examinee identifier column.
#' @param max_score Optional named numeric vector giving the maximum score of
#'   each item. By default the maximum observed score (over both raters) is
#'   used.
#'
#' @return An object of class `ca_data`: a list with elements `first`,
#'   `second`, `item` (factor), `ability` (numeric matrix, one column per
#'   plausible value), `cluster`, `id`, `max_score`, and bookkeeping
#'   information.
#'
#' @examples
#' d <- simulate_agree(n_schools = 20, n_per_school = 10, seed = 1)
#' x <- ca_data(d, first = "first", second = "second", item = "item",
#'              ability = paste0("PV", 1:5), cluster = "school", id = "student")
#' x
#' @export
ca_data <- function(data, first, second, item, ability, cluster = NULL,
                    id = NULL, max_score = NULL) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.", call. = FALSE)
  need <- c(first, second, item, ability, cluster, id)
  miss <- setdiff(need, names(data))
  if (length(miss)) {
    stop("Column(s) not found in 'data': ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  keep <- stats::complete.cases(data[, need, drop = FALSE])
  n_drop <- sum(!keep)
  data <- data[keep, , drop = FALSE]
  if (!nrow(data)) stop("No complete rows remain.", call. = FALSE)

  s1 <- as.numeric(data[[first]])
  s2 <- as.numeric(data[[second]])
  if (any(s1 < 0 | s2 < 0) || any(s1 != round(s1)) || any(s2 != round(s2))) {
    stop("Scores must be non-negative integers.", call. = FALSE)
  }
  it <- factor(data[[item]])
  ab <- as.matrix(data[, ability, drop = FALSE])
  storage.mode(ab) <- "double"
  colnames(ab) <- ability

  obs_max <- tapply(pmax(s1, s2), it, max)
  if (is.null(max_score)) {
    ms <- obs_max
  } else {
    if (is.null(names(max_score)) || !all(levels(it) %in% names(max_score))) {
      stop("'max_score' must be named and cover every item.", call. = FALSE)
    }
    ms <- max_score[levels(it)]
    if (any(obs_max > ms)) {
      stop("Observed scores exceed 'max_score' for some items.", call. = FALSE)
    }
  }
  ms <- stats::setNames(as.numeric(ms), levels(it))

  idv <- if (is.null(id)) seq_len(nrow(data)) else data[[id]]
  cl <- if (!is.null(cluster)) data[[cluster]] else idv

  out <- list(first = s1, second = s2, item = it, ability = ab,
              cluster = factor(cl), id = idv, max_score = ms,
              n = length(s1), M = ncol(ab), n_dropped = n_drop,
              vars = list(first = first, second = second, item = item,
                          ability = ability, cluster = cluster, id = id))
  class(out) <- "ca_data"
  out
}

#' @export
print.ca_data <- function(x, ...) {
  ag <- mean(x$first == x$second)
  cat("Double-scored constructed responses (condagree)\n")
  cat(sprintf("  Responses: %d   Items: %d   Clusters: %d   Examinees: %d\n",
              x$n, nlevels(x$item), nlevels(x$cluster), length(unique(x$id))))
  cat(sprintf("  Conditioning variable(s): %s\n",
              paste(colnames(x$ability), collapse = ", ")))
  cat(sprintf("  Global exact agreement: %.4f\n", ag))
  if (x$n_dropped) cat(sprintf("  Rows dropped (incomplete): %d\n", x$n_dropped))
  invisible(x)
}

## ---- internal helpers -------------------------------------------------

# Student-level (or row-level) pooled mean and SD of the conditioning
# variable, used to standardize ability across plausible values.
.ab_scale <- function(x) {
  u <- !duplicated(x$id)
  v <- as.vector(x$ability[u, , drop = FALSE])
  c(mean = mean(v), sd = stats::sd(v))
}

# Rubin's rules for a vector of scalar estimates.
.rubin <- function(est, var, ncom = Inf) {
  est <- as.matrix(est)
  var <- as.matrix(var)
  M <- nrow(est)
  qbar <- colMeans(est)
  ubar <- colMeans(var)
  if (M == 1L) {
    return(list(est = qbar, se = sqrt(ubar), df = rep(ncom, length(qbar)),
                fmi = rep(0, length(qbar))))
  }
  b <- apply(est, 2, stats::var)
  tt <- ubar + (1 + 1 / M) * b
  r <- (1 + 1 / M) * b / pmax(ubar, .Machine$double.eps)
  df <- (M - 1) * (1 + 1 / pmax(r, .Machine$double.eps))^2
  if (is.finite(ncom)) {
    # Barnard and Rubin (1999) small-sample degrees of freedom
    lam <- (1 + 1 / M) * b / pmax(tt, .Machine$double.eps)
    dobs <- (ncom + 1) / (ncom + 3) * ncom * (1 - lam)
    df <- 1 / (1 / df + 1 / dobs)
  }
  list(est = qbar, se = sqrt(tt), df = df, fmi = (r + 2 / (df + 3)) / (r + 1))
}

# D1 statistic (Li, Raghunathan & Rubin, 1991) for a k-dimensional
# parameter estimated in each of M completed data sets.
.d1 <- function(Q, U, ncom = Inf) {
  M <- length(Q)
  k <- length(Q[[1]])
  qbar <- Reduce(`+`, Q) / M
  ubar <- Reduce(`+`, U) / M
  ui <- solve(ubar)
  if (M == 1L) {
    w <- drop(t(qbar) %*% ui %*% qbar)
    return(c(stat = w / k, df1 = k, df2 = ncom,
             p = stats::pf(w / k, k, ncom, lower.tail = FALSE)))
  }
  B <- Reduce(`+`, lapply(Q, function(q) tcrossprod(q - qbar))) / (M - 1)
  r <- (1 + 1 / M) * sum(diag(B %*% ui)) / k
  D <- drop(t(qbar) %*% ui %*% qbar) / (k * (1 + r))
  t_ <- k * (M - 1)
  nu <- if (t_ > 4) {
    4 + (t_ - 4) * (1 + (1 - 2 / t_) / r)^2
  } else {
    t_ * (1 + 1 / k) * (1 + 1 / r)^2 / 2
  }
  nu <- min(nu, ncom)
  c(stat = D, df1 = k, df2 = nu, p = stats::pf(D, k, nu, lower.tail = FALSE))
}
