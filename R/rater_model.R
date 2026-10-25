#' Latent-category rater error model conditioned on ability
#'
#' Separates the \emph{score-opportunity} component of an ability gradient
#' in rater disagreement from a \emph{residual ability} component by
#' treating the correct score of each response as a latent category.
#'
#' For a response to item \eqn{i} (maximum \eqn{m_i}) by an examinee with
#' standardized ability \eqn{z}, the latent correct category \eqn{C} follows
#' a generalized partial credit model with item-specific slope \eqn{a_i}
#' and step parameters \eqn{b_{il}},
#' \eqn{\log\{P(C = l)/P(C = l - 1)\} = a_i z - b_{il}}. Each rater
#' \eqn{r = 1, 2} reproduces the latent category with probability
#' \eqn{1 - e(C, z)} and otherwise assigns an adjacent category, where
#' \deqn{\mathrm{logit}\, e(C, z) = \alpha_{p(C)} + \gamma z,}
#' and \eqn{p(C)} indicates whether \eqn{C} is the lowest (zero), an
#' intermediate (partial), or the highest (full) category. Errors from an
#' intermediate category are downward with probability
#' \eqn{\mathrm{logit}^{-1}(\delta_r)}. Thus \eqn{\alpha_p} captures the
#' score-opportunity effect (some categories are intrinsically harder to
#' score), \eqn{\gamma} is the residual ability effect on scoring error,
#' and \eqn{\delta_2 - \delta_1} is directional severity of the second rater
#' relative to the first. The two raters are exchangeable apart from
#' \eqn{\delta_r}, which identifies the model with two ratings per response.
#'
#' Parameters are estimated by maximum likelihood with the BHHH algorithm,
#' using exact per-response scores obtained from the Fisher identity, with
#' a weak normal penalty (SD 10) on step and error-intercept parameters to
#' keep estimates finite when a category is empty or virtually error-free. Standard errors are cluster-robust sandwich estimates with the
#' observed information as bread. With
#' several plausible values, the model is fitted to each and results are
#' pooled with Rubin's rules. Confidence intervals and tests use a \eqn{t}
#' reference distribution with the complete-data degrees of freedom set to
#' the number of clusters minus one and, with plausible values, the
#' small-sample degrees of freedom of Barnard and Rubin (1999).
#'
#' The function also reports a model-based decomposition: the implied
#' disagreement curve (averaged over items) under the fitted model and under
#' a counterfactual with \eqn{\gamma = 0} (error rates fixed at their values
#' for average ability). The share of the high-minus-low gradient that
#' remains under the counterfactual is the share attributable to score
#' opportunity. Its standard error is obtained by the delta method (numerical
#' derivatives), or optionally by drawing parameters from their estimated
#' sampling distribution. The delta method is preferred because weakly
#' identified nuisance parameters near the boundary (e.g., an almost
#' error-free zero category) make normal parameter draws unstable.
#'
#' After BHHH iterations, the estimates are polished with Newton steps that
#' use the observed information (a numerical Hessian of the analytic
#' penalized gradient). The same matrix is the bread of the cluster-robust
#' sandwich covariance. Convergence is declared when the Newton decrement
#' \eqn{g^\top I^{-1} g} falls below \eqn{10^{-6}}.
#'
#' @param x A `ca_data` object.
#' @param alpha_by Grouping of the error intercepts: `"position"` (default;
#'   one intercept per category position, shared by all items), `"maxscore"`
#'   (separate intercepts for items with different maximum scores), or
#'   `"item"` (item-specific intercepts for each category position). The
#'   finer groupings guard against confounding of the residual ability
#'   effect by item composition, for example when items that are harder to
#'   score are also those on which able examinees earn partial credit.
#' @param contrast Probabilities of the pooled conditioning variable that
#'   define the low and high points of the gradient (default 10th and 90th
#'   percentiles).
#' @param grid Number of grid points for the implied curves.
#' @param se_decomp `"delta"` (default) or `"draws"`: method for the
#'   standard errors of the decomposition.
#' @param ndraws Number of parameter draws when `se_decomp = "draws"`.
#' @param maxit,tol Maximum BHHH iterations and convergence tolerance.
#' @param seed Optional random seed (affects only `se_decomp = "draws"`).
#'
#' @return An object of class `rater_model` with pooled `coefficients`,
#'   the pooled error intercepts (`alpha`), the number of parameters and the
#'   AIC (based on the mean log-likelihood over plausible values; useful for
#'   comparing `alpha_by` groupings),
#'   tests of \eqn{\gamma = 0} and \eqn{\delta_2 = \delta_1},
#'   `decomposition`, implied `curves`, a `fit_check` comparing observed and
#'   implied disagreement by ability quintile, and per-plausible-value fits.
#'
#' @references
#' Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with
#' multiple imputation. \emph{Biometrika, 86}(4), 948--955.
#' \doi{10.1093/biomet/86.4.948}
#'
#' Patz, R. J., Junker, B. W., Johnson, M. S., & Mariano, L. T. (2002). The
#' hierarchical rater model for rated test items and its application to
#' large-scale educational assessment data. \emph{Journal of Educational and
#' Behavioral Statistics, 27}(4), 341--384. \doi{10.3102/10769986027004341}
#'
#' @examples
#' d <- simulate_agree(n_schools = 30, n_per_school = 10, gamma = 0.4,
#'                     alpha = c(zero = -4.5, partial = -2.2, full = -3.5),
#'                     reliability = 1, n_pv = 1, seed = 5)
#' x <- ca_data(d, "first", "second", "item", "PV1", cluster = "school",
#'              id = "student")
#' rm <- rater_model(x)
#' rm
#' @export
rater_model <- function(x, alpha_by = c("position", "maxscore", "item"),
                        contrast = c(0.10, 0.90), grid = 41,
                        se_decomp = c("delta", "draws"), ndraws = 500,
                        maxit = 200, tol = 1e-8, seed = NULL) {
  if (!inherits(x, "ca_data")) stop("'x' must be a 'ca_data' object.", call. = FALSE)
  se_decomp <- match.arg(se_decomp)
  alpha_by <- match.arg(alpha_by)
  if (!is.null(seed)) set.seed(seed)
  S <- .rm_setup(x, alpha_by)
  sc <- .ab_scale(x)
  u <- !duplicated(x$id)
  pz <- (as.vector(x$ability[u, , drop = FALSE]) - sc["mean"]) / sc["sd"]
  zq <- stats::quantile(pz, contrast, names = FALSE)
  zg <- seq(stats::quantile(pz, 0.02), stats::quantile(pz, 0.98), length.out = grid)
  wj <- as.numeric(table(S$j)) / S$n

  fits <- vector("list", x$M)
  start <- NULL
  est <- NULL
  vv <- NULL
  est_a <- NULL
  vv_a <- NULL
  curves_full <- matrix(NA_real_, x$M, grid)
  curves_opp <- matrix(NA_real_, x$M, grid)
  for (m in seq_len(x$M)) {
    z <- (x$ability[, m] - sc["mean"]) / sc["sd"]
    f <- .rm_fit(S, z, start = start, maxit = maxit, tol = tol)
    start <- f$par
    ip <- S$idx
    g <- f$par[ip$gamma]
    sev <- if (S$has_partial) f$par[ip$delta[2]] - f$par[ip$delta[1]] else NA_real_
    vsev <- if (S$has_partial) {
      f$V[ip$delta[2], ip$delta[2]] + f$V[ip$delta[1], ip$delta[1]] -
        2 * f$V[ip$delta[1], ip$delta[2]]
    } else NA_real_
    # decomposition with parameter draws
    dfun <- function(p) {
      full <- .rm_curve(S, p, c(zq, zg), wj, zero_gamma = FALSE)
      opp <- .rm_curve(S, p, c(zq, zg), wj, zero_gamma = TRUE)
      dfull <- full[2] - full[1]
      dopp <- opp[2] - opp[1]
      c(dfull, dopp, dopp / dfull, full[-(1:2)], opp[-(1:2)])
    }
    pt <- dfun(f$par)
    dfun2 <- function(p) {
      full <- .rm_curve(S, p, zq, wj, zero_gamma = FALSE)
      opp <- .rm_curve(S, p, zq, wj, zero_gamma = TRUE)
      c(full[2] - full[1], opp[2] - opp[1], (opp[2] - opp[1]) / (full[2] - full[1]))
    }
    if (se_decomp == "draws") {
      dr <- .mvdraw(f$par, f$V, ndraws)
      dd <- apply(dr, 1, dfun2)
      dv <- apply(dd, 1, stats::var)
    } else {
      Jd <- vapply(seq_along(f$par), function(k) {
        h <- 1e-4 * max(1, abs(f$par[k]))
        e <- replace(numeric(length(f$par)), k, h)
        (dfun2(f$par + e) - dfun2(f$par - e)) / (2 * h)
      }, numeric(3))
      dv <- diag(Jd %*% f$V %*% t(Jd))
    }
    curves_full[m, ] <- pt[3 + seq_len(grid)]
    curves_opp[m, ] <- pt[3 + grid + seq_len(grid)]
    ia <- if (alpha_by == "position") ip$alpha else integer(0)
    th <- c(f$par[ia], unname(g), f$par[ip$delta], unname(sev), pt[1:3])
    names(th) <- c(S$names[ia], "gamma", S$names[ip$delta], "severity",
                   "grad_total", "grad_opportunity", "share_opportunity")
    tv <- c(diag(f$V)[ia], f$V[ip$gamma, ip$gamma], diag(f$V)[ip$delta],
            vsev, dv)
    est <- rbind(est, th)
    vv <- rbind(vv, tv)
    est_a <- rbind(est_a, f$par[ip$alpha])
    vv_a <- rbind(vv_a, diag(f$V)[ip$alpha])
    fits[[m]] <- f
  }
  pool <- .rubin(est, vv, ncom = nlevels(x$cluster) - 1)
  tt <- pool$est / pool$se
  coefs <- data.frame(parameter = colnames(est), estimate = pool$est, se = pool$se,
                      df = pool$df,
                      lower = pool$est - stats::qt(0.975, pmax(pool$df, 1)) * pool$se,
                      upper = pool$est + stats::qt(0.975, pmax(pool$df, 1)) * pool$se,
                      p = 2 * stats::pt(-abs(tt), pmax(pool$df, 1)),
                      row.names = NULL, stringsAsFactors = FALSE)
  coefs$p[grepl("^alpha", coefs$parameter)] <- NA_real_
  coefs$p[coefs$parameter %in% c("share_opportunity")] <- NA_real_

  # fit check on the first plausible value
  z1 <- (x$ability[, 1] - sc["mean"]) / sc["sd"]
  pdis <- .rm_rowdis(S, fits[[1]]$par, z1)
  qb <- cut(z1, stats::quantile(z1, 0:5 / 5), include.lowest = TRUE, labels = FALSE)
  obs <- as.numeric(x$first != x$second)
  fc <- data.frame(quintile = 1:5,
                   observed = as.numeric(tapply(obs, qb, mean)),
                   implied = as.numeric(tapply(pdis, qb, mean)),
                   n = as.integer(table(qb)))

  pa <- .rubin(est_a, vv_a, ncom = nlevels(x$cluster) - 1)
  alpha_tab <- data.frame(parameter = S$names[S$idx$alpha], estimate = pa$est,
                          se = pa$se, row.names = NULL, stringsAsFactors = FALSE)
  ll_mean <- mean(vapply(fits, `[[`, numeric(1), "loglik"))
  out <- list(coefficients = coefs, alpha = alpha_tab, alpha_by = alpha_by,
              npar = S$P, aic = -2 * ll_mean + 2 * S$P,
              curves = data.frame(z = zg, ability = sc["mean"] + sc["sd"] * zg,
                                  disagree_full = colMeans(curves_full),
                                  disagree_opportunity = colMeans(curves_opp),
                                  row.names = NULL),
              contrast_z = zq, fit_check = fc, fits = fits, M = x$M,
              scale = sc, converged = vapply(fits, `[[`, logical(1), "converged"),
              loglik = vapply(fits, `[[`, numeric(1), "loglik"))
  class(out) <- "rater_model"
  out
}

#' @export
print.rater_model <- function(x, digits = 3, ...) {
  cat(sprintf("Latent-category rater error model (M = %d; alpha_by = %s; converged: %s)\n",
              x$M, x$alpha_by, paste(unique(x$converged), collapse = "/")))
  cat(sprintf("Parameters: %d; AIC: %.1f\n", x$npar, x$aic))
  print(format(x$coefficients, digits = digits), row.names = FALSE)
  cat("\nObserved vs implied disagreement by ability quintile (PV 1):\n")
  print(format(x$fit_check, digits = digits), row.names = FALSE)
  invisible(x)
}

#' @export
plot.rater_model <- function(x, xlab = "Ability", ylab = "Implied disagreement",
                             ...) {
  cv <- x$curves
  graphics::plot(cv$ability, cv$disagree_full, type = "l", lwd = 2, las = 1,
                 ylim = range(c(0, cv$disagree_full, cv$disagree_opportunity)),
                 xlab = xlab, ylab = ylab, ...)
  graphics::lines(cv$ability, cv$disagree_opportunity, lwd = 2, lty = 2)
  graphics::legend("topleft", c("Fitted model", "Score opportunity only"),
                   lty = 1:2, lwd = 2, bty = "n")
  invisible(x)
}

## ---- internal --------------------------------------------------------

.rm_setup <- function(x, alpha_by = "position") {
  ms <- x$max_score
  if (any(ms < 1)) {
    stop("All items must have a maximum score of at least 1.", call. = FALSE)
  }
  j <- as.integer(x$item)
  J <- nlevels(x$item)
  mj <- as.integer(ms)
  K <- max(mj) + 1L
  n <- x$n
  cc <- matrix(0:(K - 1), n, K, byrow = TRUE)
  valid <- cc <= mj[j]
  pos <- ifelse(cc == 0, 1L, ifelse(cc == mj[j], 3L, 2L))
  pos[!valid] <- NA_integer_
  has_partial <- any(mj >= 2)
  nb <- sum(mj)
  bitem <- rep(seq_len(J), mj)
  blev <- unlist(lapply(mj, seq_len))
  plab <- c("zero", "partial", "full")
  grp <- matrix(NA_character_, J, 3)
  for (jj in seq_len(J)) {
    pres <- if (mj[jj] >= 2) 1:3 else c(1L, 3L)
    key <- switch(alpha_by,
                  position = plab[pres],
                  maxscore = paste0("max", mj[jj], "_", plab[pres]),
                  item = paste0(levels(x$item)[jj], "_", plab[pres]))
    grp[jj, pres] <- key
  }
  ulab <- unique(as.vector(t(grp)))
  ulab <- ulab[!is.na(ulab)]
  if (alpha_by == "position") ulab <- intersect(plab, ulab)
  agrid <- matrix(match(grp, ulab), J, 3)
  apos <- matrix(NA_integer_, n, K)
  for (k in seq_len(K)) {
    v <- valid[, k]
    apos[v, k] <- agrid[cbind(j[v], pos[v, k])]
  }
  alab <- paste0("alpha_", ulab)
  na <- length(alab)
  P <- J + nb + na + 1 + if (has_partial) 2 else 0
  idx <- list(a = seq_len(J), b = J + seq_len(nb),
              alpha = J + nb + seq_len(na), gamma = J + nb + na + 1,
              delta = if (has_partial) J + nb + na + 1 + 1:2 else integer(0))
  nm <- c(paste0("a_", levels(x$item)), paste0("b_", levels(x$item)[bitem], "_", blev),
          alab, "gamma", if (has_partial) c("delta_first", "delta_second"))
  list(n = n, J = J, K = K, j = j, mj = mj, cc = cc, valid = valid, pos = pos,
       rows = split(seq_len(n), factor(j, levels = seq_len(J))),
       pos1 = !is.na(pos) & pos == 1L, pos2 = !is.na(pos) & pos == 2L,
       pos3 = !is.na(pos) & pos == 3L,
       s1 = x$first, s2 = x$second, cluster = x$cluster, has_partial = has_partial,
       bitem = bitem, blev = blev, P = P, idx = idx, names = nm,
       alab = alab, agrid = agrid, apos = apos, alpha_by = alpha_by)
}

# Map position codes to alpha index (handles absence of partial category).
.rm_alpha <- function(S, par) {
  par[S$idx$alpha]
}

.rm_core <- function(S, par, z, scores = TRUE) {
  J <- S$J
  K <- S$K
  a <- par[S$idx$a]
  bm <- matrix(0, J, K)
  bm[cbind(S$bitem, S$blev + 1L)] <- par[S$idx$b]
  cumB <- t(apply(bm, 1, cumsum))
  az <- a[S$j] * z
  eta <- az * S$cc - cumB[S$j, , drop = FALSE]
  eta[!S$valid] <- -Inf
  mx <- eta[, 1]
  for (k in 2:K) mx <- pmax(mx, eta[, k])
  pc <- exp(eta - mx)
  pc <- pc / rowSums(pc)
  al <- .rm_alpha(S, par)
  g <- par[S$idx$gamma]
  e <- stats::plogis(matrix(al[S$apos], S$n) + g * z)
  e[!S$valid] <- 0
  dn <- function(r) {
    dv <- if (S$has_partial) stats::plogis(par[S$idx$delta[r]]) else 0.5
    out <- matrix(dv, S$n, K)
    out[S$cc == 0] <- 0
    out[S$pos3] <- 1
    out
  }
  d1 <- dn(1)
  d2 <- dn(2)
  ps <- function(s, d) {
    df <- s - S$cc
    out <- (df == 0) * (1 - e) + (df == -1) * (e * d) + (df == 1) * (e * (1 - d))
    out[!S$valid] <- 0
    out
  }
  L <- pc * ps(S$s1, d1) * ps(S$s2, d2)
  Ls <- rowSums(L)
  ll <- log(Ls)
  if (!scores) return(list(ll = ll, pc = pc))
  pi <- L / Ls
  Sm <- matrix(0, S$n, S$P)
  # item slope
  Ec <- rowSums(pc * S$cc)
  Ep <- rowSums(pi * S$cc)
  Sm[cbind(seq_len(S$n), S$j)] <- z * (Ep - Ec)
  # item steps: d/db_il = -P_pi(C >= l) + P_model(C >= l)
  tailp <- pc
  tailpi <- pi
  for (k in (K - 1):1) {
    tailp[, k] <- tailp[, k] + tailp[, k + 1]
    tailpi[, k] <- tailpi[, k] + tailpi[, k + 1]
  }
  dtail <- tailp - tailpi
  for (k in seq_along(S$bitem)) {
    rows <- S$rows[[S$bitem[k]]]
    Sm[rows, S$idx$b[k]] <- dtail[rows, S$blev[k] + 1]
  }
  # error part
  u <- ((S$s1 != S$cc) - e) + ((S$s2 != S$cc) - e)
  u[!S$valid] <- 0
  w <- pi * u
  a0 <- S$idx$alpha[1] - 1L
  for (k in seq_len(K)) {
    v <- which(S$valid[, k])
    ix <- cbind(v, a0 + S$apos[v, k])
    Sm[ix] <- Sm[ix] + w[v, k]
  }
  Sm[, S$idx$gamma] <- rowSums(w) * z
  if (S$has_partial) {
    part <- S$pos == 2 & S$valid
    part[is.na(part)] <- FALSE
    for (r in 1:2) {
      s <- if (r == 1) S$s1 else S$s2
      d <- if (r == 1) d1 else d2
      v <- part * (s != S$cc) * ((s < S$cc) - d)
      Sm[, S$idx$delta[r]] <- rowSums(pi * v)
    }
  }
  list(ll = ll, scores = Sm)
}

.rm_start <- function(S, z) {
  par <- numeric(S$P)
  par[S$idx$a] <- 1
  bs <- numeric(length(S$bitem))
  for (k in seq_along(S$bitem)) {
    sj <- S$s1[S$j == S$bitem[k]]
    l <- S$blev[k]
    bs[k] <- log((sum(sj == l - 1) + 0.5) / (sum(sj == l) + 0.5))
  }
  par[S$idx$b] <- bs
  dis <- mean(S$s1 != S$s2)
  par[S$idx$alpha] <- stats::qlogis(min(max(dis / 2, 0.005), 0.3))
  par
}

.rm_pen <- 1 / 100

.rm_fit <- function(S, z, start = NULL, maxit = 200, tol = 1e-8) {
  par <- if (is.null(start)) .rm_start(S, z) else start
  bi <- c(S$idx$b, S$idx$alpha)
  pll <- function(p) sum(.rm_core(S, p, z, scores = FALSE)$ll) - 0.5 * .rm_pen * sum(p[bi]^2)
  cur <- .rm_core(S, par, z)
  f0 <- sum(cur$ll) - 0.5 * .rm_pen * sum(par[bi]^2)
  conv <- FALSE
  for (it in seq_len(maxit)) {
    gr <- colSums(cur$scores)
    gr[bi] <- gr[bi] - .rm_pen * par[bi]
    A <- crossprod(cur$scores)
    diag(A)[bi] <- diag(A)[bi] + .rm_pen
    step <- tryCatch(solve(A, gr), error = function(e) qr.solve(A + diag(1e-8, S$P), gr))
    lam <- 1
    repeat {
      np <- par + lam * step
      f1 <- pll(np)
      if (is.finite(f1) && f1 >= f0 - 1e-10) break
      lam <- lam / 2
      if (lam < 1e-6) break
    }
    if (lam < 1e-6) {
      conv <- max(abs(gr)) < 1e-3
      break
    }
    par <- np
    df <- f1 - f0
    f0 <- f1
    cur <- .rm_core(S, par, z)
    if (abs(df) < tol * (abs(f0) + tol) && max(abs(lam * step)) < 1e-5) {
      conv <- TRUE
      break
    }
  }
  # Newton polishing with the observed information (numerical Hessian of the
  # analytic penalized gradient); the same matrix is the sandwich bread.
  gpen <- function(p) {
    g <- colSums(.rm_core(S, p, z)$scores)
    g[bi] <- g[bi] - .rm_pen * p[bi]
    g
  }
  info <- function(p, g0 = gpen(p)) {
    H <- vapply(seq_along(p), function(k) {
      h <- 1e-6 * max(1, abs(p[k]))
      e <- replace(numeric(length(p)), k, h)
      (gpen(p + e) - g0) / h
    }, numeric(length(p)))
    I <- -(H + t(H)) / 2
    ev <- eigen(I, symmetric = TRUE)
    if (min(ev$values) <= 1e-8 * max(ev$values)) {
      I <- ev$vectors %*% diag(pmax(ev$values, 1e-8 * max(ev$values))) %*% t(ev$vectors)
    }
    I
  }
  I <- info(par)
  par0 <- par
  dec <- Inf
  for (nt in seq_len(25)) {
    gr <- gpen(par)
    step <- solve(I, gr)
    dec <- sum(gr * step)
    if (dec < 1e-10) break
    lam <- 1
    repeat {
      np <- par + lam * step
      f1 <- pll(np)
      if (is.finite(f1) && f1 >= f0 - 1e-12) break
      lam <- lam / 2
      if (lam < 1e-8) break
    }
    if (lam < 1e-8) break
    par <- np
    f0 <- f1
  }
  cur <- .rm_core(S, par, z)
  gr <- gpen(par)
  if (max(abs(par - par0)) > 1e-3) I <- info(par, gr)
  dec <- sum(gr * solve(I, gr))
  conv <- is.finite(dec) && dec < 1e-6
  Ii <- solve(I)
  Sm <- cur$scores
  Sg <- rowsum(Sm, S$cluster)
  G <- nrow(Sg)
  V <- Ii %*% (crossprod(Sg) * G / (G - 1)) %*% Ii
  names(par) <- S$names
  dimnames(V) <- list(S$names, S$names)
  list(par = par, V = V, loglik = sum(cur$ll), iterations = it, converged = conv,
       decrement = dec)
}

# Implied item-averaged disagreement at ability points zz.
.rm_curve <- function(S, par, zz, wj, zero_gamma = FALSE) {
  if (zero_gamma) par[S$idx$gamma] <- 0
  J <- S$J
  K <- S$K
  a <- par[S$idx$a]
  bm <- matrix(0, J, K)
  bm[cbind(S$bitem, S$blev + 1L)] <- par[S$idx$b]
  cumB <- t(apply(bm, 1, cumsum))
  al <- .rm_alpha(S, par)
  g <- par[S$idx$gamma]
  dlt <- if (S$has_partial) stats::plogis(par[S$idx$delta]) else c(0.5, 0.5)
  out <- numeric(length(zz))
  for (q in seq_along(zz)) {
    tot <- 0
    for (jj in seq_len(J)) {
      m <- S$mj[jj]
      cc <- 0:m
      eta <- a[jj] * zz[q] * cc - cumB[jj, cc + 1]
      p <- exp(eta - max(eta))
      p <- p / sum(p)
      ps <- ifelse(cc == 0, 1, ifelse(cc == m, 3, 2))
      e <- stats::plogis(al[S$agrid[jj, ps]] + g * zz[q])
      same <- ifelse(ps == 2, dlt[1] * dlt[2] + (1 - dlt[1]) * (1 - dlt[2]), 1)
      agr <- (1 - e)^2 + e^2 * same
      tot <- tot + wj[jj] * (1 - sum(p * agr))
    }
    out[q] <- tot
  }
  out
}

# Implied disagreement probability for each observed response.
.rm_rowdis <- function(S, par, z) {
  pc <- .rm_core(S, par, z, scores = FALSE)$pc
  al <- .rm_alpha(S, par)
  e <- stats::plogis(matrix(al[S$apos], S$n) + par[S$idx$gamma] * z)
  e[!S$valid] <- 0
  dlt <- if (S$has_partial) stats::plogis(par[S$idx$delta]) else c(0.5, 0.5)
  same <- ifelse(S$pos == 2, dlt[1] * dlt[2] + (1 - dlt[1]) * (1 - dlt[2]), 1)
  agr <- (1 - e)^2 + e^2 * same
  agr[!S$valid] <- 0
  1 - rowSums(pc * agr)
}

# Multivariate normal draws without external dependencies.
.mvdraw <- function(mu, V, n) {
  ev <- eigen((V + t(V)) / 2, symmetric = TRUE)
  R <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), length(mu))
  sweep(matrix(stats::rnorm(n * length(mu)), n) %*% t(R), 2, mu, `+`)
}
