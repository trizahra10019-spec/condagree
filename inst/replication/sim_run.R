# Simulation study for condagree
# Usage: Rscript sim_run.R <outdir> <reps> <condition ids, comma-separated>
# Each condition is saved as <outdir>/cond_XX.rds; completed conditions are
# skipped, so the script can be restarted in chunks.
suppressMessages(library(condagree))

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "sim_out"
reps <- if (length(args) >= 2) as.integer(args[2]) else 200L
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

design <- expand.grid(
  N = c(300, 1200),
  gamma = c(0, 0.15, 0.30),
  opportunity = c(FALSE, TRUE),
  measure = c("true", "pv"),
  stringsAsFactors = FALSE)
design$id <- seq_len(nrow(design))
ids <- if (length(args) >= 3) as.integer(strsplit(args[3], ",")[[1]]) else design$id

ALPHA_OPP <- c(zero = -4.5, partial = -2.2, full = -3.5)
ALPHA_FLAT <- c(zero = -3.2, partial = -3.2, full = -3.2)
PTS <- c(-1.2816, 0, 1.2816)
H <- 0.3
REL <- 0.90
B <- 199

gen <- function(cd, n_schools, n_per, seed) {
  simulate_agree(n_schools = n_schools, n_per_school = n_per,
                 gamma = cd$gamma,
                 alpha = if (cd$opportunity) ALPHA_OPP else ALPHA_FLAT,
                 reliability = if (cd$measure == "pv") REL else 1,
                 n_pv = if (cd$measure == "pv") 5 else 1, seed = seed)
}
mk <- function(d) {
  pv <- grep("^PV", names(d), value = TRUE)
  ca_data(d, "first", "second", "item", pv, cluster = "school", id = "student")
}

one_rep <- function(cd, seed) {
  n_per <- if (cd$N == 300) 10 else 20
  d <- gen(cd, cd$N / n_per, n_per, seed)
  x <- mk(d)
  ft <- agree_test(x, adjust = c("none", "score"))$results
  rm <- rater_model(x, ndraws = 100)
  co <- rm$coefficients
  g <- co[co$parameter == "gamma", ]
  sv <- co[co$parameter == "severity", ]
  sh <- co[co$parameter == "share_opportunity", ]
  st <- severity_test(x)$joint
  ca <- cond_agree(x, at = PTS, bandwidth = H, contrast_at = PTS[c(1, 3)], B = B)
  cv <- ca$curve
  ck <- ca$contrast[ca$contrast$statistic == "kappa", ]
  c(disagree = mean(d$first != d$second),
    fe_slope = ft$slope_sd[1], fe_trend_p = ft$trend_p[1], fe_omni_p = ft$omnibus_p[1],
    naive_slope = ft$slope_sd[2], naive_trend_p = ft$trend_p[2],
    naive_omni_p = ft$omnibus_p[2],
    rm_gamma = g$estimate, rm_gamma_se = g$se, rm_gamma_lo = g$lower,
    rm_gamma_hi = g$upper, rm_gamma_p = g$p,
    rm_sev = sv$estimate, rm_sev_p = sv$p, rm_share = sh$estimate,
    rm_conv = all(rm$converged),
    sev_joint_p = unname(st["p"]),
    kap1 = cv$kappa[1], kap2 = cv$kappa[2], kap3 = cv$kappa[3],
    kap1_lo = cv$kappa_lo[1], kap2_lo = cv$kappa_lo[2], kap3_lo = cv$kappa_lo[3],
    kap1_hi = cv$kappa_hi[1], kap2_hi = cv$kappa_hi[2], kap3_hi = cv$kappa_hi[3],
    ag1 = cv$agree[1], ag2 = cv$agree[2], ag3 = cv$agree[3],
    ag1_lo = cv$agree_lo[1], ag2_lo = cv$agree_lo[2], ag3_lo = cv$agree_lo[3],
    ag1_hi = cv$agree_hi[1], ag2_hi = cv$agree_hi[2], ag3_hi = cv$agree_hi[3],
    kc = ck$difference, kc_lo = ck$lower, kc_hi = ck$upper)
}

# Population values of the kernel estimands (same bandwidth and points),
# from one very large sample per data-generating process.
pop_truth <- function(cd) {
  d <- gen(cd, 2000, 50, seed = 999000 + cd$id)
  x <- mk(d)
  ca <- cond_agree(x, at = PTS, bandwidth = H, contrast_at = PTS[c(1, 3)], B = 0)
  c(kap = ca$curve$kappa, ag = ca$curve$agree,
    kc = ca$contrast$difference[ca$contrast$statistic == "kappa"])
}

NM <- c("disagree", "fe_slope", "fe_trend_p", "fe_omni_p", "naive_slope",
        "naive_trend_p", "naive_omni_p", "rm_gamma", "rm_gamma_se", "rm_gamma_lo",
        "rm_gamma_hi", "rm_gamma_p", "rm_sev", "rm_sev_p", "rm_share", "rm_conv",
        "sev_joint_p", "kap1", "kap2", "kap3", "kap1_lo", "kap2_lo", "kap3_lo",
        "kap1_hi", "kap2_hi", "kap3_hi", "ag1", "ag2", "ag3", "ag1_lo", "ag2_lo",
        "ag3_lo", "ag1_hi", "ag2_hi", "ag3_hi", "kc", "kc_lo", "kc_hi")

for (i in ids) {
  f <- file.path(outdir, sprintf("cond_%02d.rds", i))
  if (file.exists(f)) next
  cd <- design[design$id == i, ]
  t0 <- Sys.time()
  truth <- pop_truth(cd)
  res <- t(vapply(seq_len(reps), function(r) {
    out <- tryCatch(one_rep(cd, seed = 1e5 * i + r), error = function(e) NULL)
    if (is.null(out)) rep(NA_real_, length(NM)) else unname(out)
  }, numeric(length(NM))))
  colnames(res) <- NM
  saveRDS(list(design = cd, truth = truth, results = res, reps = reps,
               seconds = as.numeric(difftime(Sys.time(), t0, units = "secs"))), f)
  cat(sprintf("condition %d done in %.1f min\n", i,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
