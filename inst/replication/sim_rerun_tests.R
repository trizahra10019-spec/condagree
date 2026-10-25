# Re-run agree_test() and severity_test() in every replication of an existing
# simulation run (same seeds, hence identical data) and replace their columns.
# Used after adding small-sample (G - 1) reference distributions.
# Usage: Rscript sim_rerun_tests.R <in folder> <out folder> <condition ids>
suppressMessages(library(condagree))
args <- commandArgs(trailingOnly = TRUE)
din <- args[1]; dout <- args[2]
ids <- as.integer(strsplit(args[3], ",")[[1]])
dir.create(dout, showWarnings = FALSE, recursive = TRUE)
ALPHA_OPP <- c(zero = -4.5, partial = -2.2, full = -3.5)
ALPHA_FLAT <- c(zero = -3.2, partial = -3.2, full = -3.2)
for (i in ids) {
  fo <- file.path(dout, sprintf("cond_%02d.rds", i))
  if (file.exists(fo)) next
  o <- readRDS(file.path(din, sprintf("cond_%02d.rds", i)))
  cd <- o$design
  n_per <- if (cd$N == 300) 10 else 20
  t0 <- Sys.time()
  for (r in seq_len(o$reps)) {
    d <- simulate_agree(n_schools = cd$N / n_per, n_per_school = n_per, gamma = cd$gamma,
                        alpha = if (cd$opportunity) ALPHA_OPP else ALPHA_FLAT,
                        reliability = if (cd$measure == "pv") 0.90 else 1,
                        n_pv = if (cd$measure == "pv") 5 else 1, seed = 1e5 * i + r)
    stopifnot(abs(mean(d$first != d$second) - o$results[r, "disagree"]) < 1e-12)
    pv <- grep("^PV", names(d), value = TRUE)
    x <- ca_data(d, "first", "second", "item", pv, cluster = "school", id = "student")
    ft <- agree_test(x, adjust = c("none", "score"))$results
    st <- severity_test(x)$joint
    o$results[r, c("fe_slope", "fe_trend_p", "fe_omni_p", "naive_slope", "naive_trend_p",
                   "naive_omni_p", "sev_joint_p")] <-
      c(ft$slope_sd[1], ft$trend_p[1], ft$omnibus_p[1], ft$slope_sd[2], ft$trend_p[2],
        ft$omnibus_p[2], unname(st["p"]))
  }
  saveRDS(o, fo)
  cat(sprintf("condition %d tests re-run in %.1f min\n", i,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
