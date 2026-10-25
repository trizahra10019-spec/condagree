# Re-fit rater_model() in every replication of an existing simulation run
# (same seeds, hence identical data) and replace its columns. Used after the
# switch from an outer-product to an observed-information sandwich.
# Usage: Rscript sim_rerun_rater_model.R <in folder> <out folder> <condition ids>
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
    rm <- tryCatch(rater_model(x), error = function(e) NULL)
    if (is.null(rm)) { o$results[r, grep("^rm_", colnames(o$results))] <- NA; next }
    co <- rm$coefficients
    g <- co[co$parameter == "gamma", ]; sv <- co[co$parameter == "severity", ]
    sh <- co[co$parameter == "share_opportunity", ]
    o$results[r, c("rm_gamma", "rm_gamma_se", "rm_gamma_lo", "rm_gamma_hi", "rm_gamma_p",
                   "rm_sev", "rm_sev_p", "rm_share", "rm_conv")] <-
      c(g$estimate, g$se, g$lower, g$upper, g$p, sv$estimate, sv$p, sh$estimate,
        all(rm$converged))
  }
  o$seconds_rm <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  saveRDS(o, fo)
  cat(sprintf("condition %d re-fitted in %.1f min\n", i, o$seconds_rm / 60))
}
