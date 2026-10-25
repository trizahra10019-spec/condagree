# Summarize the simulation study (output of sim_run.R)
# Usage: Rscript sim_summarize.R <sim folder> <output csv>
args <- commandArgs(trailingOnly = TRUE)
din <- if (length(args) >= 1) args[1] else "sim_out"
fout <- if (length(args) >= 2) args[2] else file.path(din, "sim_summary.csv")
files <- sort(list.files(din, "^cond_[0-9]+\\.rds$", full.names = TRUE))
rej <- function(p) mean(p < 0.05, na.rm = TRUE)
cov <- function(lo, hi, t) mean(lo <= t & hi >= t, na.rm = TRUE)
rows <- lapply(files, function(f) {
  o <- readRDS(f)
  r <- as.data.frame(o$results)
  d <- o$design
  tr <- o$truth
  ok <- stats::complete.cases(r)
  r <- r[ok, ]
  cov_k <- mean(c(cov(r$kap1_lo, r$kap1_hi, tr["kap1"]),
                  cov(r$kap2_lo, r$kap2_hi, tr["kap2"]),
                  cov(r$kap3_lo, r$kap3_hi, tr["kap3"])))
  cov_a <- mean(c(cov(r$ag1_lo, r$ag1_hi, tr["ag1"]),
                  cov(r$ag2_lo, r$ag2_hi, tr["ag2"]),
                  cov(r$ag3_lo, r$ag3_hi, tr["ag3"])))
  bias_k <- mean(c(mean(r$kap1) - tr["kap1"], mean(r$kap2) - tr["kap2"],
                   mean(r$kap3) - tr["kap3"]))
  data.frame(id = d$id, N = d$N, gamma = d$gamma, opportunity = d$opportunity,
             measure = d$measure, reps = nrow(r), failed = sum(!ok),
             disagree = mean(r$disagree),
             rej_total_trend = rej(r$fe_trend_p), rej_total_omnibus = rej(r$fe_omni_p),
             rej_naive_trend = rej(r$naive_trend_p),
             rej_gamma = rej(r$rm_gamma_p),
             gamma_mean = mean(r$rm_gamma), gamma_bias = mean(r$rm_gamma) - d$gamma,
             gamma_rmse = sqrt(mean((r$rm_gamma - d$gamma)^2)),
             gamma_se_ratio = mean(r$rm_gamma_se) / stats::sd(r$rm_gamma),
             gamma_cover = cov(r$rm_gamma_lo, r$rm_gamma_hi, d$gamma),
             rej_severity_model = rej(r$rm_sev_p), rej_severity_test = rej(r$sev_joint_p),
             share_mean = mean(r$rm_share), converged = mean(r$rm_conv),
             kappa_bias = unname(bias_k), kappa_cover = cov_k, agree_cover = cov_a,
             contrast_true = unname(tr["kc"]), contrast_bias = mean(r$kc) - unname(tr["kc"]),
             contrast_cover = cov(r$kc_lo, r$kc_hi, tr["kc"]),
             minutes = o$seconds / 60)
})
out <- do.call(rbind, rows)
out <- out[order(out$measure, out$opportunity, out$gamma, out$N), ]
utils::write.csv(out, fout, row.names = FALSE)
print(format(out, digits = 3), row.names = FALSE)
