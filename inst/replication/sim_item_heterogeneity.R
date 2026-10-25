# Supplementary simulation: item-specific scoring difficulty correlated with
# item location (harder items are also harder to score). Compares
# rater_model() with position-level and item-level error intercepts.
# Usage: Rscript sim_item_heterogeneity.R <out folder> <gamma> <reps>
suppressMessages(library(condagree))
args <- commandArgs(trailingOnly = TRUE)
dout <- args[1]; gam <- as.numeric(args[2]); reps <- as.integer(args[3])
dir.create(dout, showWarnings = FALSE, recursive = TRUE)
J <- 12
shift <- 1.2 * seq(-1, 1, length.out = J)
res <- t(vapply(seq_len(reps), function(r) {
  d <- simulate_agree(n_schools = 60, n_per_school = 20, gamma = gam,
                      alpha = c(zero = -4.5, partial = -2.2, full = -3.5),
                      alpha_item = shift, reliability = 1, n_pv = 1,
                      seed = 7e5 + 1000 * gam + r)
  x <- ca_data(d, "first", "second", "item", "PV1", cluster = "school", id = "student")
  out <- c()
  for (ab in c("position", "item")) {
    rm <- tryCatch(rater_model(x, alpha_by = ab), error = function(e) NULL)
    g <- if (is.null(rm)) rep(NA, 5) else {
      co <- rm$coefficients[rm$coefficients$parameter == "gamma", ]
      c(co$estimate, co$se, co$lower, co$upper, co$p)
    }
    out <- c(out, g, if (is.null(rm)) NA else all(rm$converged))
  }
  out
}, numeric(12)))
colnames(res) <- c(paste0("pos_", c("g", "se", "lo", "hi", "p", "conv")),
                   paste0("item_", c("g", "se", "lo", "hi", "p", "conv")))
saveRDS(list(gamma = gam, shift = shift, results = res), file.path(dout, sprintf("het_gamma_%.2f.rds", gam)))
cat("done", gam, "\n")
