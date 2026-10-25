# Summarize sim_item_heterogeneity.R output
# Usage: Rscript sim_item_heterogeneity_summarize.R <folder>
din <- commandArgs(trailingOnly = TRUE)[1]
rows <- lapply(list.files(din, "^het_gamma_.*\\.rds$", full.names = TRUE), function(f) {
  o <- readRDS(f); r <- as.data.frame(o$results)
  do.call(rbind, lapply(c("pos", "item"), function(m) {
    g <- r[[paste0(m, "_g")]]; ok <- !is.na(g)
    data.frame(gamma = o$gamma, model = if (m == "pos") "position" else "item",
               reps = sum(ok), mean = mean(g[ok]), bias = mean(g[ok]) - o$gamma,
               rmse = sqrt(mean((g[ok] - o$gamma)^2)),
               se_ratio = mean(r[[paste0(m, "_se")]][ok]) / stats::sd(g[ok]),
               reject = mean(r[[paste0(m, "_p")]][ok] < 0.05),
               cover = mean(r[[paste0(m, "_lo")]][ok] <= o$gamma & r[[paste0(m, "_hi")]][ok] >= o$gamma),
               converged = mean(r[[paste0(m, "_conv")]][ok] == 1))
  }))
})
out <- do.call(rbind, rows)
utils::write.csv(out, file.path(din, "het_summary.csv"), row.names = FALSE)
print(out, digits = 3)
