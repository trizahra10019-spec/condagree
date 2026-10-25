# Supplementary figure from saved PIRLS results
# Usage: Rscript supplement_figures.R <pirls output folder>
suppressMessages(library(jsonlite))
dout <- commandArgs(trailingOnly = TRUE)[1]
r <- fromJSON(file.path(dout, "pirls_results.json"))
pv <- r$ca_bands$curve; mc <- r$ca_bands_mc$curve; fc <- r$rater_model_pv$fit_check
png(file.path(dout, "figureS1.png"), width = 2400, height = 1100, res = 300)
op <- par(mfrow = c(1, 2), mar = c(4.2, 4.4, 1.6, 0.8), cex = 0.8)
yl <- range(c(pv$kappa_lo, pv$kappa_hi, mc$kappa_lo, mc$kappa_hi))
plot(1:5 - 0.08, pv$kappa, ylim = yl, pch = 19, xaxt = "n", las = 1,
     xlab = "Achievement quintile", ylab = "Conditional kappa", xlim = c(0.6, 5.4))
arrows(1:5 - 0.08, pv$kappa_lo, 1:5 - 0.08, pv$kappa_hi, angle = 90, code = 3, length = 0.03)
points(1:5 + 0.08, mc$kappa, pch = 1)
arrows(1:5 + 0.08, mc$kappa_lo, 1:5 + 0.08, mc$kappa_hi, angle = 90, code = 3, length = 0.03, lty = 2)
axis(1, at = 1:5, labels = paste0("Q", 1:5))
legend("bottomleft", c("Plausible values", "Multiple-choice score"), pch = c(19, 1),
       lty = 1:2, bty = "n")
mtext("A", side = 3, adj = 0, line = 0.3, font = 2)
plot(fc$quintile, 100 * fc$observed, pch = 19, las = 1, xaxt = "n", xlab = "Achievement quintile (PV 1)",
     ylab = "Disagreement (%)", ylim = c(0, 100 * max(fc$observed, fc$implied) * 1.1))
lines(fc$quintile, 100 * fc$implied, lwd = 2)
axis(1, at = 1:5, labels = paste0("Q", 1:5))
legend("topleft", c("Observed", "Implied by rater_model()"), pch = c(19, NA), lty = c(NA, 1),
       lwd = c(NA, 2), bty = "n")
mtext("B", side = 3, adj = 0, line = 0.3, font = 2)
par(op); invisible(dev.off())
