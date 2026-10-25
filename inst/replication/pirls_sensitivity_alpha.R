# Sensitivity of the residual ability effect to the grouping of the error
# intercepts in rater_model() (position, maximum score x position, item x
# position), for plausible values and the multiple-choice score.
# Usage: Rscript pirls_sensitivity_alpha.R <data folder> <pirls output folder>
suppressMessages({library(condagree); library(haven); library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE)
din <- args[1]; dout <- args[2]
fr <- function(p) list.files(din, p, full.names = TRUE, ignore.case = TRUE)[1]
r <- read_sav(fr("asridnr3\\.sav$")); a <- read_sav(fr("asaidnr3\\.sav$"))
rel <- grep("^R[123]R[A-Z]", names(r), value = TRUE)
long <- do.call(rbind, lapply(rel, function(v) data.frame(
  IDSTUD = as.numeric(r$IDSTUD), IDSCHOOL = as.numeric(r$IDSCHOOL), item = substr(v, 4, nchar(v)),
  second = as.numeric(zap_labels(r[[v]])), first = as.numeric(zap_labels(r[[sub("^(R[123])R", "\\11", v)]])))))
long <- long[!is.na(long$first) & !is.na(long$second) & long$first <= 3 & long$second <= 3, ]
pv <- paste0("ASRREA0", 1:5)
ach <- data.frame(IDSTUD = as.numeric(a$IDSTUD), IDBOOK = as.numeric(a$IDBOOK), sapply(pv, function(v) as.numeric(a[[v]])))
mc <- grep("^R[123]1[A-Z0-9]+M$", names(a), value = TRUE)
sc <- sapply(mc, function(v) { lab <- attr(a[[v]], "labels"); key <- lab[grepl("\\*", names(lab))]
  as.numeric(zap_labels(a[[v]])) == as.numeric(key) })
adm <- sapply(mc, function(v) tapply(!is.na(a[[v]]), a$IDBOOK, any))
admin <- adm[as.character(as.numeric(a$IDBOOK)), , drop = FALSE]
sc[is.na(sc) & admin] <- FALSE; sc[!admin] <- NA
ach$mc_z <- ave(rowMeans(sc, na.rm = TRUE), ach$IDBOOK, FUN = function(v) (v - mean(v)) / stats::sd(v))
dat <- merge(long, ach[, c("IDSTUD", pv, "mc_z")], by = "IDSTUD")
x_pv <- ca_data(dat, "first", "second", "item", pv, cluster = "IDSCHOOL", id = "IDSTUD")
x_mc <- ca_data(dat, "first", "second", "item", "mc_z", cluster = "IDSCHOOL", id = "IDSTUD")
out <- list()
for (ab in c("position", "maxscore", "item")) for (cv in c("pv", "mc")) {
  t0 <- Sys.time()
  rm <- rater_model(if (cv == "pv") x_pv else x_mc, alpha_by = ab)
  co <- rm$coefficients
  g <- co[co$parameter == "gamma", ]; sh <- co[co$parameter == "share_opportunity", ]
  sv <- co[co$parameter == "severity", ]
  out[[paste(ab, cv)]] <- data.frame(alpha_by = ab, conditioning = cv,
    n_alpha = nrow(rm$alpha), gamma = g$estimate, se = g$se, lower = g$lower, upper = g$upper,
    p = g$p, share = sh$estimate, share_lo = sh$lower, share_hi = sh$upper,
    severity = sv$estimate, severity_p = sv$p, converged = all(rm$converged),
    loglik = mean(rm$loglik), npar = rm$npar, aic = rm$aic, minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  print(out[[paste(ab, cv)]])
}
res <- do.call(rbind, out); rownames(res) <- NULL
write.csv(res, file.path(dout, "pirls_alpha_sensitivity.csv"), row.names = FALSE)
