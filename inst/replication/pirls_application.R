# Empirical illustration: PIRLS 2011 Indonesia reliability-scoring data
#
# Data are NOT redistributed. Download the PIRLS 2011 international database
# (SPSS format) from the IEA (https://www.iea.nl/data-tools/repository/pirls)
# and place the Indonesian files asridnr3.sav (reliability scoring) and
# asaidnr3.sav (achievement) in the folder given as first argument.
#
# Usage: Rscript pirls_application.R <data folder> <output folder> [B]
suppressMessages({
  library(condagree)
  library(haven)
  library(jsonlite)
})
args <- commandArgs(trailingOnly = TRUE)
din <- if (length(args) >= 1) args[1] else "."
dout <- if (length(args) >= 2) args[2] else "pirls_out"
B <- if (length(args) >= 3) as.integer(args[3]) else 999L
dir.create(dout, showWarnings = FALSE, recursive = TRUE)
set.seed(20261024)

f_rel <- list.files(din, "asridnr3\\.sav$", full.names = TRUE, ignore.case = TRUE)[1]
f_ach <- list.files(din, "asaidnr3\\.sav$", full.names = TRUE, ignore.case = TRUE)[1]
r <- read_sav(f_rel)
a <- read_sav(f_ach)

## 1. Long format of double-scored constructed responses ------------------
# R1R..., R2R..., R3R... = second (reliability) score; R11..., R21..., R31...
# = first (operational) score of the same response.
rel_vars <- grep("^R[123]R[A-Z]", names(r), value = TRUE)
first_vars <- sub("^(R[123])R", "\\11", rel_vars)
long <- do.call(rbind, lapply(seq_along(rel_vars), function(k) {
  data.frame(IDSTUD = as.numeric(r$IDSTUD), IDSCHOOL = as.numeric(r$IDSCHOOL),
             IDBOOK = as.numeric(r$IDBOOK),
             item = substr(rel_vars[k], 4, nchar(rel_vars[k])),
             second = as.numeric(zap_labels(r[[rel_vars[k]]])),
             first = as.numeric(zap_labels(r[[first_vars[k]]])))
}))
long <- long[!is.na(long$first) & !is.na(long$second), ]
long <- long[long$first <= 3 & long$second <= 3, ]

## 2. Conditioning variables ----------------------------------------------
pv <- paste0("ASRREA0", 1:5)
ach <- data.frame(IDSTUD = as.numeric(a$IDSTUD), IDBOOK = as.numeric(a$IDBOOK),
                  sapply(pv, function(v) as.numeric(a[[v]])))
# Multiple-choice score: items keyed from value labels (key marked "*");
# administered-but-missing responses (omitted/not reached) scored 0;
# proportion correct standardized within booklet (booklets are randomly
# rotated, so booklet groups are randomly equivalent).
mc <- grep("^R[123]1[A-Z0-9]+M$", names(a), value = TRUE)
sc <- sapply(mc, function(v) {
  lab <- attr(a[[v]], "labels")
  key <- lab[grepl("\\*", names(lab))]
  as.numeric(zap_labels(a[[v]])) == as.numeric(key)
})
adm <- sapply(mc, function(v) tapply(!is.na(a[[v]]), a$IDBOOK, any))
admin <- adm[as.character(as.numeric(a$IDBOOK)), , drop = FALSE]
sc[is.na(sc) & admin] <- FALSE
sc[!admin] <- NA
ach$mc_n <- rowSums(!is.na(sc))
ach$mc_pc <- rowMeans(sc, na.rm = TRUE)
ach$mc_z <- ave(ach$mc_pc, ach$IDBOOK, FUN = function(v) (v - mean(v, na.rm = TRUE)) / stats::sd(v, na.rm = TRUE))
dat <- merge(long, ach[, c("IDSTUD", pv, "mc_z", "mc_n")], by = "IDSTUD")

stopifnot(nrow(dat) == nrow(long))
x_pv <- ca_data(dat, "first", "second", "item", pv, cluster = "IDSCHOOL", id = "IDSTUD")
x_mc <- ca_data(dat, "first", "second", "item", "mc_z", cluster = "IDSCHOOL", id = "IDSTUD")

## 3. Analyses ------------------------------------------------------------
ca_k <- cond_agree(x_pv, B = B)                 # kernel curve, all items
ca_b <- cond_agree(x_pv, bands = 5, B = B)      # quintile bands
ca_q <- cond_agree(x_pv, bands = 5, B = B, weights = "quadratic")
at_pv <- agree_test(x_pv, adjust = c("none", "score"), unit = 100)
rm_pv <- rater_model(x_pv, alpha_by = "maxscore")   # AIC-preferred grouping (see pirls_sensitivity_alpha.R)
sv_pv <- severity_test(x_pv)

ca_mc <- cond_agree(x_mc, bands = 5, B = B)
at_mc <- agree_test(x_mc, adjust = c("none", "score"))
rm_mc <- rater_model(x_mc, alpha_by = "maxscore")
sv_mc <- severity_test(x_mc)

## 4. Item-level summary ----------------------------------------------------
u <- !duplicated(dat$IDSTUD)
med <- stats::median(rowMeans(dat[u, pv]))
dat$hi <- rowMeans(dat[, pv]) > med
kap1 <- function(s1, s2) {
  lv <- sort(unique(c(s1, s2)))
  if (length(lv) < 2) return(NA_real_)
  t <- table(factor(s1, lv), factor(s2, lv)) / length(s1)
  pe <- sum(rowSums(t) * colSums(t))
  (sum(diag(t)) - pe) / (1 - pe)
}
items <- do.call(rbind, lapply(split(dat, dat$item), function(d) {
  data.frame(item = d$item[1], max_score = max(d$first, d$second), n = nrow(d),
             disagreements = sum(d$first != d$second),
             agree = mean(d$first == d$second), kappa = kap1(d$first, d$second),
             agree_low = mean(d$first[!d$hi] == d$second[!d$hi]),
             agree_high = mean(d$first[d$hi] == d$second[d$hi]))
}))
items <- items[order(items$item), ]
rownames(items) <- NULL

## 5. Figure ----------------------------------------------------------------
fig <- function(file) {
  grDevices::png(file, width = 2400, height = 1100, res = 300)
  op <- graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.4, 1.6, 0.8), cex = 0.8)
  plot(ca_k, xlab = "Reading achievement (PIRLS scale)",
       ylab = "Conditional kappa", main = "")
  graphics::mtext("A", side = 3, adj = 0, line = 0.3, font = 2)
  plot(rm_pv, xlab = "Reading achievement (PIRLS scale)",
       ylab = "Implied probability of disagreement", main = "")
  graphics::mtext("B", side = 3, adj = 0, line = 0.3, font = 2)
  graphics::par(op)
  grDevices::dev.off()
}
fig(file.path(dout, "figure1.png"))

## 6. Save -----------------------------------------------------------------
res <- list(
  n_students = length(unique(dat$IDSTUD)), n_schools = length(unique(dat$IDSCHOOL)),
  n_responses = nrow(dat), n_items = length(unique(dat$item)),
  items_by_max = as.list(unclass(table(items$max_score))),
  n_per_item = as.list(unclass(summary(items$n))),
  n_disagree = sum(dat$first != dat$second),
  mc_items_per_student = as.list(unclass(summary(ach$mc_n[ach$IDSTUD %in% dat$IDSTUD]))),
  cor_mc_pv = stats::cor(dat$mc_z[u], rowMeans(dat[u, pv])),
  ca_kernel = ca_k[c("curve", "contrast", "global", "bandwidth")],
  ca_bands = ca_b[c("curve", "contrast", "global")],
  ca_bands_quadratic = ca_q[c("curve", "contrast", "global")],
  agree_test_pv = at_pv$results, agree_test_pv_scale = at_pv$scale,
  rater_model_pv = rm_pv[c("coefficients", "alpha", "fit_check", "converged", "loglik", "contrast_z", "scale", "npar", "aic")],
  severity_pv = list(coef = sv_pv$coefficients, joint = as.list(sv_pv$joint),
                     n = sv_pv$n_discordant, prop = sv_pv$prop_second_lower),
  ca_bands_mc = ca_mc[c("curve", "contrast", "global")],
  agree_test_mc = at_mc$results,
  rater_model_mc = rm_mc[c("coefficients", "alpha", "fit_check", "converged", "loglik", "npar", "aic")],
  severity_mc = list(coef = sv_mc$coefficients, joint = as.list(sv_mc$joint)),
  items = items, B = B,
  session = utils::sessionInfo()$R.version$version.string)
saveRDS(list(res = res, rm_curves = rm_pv$curves, ca_kernel_full = ca_k), file.path(dout, "pirls_results.rds"))
write_json(res, file.path(dout, "pirls_results.json"), digits = NA, auto_unbox = TRUE, pretty = TRUE)
utils::write.csv(items, file.path(dout, "pirls_items.csv"), row.names = FALSE)
print(ca_b); print(at_pv); print(rm_pv); print(sv_pv)
cat("\n--- Sensitivity: multiple-choice score ---\n")
print(ca_mc); print(at_mc); print(rm_mc); print(sv_mc)
