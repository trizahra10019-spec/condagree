# Usage: Rscript verify_pirls_independent.R <data folder> <pirls output folder>
# Independent re-computation of key PIRLS numbers (no condagree functions)
suppressMessages({library(haven); library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE)  # <data folder> <pirls output folder>
fr <- function(p) list.files(args[1], p, full.names = TRUE, ignore.case = TRUE)[1]
r <- read_sav(fr("asridnr3\\.sav$")); a <- read_sav(fr("asaidnr3\\.sav$"))
J <- fromJSON(file.path(args[2], "pirls_results.json"))
rv <- grep("^R[123]R[A-Z]", names(r), value=TRUE)
L <- do.call(rbind, lapply(rv, function(v) data.frame(id=as.numeric(r$IDSTUD), sch=as.numeric(r$IDSCHOOL), item=v,
  s2=as.numeric(zap_labels(r[[v]])), s1=as.numeric(zap_labels(r[[sub("^(R[123])R","\\11",v)]])))))
L <- L[complete.cases(L),]
pv <- sapply(1:5, function(k) as.numeric(a[[paste0("ASRREA0",k)]])[match(L$id, as.numeric(a$IDSTUD))])
chk <- function(name, mine, pkg) cat(sprintf("%-38s independent=%.6f  package=%.6f  %s\n", name, mine, pkg, ifelse(abs(mine-pkg)<1e-6,"OK","DIFF")))
chk("n responses", nrow(L), J$n_responses)
chk("global exact agreement", mean(L$s1==L$s2), J$ca_bands$global[1])
# item-stratified pooled kappa
ks <- function(d){ po <- mean(d$s1==d$s2); pe <- sum(sapply(split(d,d$item), function(x){ lv<-union(x$s1,x$s2)
  sum(sapply(lv, function(k) mean(x$s1==k)*mean(x$s2==k)))*nrow(x)}))/nrow(d); (po-pe)/(1-pe) }
chk("global item-stratified kappa", ks(L), J$ca_bands$global[2])
# quintile kappa pooled over PVs (cut on pooled student-level PV distribution)
u <- !duplicated(L$id); br <- quantile(as.vector(pv[u,]), 0:5/5); br[1] <- -Inf; br[6] <- Inf
kq <- sapply(1:5, function(m){ b <- cut(pv[,m], br, labels=FALSE); sapply(1:5, function(q) ks(L[b==q,])) })
aq <- sapply(1:5, function(m){ b <- cut(pv[,m], br, labels=FALSE); sapply(1:5, function(q) mean(L$s1[b==q]==L$s2[b==q])) })
for (q in c(1,5)) { chk(paste0("kappa Q",q," (mean over PVs)"), mean(kq[q,]), J$ca_bands$curve$kappa[q]); chk(paste0("exact Q",q), mean(aq[q,]), J$ca_bands$curve$agree[q]) }
# trend slope via glm with item fixed effects, averaged over PVs
sc <- c(mean(as.vector(pv[u,])), sd(as.vector(pv[u,])))
y <- as.numeric(L$s1!=L$s2); inf <- names(which(tapply(y,L$item,sum)>0)); k <- L$item %in% inf
sl <- sapply(1:5, function(m) coef(glm(y[k] ~ factor(L$item[k]) + I((pv[k,m]-sc[1])/sc[2]), family=binomial))[length(inf)+1])
chk("trend slope per SD (mean over PVs)", mean(sl), J$agree_test_pv$slope_sd[1])
chk("odds ratio per 100 points", exp(mean(sl)*100/sc[2]), J$agree_test_pv$or_unit[1])
chk("n disagreements", sum(y), J$n_disagree)
