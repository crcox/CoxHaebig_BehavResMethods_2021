library('dplyr')

ttest_vec <- function(x) {
  c(
    tval = as.vector(x[['statistic']]),
    df = as.vector(x[['parameter']]),
    p = as.vector(x[['p.value']]),
    ci.l = as.vector(x[['conf.int']][1]),
    ci.u = as.vector(x[['conf.int']][2]),
    stderr = as.vector(x[['stderr']])
  )
}

# Load response statistics ----
load('./data/assoc-resp-stats.Rdata')

# Begin ----
assoc_resp_stats_avg <- assoc_resp_stats %>%
  group_by(CUE, COND, RESP_ID) %>%
  summarize(
    Nletters = mean(Nletters, na.rm = TRUE),
    Nphon = mean(Nphon, na.rm = TRUE),
    Nsyll = mean(Nsyll, na.rm = TRUE),
    AoA_Kup = mean(AoA_Kup, na.rm = TRUE),
    Lg10WF = mean(Lg10WF, na.rm = TRUE),
    Lg10CD = mean(Lg10CD, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(RESP_ID = as.factor(RESP_ID))

str(assoc_resp_stats_avg)
head(assoc_resp_stats_avg)

run_anova <- function(dv, .data) {
#    f <- as.formula(paste(dv, "RESP_ID + Error(CUE / RESP_ID)", sep = "~"))
#    return(aov(f, data = .data))
    return(eval(substitute(ez::ezANOVA(.data,
                           dv = dv,
                           wid = .("CUE"),
                           within = .("COND", "RESP_ID"),
                           type = 2), list(dv = dv))))
}


x <- c("Nletters", "Nphon", "Nsyll", "AoA_Kup", "Lg10WF", "Lg10CD")
anova_list <- lapply(x, run_anova, .data = assoc_resp_stats_avg)
names(anova_list) <- x
lapply(anova_list, apaTables::apa.ezANOVA.table)
lapply(anova_list, print)

means_tbl <- lapply(x, function(label, data) tapply(data[[label]], data[c('COND', 'RESP_ID')], mean), data = assoc_resp_stats_avg)
names(means_tbl) <- x
print(means_tbl)

stdev_tbl <- lapply(x, function(label, data) {tapply(data[[label]], data[c('COND', 'RESP_ID')], sd)}, data = assoc_resp_stats_avg)
names(stdev_tbl) <- x
print(stdev_tbl)


run_ttest <- function(dv, iv, data) {
    f <- formula(paste(dv, iv, sep = "~"))
    return(t.test(formula = f, data = data, paired = TRUE))
}

x <- c("Nletters", "Nphon", "Nsyll", "AoA_Kup", "Lg10WF", "Lg10CD")
ttest_list <- lapply(x, run_ttest, iv = "COND", data = assoc_resp_stats_avg)
tt_cond_paired <- t(vapply(ttest_list, ttest_vec, numeric(6)))
rownames(tt_cond_paired) <- x
print(tt_cond_paired, digits = 3)


# Paired t-tests and plots ----
# Itemwise analysis comparing adult- and child-oriented responses within cues
# for each response index.
f <- function(x, val_fun = mean) {
  y <- as.data.frame(tapply(x[[1]], list(x[["CUE"]], x[["COND"]]), val_fun))
  return(t.test(Pair(adult, child) ~ 1, data = y))
}
quick_plot <- function(tt, ylab = "", ylim = NULL) {
  ci <- sapply(tt, `[[`, "conf.int")

  if (is.null(ylim)) {
    if (max(ci) < 0)
      ylim <- c(min(ci), 0)
    else if (min(ci) > 0)
      ylim <- c(0, max(ci))
  }

  plot(rep(1:3, each = 2), ci, type = 'n', xlim = c(0.8, 3.2), xlab = "Response Index", ylim = ylim, ylab = ylab, axes = FALSE)
  segments(
    x0 = 1:3,
    y0 = ci[1, ],
    x1 = 1:3,
    y1 = ci[2, ]
  )
  points(1:3, sapply(tt, `[[`, "estimate"), pch = 21, bg = 'white')
  axis(1, at = 1:3)
  axis(2)
}

library('svglite')
svg("./condition-effect-within-cues.svg", height = 4, width = 5)
pp <- par(mfrow = c(1, 4), mar = c(5.1, 4.1, 1.1, 1.1))
d <- get_all_vars(Nletters ~ CUE + COND + RESP_ID, data = assoc_resp_stats_avg)
quick_plot(lapply(split(d, d$RESP_ID), f), "number of letters", ylim = c(0, 0.2))

d <- get_all_vars(AoA_Kup ~ CUE + COND + RESP_ID, data = assoc_resp_stats_avg)
quick_plot(lapply(split(d, d$RESP_ID), f), "age of acquisition", ylim = c(0, 0.5))

d <- get_all_vars(Lg10WF ~ CUE + COND + RESP_ID, data = assoc_resp_stats_avg)
quick_plot(lapply(split(d, d$RESP_ID), f), "log10 word frequency", ylim = c(-.12, 0))

d <- get_all_vars(Lg10CD ~ CUE + COND + RESP_ID, data = assoc_resp_stats_avg)
quick_plot(lapply(split(d, d$RESP_ID), f), "log10 contextual diversity", ylim = c(-.10, 0))
dev.off()

# Plots
