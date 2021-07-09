library('dplyr')
source('./R/merge_kuperman.R')
source('./R/merge_SUBTLEX.R')

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
                           within = .("RESP_ID"),
                           between = .("COND"),
                           type = 2), list(dv = dv))))
}


x <- c("Nletters", "Nphon", "Nsyll", "AoA_Kup", "Lg10WF", "Lg10CD")
anova_list <- lapply(x, run_anova, .data = assoc_resp_stats_avg)
names(anova_list) <- x
lapply(anova_list, apaTables::apa.ezANOVA.table)

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
