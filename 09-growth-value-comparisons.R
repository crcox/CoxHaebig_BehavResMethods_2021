library('apaTables')
library('ez')

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

trim_prefix <- function(x, prefix) {
    return(trimws(x, which = "left", whitespace = prefix))
}

trim_suffix <- function(x, suffix) {
    return(trimws(x, which = "right", whitespace = suffix))
}

add_prefix <- function(x, prefix) {
    str_prepend <- function(prefix, ...) paste(prefix, ..., sep = "")
    return(vapply(x, FUN = str_prepend, FUN.VALUE = character(1), prefix = prefix, USE.NAMES = FALSE))
}

add_suffix <- function(x, suffix) {
    str_append <- function(suffix, ...) paste(..., suffix, sep = "")
    return(vapply(x, FUN = str_append, FUN.VALUE = character(1), suffix = suffix, USE.NAMES = FALSE))
}

load_to_list <- function(files) {
    X <- new.env()
    lapply(files, load, envir = X)
    return(as.list(X))
}

# Load processed word associations ----
load('./data/cdi-metadata-preproc.Rdata')
growth_values <- load_to_list(c(
  "./network/growthvalues-adult-preproc.Rdata",
  "./network/growthvalues-child-preproc.Rdata"
))
names(growth_values) <- trim_suffix(trim_prefix(names(growth_values), "growthvalues_"), "_preproc")

# Combine into one dataframe ----
growth_values$adult$cond <- factor(1, 1:2, c("adult", "child"))
growth_values$child$cond <- factor(2, 1:2, c("adult", "child"))
d <- droplevels(subset(do.call('rbind', growth_values), learned == TRUE & month < 30))

# Report means and one-way t-tests against zero ----
means_tbl <- tapply(d$zscore, d[c("model", "cond")], mean)
print(means_tbl, digits = 3)

stdev_tbl <- tapply(d$zscore, d[c("model", "cond")], sd)
print(means_tbl, digits = 3)

barplot(t(means_tbl[c(1,3,2), ]), beside = TRUE)

ttests_list <- lapply(split(d$zscore, list(d$model, d$cond)), t.test)
ttests_tbl <- t(vapply(ttests_list, ttest_vec, numeric(6)))
print(ttests_tbl, digits = 3)

eza <- ezANOVA(data = d, wid = .(word), dv = .(zscore), within = .(model, cond))
apa.ezANOVA.table(eza)

# Simple effects by model ----
## Paired t-test
simple_effects_by_model <- lapply(split(d, d$model), function(x) {t.test(zscore ~ cond, data = x, paired = TRUE)})
tt_by_model <- t(vapply(simple_effects_by_model, ttest_vec, numeric(6)))
print(tt_by_model, 3)

w <- pivot_wider(d, id_cols = c("model", "word", "month", "aoa"), names_from = "cond", values_from = "zscore")
lapply(split(w, w$model), function(x) abs(mean(x$adult - x$child) / sd(x$adult - x$child)))

# Simple effects by condition ----
## Repeated measures ANOVA
apaTables::apa.ezANOVA.table(
  ez::ezANOVA(subset(d, cond == "adult"),
              dv = .(zscore),
              wid = .(word),
              within = .(model)
  )
)

apaTables::apa.ezANOVA.table(
  ez::ezANOVA(subset(d, cond == "child"),
              dv = .(zscore),
              wid = .(word),
              within = .(model)
  )
)

## Simple effects diving down another level.
t.test(zscore ~ model, data = droplevels(subset(d, model %in% c("preferential_acquisition", "lure_of_the_associates") & cond == "adult")), paired = TRUE)
t.test(zscore ~ model, data = droplevels(subset(d, model %in% c("preferential_acquisition", "lure_of_the_associates") & cond == "child")), paired = TRUE)
