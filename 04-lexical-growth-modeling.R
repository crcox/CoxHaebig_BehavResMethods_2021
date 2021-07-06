library('dplyr')
library('netgrowr')
library('parallel')
source('./R/local_utils.R')

append_formula <- function(f, x) {
    update(f, paste("~.", x, sep = "+"))
}

combn_and_paste <- function(x, m) {
  return(apply(combn(x, m), MARGIN = 2, paste, collapse = "+"))
}

model_comp_helper <- function(full, restricted, M) {
  return(netgrowr::model_comparison(M[[full]], M[[restricted]]))
}
# Load data ----
load('./network/modelvars.Rdata')

# Define cluster for parallel computing ----
ncores <- parallel::detectCores() - 1
cl <- parallel::makeCluster(ncores)
parallel::clusterSetRNGStream(cl)
invisible(parallel::clusterEvalQ(cl, {
  library('dplyr')
  library('netgrowr')
  source('./R/local_utils.R')
}))


# Empty (random selection) model ----
# +++ All words assigned equal probability
fE <- aoa ~ 1
ME <- netgrowr::mle_network_growth(fE, data = na.omit(modelvars), split_by = "month", label_with = "num_item_id")

# Psycholinguistic baseline model ----
# +++ number of phonemes
# +++ phonotactic probability (biphone)
# +++ phonological neighborhood density (KU child corpus)
# +++ CHILDES Frequency
f0 <- update(fE, ~ . + I(Z(log(nphon + 1))) + I(Z(log(CHILDES_Freq + 1))) + Z(BiphonProb.avg) + Z(PNDC.avg))
M0 <- netgrowr::mle_network_growth(f0, data = na.omit(modelvars), split_by = "month", label_with = "num_item_id")
round(netgrowr::model_comparison(M0, ME), 3)

# Networks over control ----
network_gv1 <- c(
    "preferential_attachment_adult",
    "preferential_attachment_child",
    "preferential_attachment_childes",
    "preferential_acquisition_adult",
    "preferential_acquisition_child",
    "preferential_acquisition_childes",
    "lure_of_the_associates_adult",
    "lure_of_the_associates_child",
    "lure_of_the_associates_childes"
)
f1 <- lapply(network_gv1, FUN = append_formula, f = f0)
M1 <- parLapply(cl, f1, netgrowr::mle_network_growth, data = modelvars, split_by = "month", label_with = "num_item_id")
names(M1) <- network_gv1

# Models that combine networks growing by pref acq. ----
network_gv2 <- combn_and_paste(c(
      "preferential_acquisition_adult",
      "preferential_acquisition_child",
      "preferential_acquisition_childes"
    ),
    m = 2
)
f2n <- lapply(network_gv2, FUN = append_formula, f = f0)
M2n <- parallel::parLapply(
  cl = cl,
  X = f2n,
  fun = netgrowr::mle_network_growth,
  data = modelvars,
  split_by = "month",
  label_with = "num_item_id"
)
names(M2n) <- network_gv2

network_gv3 <- combn_and_paste(c(
      "preferential_acquisition_adult",
      "preferential_acquisition_child",
      "preferential_acquisition_childes"
    ),
    m = 3
)
f3n <- lapply(network_gv3, FUN = append_formula, f = f0)
M3n <- parallel::parLapply(
  cl = cl,
  X = f3n,
  fun = netgrowr::mle_network_growth,
  data = modelvars,
  split_by = "month",
  label_with = "num_item_id"
)
names(M3n) <- network_gv3

# Models that combine growth values derived from different growth models for the same network ----
adult_gv2 <- combn_and_paste(c(
      "preferential_attachment_adult",
      "preferential_acquisition_adult",
      "lure_of_the_associates_adult"
    ),
    m = 2
)
child_gv2 <- combn_and_paste(c(
      "preferential_attachment_child",
      "preferential_acquisition_child",
      "lure_of_the_associates_child"
    ),
    m = 2
)
childes_gv2 <- combn_and_paste(c(
      "preferential_attachment_childes",
      "preferential_acquisition_childes",
      "lure_of_the_associates_childes"
    ),
    m = 2
)
f2g <- lapply(
  c(
    adult_gv2,
    child_gv2,
    childes_gv2
  ),
  FUN = append_formula,
  f = f0
)
M2g <- parallel::parLapply(
  cl = cl,
  X = f2g,
  fun = netgrowr::mle_network_growth,
  data = modelvars,
  split_by = "month",
  label_with = "num_item_id"
)
names(M2g) <- c(adult_gv2, child_gv2, childes_gv2)

# Comparison 1: Each growth model (and combinations) over baseline ----
X <- list()
X[[1]] <- t(vapply(c(M1, M2n, M3n), netgrowr::model_comparison, numeric(6), restricted = M0))
round(X[[1]], 3)

# Comparison 2: Preferential Acq. vs. LOA ----
comparisons <- as.data.frame(rbind(
  c(full = 'preferential_acquisition_adult+lure_of_the_associates_adult', restricted = 'preferential_acquisition_adult'),
  c(full = 'preferential_acquisition_adult+lure_of_the_associates_adult', restricted = 'lure_of_the_associates_adult'),
  c(full = 'preferential_acquisition_child+lure_of_the_associates_child', restricted = 'preferential_acquisition_child'),
  c(full = 'preferential_acquisition_child+lure_of_the_associates_child', restricted = 'lure_of_the_associates_child'),
  c(full = 'preferential_acquisition_childes+lure_of_the_associates_childes', restricted = 'preferential_acquisition_childes'),
  c(full = 'preferential_acquisition_childes+lure_of_the_associates_childes', restricted = 'lure_of_the_associates_childes')
))
X[[2]] <- t(mapply(model_comp_helper, comparisons[["full"]], comparisons[["restricted"]], MoreArgs = list(M = c(M1, M2g))))
rownames(X[[2]]) <- apply(comparisons, 1, paste, collapse="|")
round(X[[2]], 3)

# Comparison 3: Adult vs. Child vs. CHILDES (Preferential Acq. only) ----
comparisons <- as.data.frame(rbind(
  c(full = 'preferential_acquisition_adult+preferential_acquisition_child', restricted = 'preferential_acquisition_adult'),
  c(full = 'preferential_acquisition_adult+preferential_acquisition_child', restricted = 'preferential_acquisition_child'),
  c(full = 'preferential_acquisition_child+preferential_acquisition_childes', restricted = 'preferential_acquisition_child'),
  c(full = 'preferential_acquisition_child+preferential_acquisition_childes', restricted = 'preferential_acquisition_childes'),
  c(full = 'preferential_acquisition_adult+preferential_acquisition_childes', restricted = 'preferential_acquisition_adult'),
  c(full = 'preferential_acquisition_adult+preferential_acquisition_childes', restricted = 'preferential_acquisition_childes')
))
X[[3]] <- t(mapply(model_comp_helper, comparisons[["full"]], comparisons[["restricted"]], MoreArgs = list(M = c(M1, M2n))))
rownames(X[[3]]) <- gsub("preferential_acquisition_", "", apply(comparisons, 1, paste, collapse="|"))
round(X[[3]], 3)

# Comparison 4: Adult & Child vs. CHILDES (Preferential Acq. only) ----
comparisons <- as.data.frame(rbind(
  c(full = 'preferential_acquisition_adult+preferential_acquisition_child+preferential_acquisition_childes', restricted = 'preferential_acquisition_adult+preferential_acquisition_child'),
  c(full = 'preferential_acquisition_adult+preferential_acquisition_child+preferential_acquisition_childes', restricted = 'preferential_acquisition_adult+preferential_acquisition_childes'),
  c(full = 'preferential_acquisition_adult+preferential_acquisition_child+preferential_acquisition_childes', restricted = 'preferential_acquisition_child+preferential_acquisition_childes')
))
X[[4]] <- t(mapply(model_comp_helper, comparisons[["full"]], comparisons[["restricted"]], MoreArgs = list(M = c(M2n, M3n))))
rownames(X[[4]]) <- gsub("preferential_acquisition_", "", apply(comparisons, 1, paste, collapse="|"))
round(X[[4]], 3)

model_comparisons <- do.call('rbind', X)
model_comparisons <- cbind(
  model_comparisons,
  p_fdr  = p.adjust(as.vector(model_comparisons[, 'p']), method = 'fdr'),
  p_bonf = p.adjust(as.vector(model_comparisons[, 'p']), method = 'bonferroni'),
  p_holm = p.adjust(as.vector(model_comparisons[, 'p']), method = 'holm')
)
save(model_comparisons, file = './network/model-comparisons.Rdata')
write.csv(model_comparisons, file = './network/model-comparisons.csv')
