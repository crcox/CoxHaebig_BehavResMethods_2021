library('dplyr')
library('netgrowr')
library('parallel')
source('./R/utils.R')

append_formula <- function(f, x) {
    update(f, paste("~.", x, sep = "+"))
}

combn_and_paste <- function(x, m) {
  return(apply(combn(x, m), MARGIN = 2, paste, collapse = "+"))
}

model_comp_helper <- function(full, restricted, M) {
  return(netgrowr::model_comparison(M[[full]], M[[restricted]]))
}

load_variable <- function(file, variable) {
  e <- new.env()
  load(file, envir = e)
  return(e[[variable]])
}

# Load data ----
# Month 30 is the highest month we have AoA data for.
# If we model month 30, all unknown words will be learned.
# To avoid this, month 30 is not modeled.
modelvars <- subset(load_variable("./network/modelvars.Rdata", "modelvars"),
                    month < 30)

# Define models by constructing formulas ----

# Empty (random selection) model
# +++ All words assigned equal probability
fE <- formula(aoa ~ 1)

# Psycholinguistic baseline model
# +++ number of phonemes
# +++ CHILDES Frequency
# +++ phonotactic probability (biphone)
# +++ phonological neighborhood density (KU child corpus)
f0 <- update(fE, ~ . + I(Z(log(nphon + 1))) + I(Z(log(CHILDES_Freq + 1))) + Z(BiphonProb.avg) + Z(PNDC.avg))
M0 <- netgrowr::mle_network_growth(f0, data = na.omit(modelvars), split_by = "month", label_with = "num_item_id")

# Networks over control
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

# Models that combine networks growing by pref acq.
network_gv2 <- combn_and_paste(c(
      "preferential_acquisition_adult",
      "preferential_acquisition_child",
      "preferential_acquisition_childes"
    ),
    m = 2
)
f2n <- lapply(network_gv2, FUN = append_formula, f = f0)

network_gv3 <- combn_and_paste(c(
      "preferential_acquisition_adult",
      "preferential_acquisition_child",
      "preferential_acquisition_childes"
    ),
    m = 3
)
f3n <- lapply(network_gv3, FUN = append_formula, f = f0)

# Models that combine growth values derived from different growth models for the same network
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

# Define cluster for parallel computing and optimize models ----

# The optimization process is not guaranteed to find a global minimum.
# Different runs of the optimization will yield numeric differences, and
# rarely very bad solutions are obtained. Thus, the optimization process is
# repeated 500 times for every model fit. The minimum negative log likelihood
# for each model over the 500 optimization attempts is retained for the model
# comparisons.
#
# In each repetition, 22 models are optimized. The process can be greatly
# accellerated with parallel computing. We define the 22 formulas defining
# each model in a list, then use parallel::parLapply to fit the models in
# parallel.
#
# This script producees 500 Rdata files in the network/ directory. These
# files are processed and aggregated into a final solution in the next
# script: 05-aggregate-growth-models.R.

ncores <- parallel::detectCores() - 1
cl <- parallel::makeCluster(ncores)
parallel::clusterSetRNGStream(cl)
invisible(parallel::clusterEvalQ(cl, {
  library('dplyr')
  library('netgrowr')
  source('./R/utils.R')
}))


for (repetition in 1:500) {
    M <- parallel::parLapply(
      cl = cl,
      X = c(f1, f2n, f3n, f2g),
      fun = netgrowr::mle_network_growth,
      data = modelvars,
      split_by = "month",
      label_with = "num_item_id"
    )
    names(M) <- c(network_gv1, network_gv2, network_gv3, adult_gv2, child_gv2, childes_gv2)
    M1 <- M[network_gv1]
    M2n <- M[network_gv2]
    M3n <- M[network_gv3]
    M2g <- M[c(adult_gv2, child_gv2, childes_gv2)]

    # Comparison 1: Each growth model (and combinations) over baseline ----
    X <- list()
    X[[1]] <- t(vapply(c(M1, M2n, M3n), netgrowr::model_comparison, numeric(7), restricted = M0))

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
    rownames(X[[2]]) <- apply(comparisons, 1, paste, collapse = "|")

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
    rownames(X[[3]]) <- gsub("preferential_acquisition_", "", apply(comparisons, 1, paste, collapse = "|"))

    # Comparison 4: Adult & Child vs. CHILDES (Preferential Acq. only) ----
    comparisons <- as.data.frame(rbind(
      c(full = 'preferential_acquisition_adult+preferential_acquisition_child+preferential_acquisition_childes', restricted = 'preferential_acquisition_adult+preferential_acquisition_child'),
      c(full = 'preferential_acquisition_adult+preferential_acquisition_child+preferential_acquisition_childes', restricted = 'preferential_acquisition_adult+preferential_acquisition_childes'),
      c(full = 'preferential_acquisition_adult+preferential_acquisition_child+preferential_acquisition_childes', restricted = 'preferential_acquisition_child+preferential_acquisition_childes')
    ))
    X[[4]] <- t(mapply(model_comp_helper, comparisons[["full"]], comparisons[["restricted"]], MoreArgs = list(M = c(M2n, M3n))))
    rownames(X[[4]]) <- gsub("preferential_acquisition_", "", apply(comparisons, 1, paste, collapse = "|"))

    model_comparisons <- do.call('rbind', X)
    model_comparisons <- cbind(
      model_comparisons,
      p_fdr  = p.adjust(as.vector(model_comparisons[, 'p']), method = 'fdr'),
      p_bonf = p.adjust(as.vector(model_comparisons[, 'p']), method = 'bonferroni'),
      p_holm = p.adjust(as.vector(model_comparisons[, 'p']), method = 'holm')
    )
    save(model_comparisons, file = sprintf('./network/model-comparisons-%03d.Rdata', repetition - 1))
}
parallel::stopCluster(cl)
