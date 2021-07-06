library('dplyr')
library('tidyr')
library('netgrowr')

load_to_list <- function(files) {
    X <- new.env()
    lapply(files, load, envir = X)
    return(as.list(X))
}

trim_prefix <- function(x, prefix) {
    return(trimws(x, which = "left", whitespace = prefix))
}

add_prefix <- function(x, prefix) {
    str_prepend <- function(prefix, ...) paste(prefix, ..., sep = "")
    return(vapply(x, FUN = str_prepend, FUN.VALUE = character(1), prefix = prefix, USE.NAMES = FALSE))
}

# Load id key ----
load('./data/cdi-metadata-preproc.Rdata')
assocnet <- load_to_list(c(
  "./network/assocnet-adult-preproc.Rdata",
  "./network/assocnet-child-preproc.Rdata",
  "./network/assocnet-childes-preproc.Rdata"
))
names(assocnet) <- trim_prefix(names(assocnet), prefix = "assocnet_")

growthvalues <- lapply(
    assocnet,
    netgrowr::growth_values,
    aoa_tbl = cdi_metadata_preproc[c('lemma', 'aoa_produces')],
    growth_models = c("preferential_attachment", "lure_of_the_associates", "preferential_acquisition")
)

# Save growth values ----
names(growthvalues) <- add_prefix(names(growthvalues), "growthvalues_")
save(growthvalues_adult_preproc, file = './network/growthvalues-adult-preproc.Rdata', envir = list2env(growthvalues))
save(growthvalues_child_preproc, file = './network/growthvalues-child-preproc.Rdata', envir = list2env(growthvalues))
save(growthvalues_childes_preproc, file = './network/growthvalues-childes-preproc.Rdata', envir = list2env(growthvalues))

# Save wide-form data for lexical growth modeling ----
# Specify and incorporate phonological baseline variables
phono_baseline <- select(
    cdi_metadata_preproc,
    word = lemma,
    num_item_id,
    nphon,
    CHILDES_Freq,
    BiphonProb.avg,
    PNDC.avg
)

growthvalues$growthvalues_adult_preproc$network <- "adult"
growthvalues$growthvalues_child_preproc$network <- "child"
growthvalues$growthvalues_childes_preproc$network <- "childes"

modelvars <- tidyr::pivot_wider(
    do.call('rbind', growthvalues),
    id_cols = c('word', 'month', 'aoa', 'learned'),
    names_from = c(model, network),
    values_from = 'value'
) %>%
    left_join(phono_baseline, by = "word")

save(modelvars, file = "./network/modelvars.Rdata")
