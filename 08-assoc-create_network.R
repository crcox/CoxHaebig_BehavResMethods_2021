library('dplyr')
library('tidyr')
library('netbuildr')
# netbuildr can be obtained from https://github.com/crcox/netbuilder
# install with:
#     devtools::install_github("crcox/netbuilder")

trim_prefix <- function(x, prefix) {
    return(trimws(x, which = "left", whitespace = prefix))
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
associations <- load_to_list(c(
  "./data/associations-adult-preproc.Rdata",
  "./data/associations-child-preproc.Rdata"
))
names(associations) <- trim_prefix(names(associations), "associations_")

# Create associative networks ----
assocnet <- lapply(
  associations,
  FUN = function(d) netbuildr::create_unweighted_network(d[['lemma']], d[['RESPONSE']])
)
assocnet <- lapply(
  assocnet,
  FUN = function(x, ix) x[ix, ix],
  ix = as.character(cdi_metadata_preproc$lemma)
)

# Save association networks ----
names(assocnet) <- add_prefix(names(assocnet), prefix = "assocnet_")

if (!dir.exists("./network")) dir.create("./network")
save(assocnet_adult_preproc, file = './network/assocnet-adult-preproc.Rdata', envir = list2env(assocnet))
save(assocnet_child_preproc, file = './network/assocnet-child-preproc.Rdata', envir = list2env(assocnet))
