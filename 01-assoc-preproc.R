library(dplyr)
library(tidyr)

load_to_list <- function(...) {
    files <- rlang::list2(...)
    purrr::map(files, function(filename) {
        e <- new.env()
        nm <- load(filename, envir = e)
        e[[nm]]
    })
}

save_as <- function(..., filename) {
    dat <- rlang::list2(...)
    save(list = names(dat), file = filename, envir = list2env(dat))
}

inner_join_assoc_cdi <- function(d, id_tbl) {
    d <- dplyr::inner_join(d, id_tbl)
    d$CUE <- as.factor(d$CUE)
    return(subset(d, RESPONSE != ""))
}

count_participants_by_cue <- function(cues, participant_id) {
  X <- lapply(
    split(
      as.factor(participant_id),
      as.factor(cues),
      drop = TRUE
    ),
    FUN = droplevels)
  return(vapply(X = X, FUN = nlevels, FUN.VALUE = numeric(1)))
}

sample_participants_by_cue <- function(d, n) {
    tidyr::pivot_wider(
        data = d,
        id = tidyselect::everything(),
        names_prefix = "R",
        names_from = "RESP_ID",
        values_from = "RESPONSE"
    ) %>%
    dplyr::group_by(CUE) %>%
    dplyr::slice_sample(n = n) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_longer(
        cols = c("R1", "R2", "R3"),
        names_to = "RESP_ID",
        names_prefix = "R",
        values_to = "RESPONSE"
    )
}

# Load metadata and associations ----
load("./data/cdi-metadata-preproc.Rdata")
associations <- load_to_list(
    adult = "./data/associations-adult.Rdata",
    child = "./data/associations-child.Rdata"
)

# Drop cues to match metadata-preproc ----
# Also, add CDI item id and lemma columns to word associations, and remove any
# empty responses.
associations <- lapply(
    associations,
    FUN = inner_join_assoc_cdi,
    id_tbl = dplyr::select(cdi_metadata_preproc, CUE = cue_CoxHae, num_item_id, lemma)
)

# Check that at least 95 participants responded to all retained cues ----
pp_count <- vapply(
  associations,
  FUN = function(d) count_participants_by_cue(d$CUE, d$PP_ID),
  FUN.VALUE = numeric(nlevels(associations$adult$CUE))
)
summary(pp_count)
stopifnot(all(pmin(pp_count[, 1], pp_count[, 2]) > 95))


# Select 100 participants (300 responses) per cue ----
# For cues with fewer than 100 participants, include all.
# For cues with more than 100 participants, exclude some to get to 100.
associations <- lapply(
    associations,
    FUN = sample_participants_by_cue,
    n = 100
)

# Save processed associations ----
names(associations) <- add_prefix(names(associations), prefix = "associations_")
names(associations) <- add_suffix(names(associations), suffix = "_preproc")

save_as(
    associations_adult_preproc = associations$adult,
    filename = "./data/associations-adult-preproc.Rdata"
)
save_as(
    associations_child_preproc = associations$child,
    filename = "./data/associations-child-preproc.Rdata"
)
