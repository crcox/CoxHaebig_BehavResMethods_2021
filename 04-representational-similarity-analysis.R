library('dplyr')
library('tidyr')

load_to_list <- function(files) {
    X <- new.env()
    lapply(files, load, envir = X)
    return(as.list(X))
}

trim_prefix <- function(x, prefix) {
    return(trimws(x, which = "left", whitespace = prefix))
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

generate_response_profiles <- function(d) {
    xtabs(~ RESPONSE + CUE, data = d)
}

correlate_binary_response_profiles <- function(x) {
    as.dist(cor(x > 0))
}

pivot_responses_wide <- function(d) {
    tidyr::pivot_wider(
        data = d,
        id = tidyselect::everything(),
        names_prefix = "R",
        names_from = "RESP_ID",
        values_from = "RESPONSE"
    )
}

pivot_responses_long <- function(d) {
    tidyr::pivot_longer(
        data = d,
        cols = c("R1", "R2", "R3"),
        names_prefix = "R",
        names_to = "RESP_ID",
        values_to = "RESPONSE"
    )
}

#' Split the data by chunk
#'
#' @param x a data frame
#' @param by a string indicating a column in \code{x} that will be used to
#'     split the data into chunks.
#' @return A list of lists. The top level will have as many elements as unique
#'     values in \code{by}, and the lower level will have as many elements as
#'     returned by \code{random_split}.
generate_splits <- function(d, by) {
    lapply(split(d, d[[by]]), random_split)
}

#' Create two sets of rows at random
#'
#' Every row in the input will be assigned to one of the sets, and the sets
#' will be as equal in size as possible.
#'
#' @param x a data frame
#' @return A list containing two data frames.
random_split <- function(x) {
    split(x, sample(nrow(x)) > floor(nrow(x) / 2))
}

#' Simulate null distribution of representational similarity
#'
#' Randomly split and recombine the data from adult and child conditions to
#' create two mixed datasets. Then compute response profiles and their
#' pairwise correlations. Finally, assess representational similarity between
#' the two correlation matrices.
#'
#' @param wide_assoc A list of two data frames, where each row contains a cue
#'     and three responses from one participant.
#' @return a correlation coefficient.
#'
#' @details
#' This function should be run many times to simulate a distribution of
#' correlation values under the null hypothesis that adult and child responses
#' are sampled from the same population.
simulate_null_cor <- function(wide_assoc) {
    splits  <- lapply(wide_assoc, generate_splits, by = "CUE")
    mixed <- list(
        do.call('rbind', lapply(c(splits$adult, splits$child), `[[`, 1)),
        do.call('rbind', lapply(c(splits$adult, splits$child), `[[`, 2))
    )
    repsim <- lapply(
        mixed,
        function(d) {
            correlate_binary_response_profiles(
                generate_response_profiles(
                    pivot_responses_long(d)
                )
            )
        }
    )
    return(cor(repsim[[1]], repsim[[2]], method = "spearman"))
}

# Load association data ----
associations <- load_to_list(c("./data/associations-adult.Rdata", "./data/associations-child.Rdata"))
names(associations) <- trim_prefix(names(associations), "associations_")
associations$adult$RESPONSE <- as.character(associations$adult$RESPONSE)
associations$child$RESPONSE <- as.character(associations$child$RESPONSE)
associations$adult$CUE <- as.character(associations$adult$CUE)
associations$child$CUE <- as.character(associations$child$CUE)
associations$adult$PP_ID <- as.character(associations$adult$PP_ID)
associations$child$PP_ID <- as.character(associations$child$PP_ID)

associations <- lapply(
    associations,
    function (d) {
        droplevels(subset(d, CUE !="so big!"))
    }
)

stopifnot(all.equal(
    levels(associations[['adult']]['CUE']),
    levels(associations[['child']]['CUE'])
))

# Check that at least 95 participants responded to all retained cues ----
pp_count <- vapply(
  associations,
  FUN = function(d) count_participants_by_cue(d$CUE, d$PP_ID),
  FUN.VALUE = numeric(length(unique(associations$adult$CUE)))
)
stopifnot(all(pmin(pp_count[, 1], pp_count[, 2]) > 95))

# Select 100 participants (300 responses) per cue ----
associations <- lapply(associations, sample_participants_by_cue, n = 100)

# Compute true representational similarity matrices ----
repsim <- lapply(
    associations,
    function(d) {
        correlate_binary_response_profiles(
            generate_response_profiles(d)
        )
    }
)
repsim_cor <- cor(repsim$adult, repsim$child, method = "spearman")

# Begin simulation ----
# Establish cluster
library('parallel')
cl <- parallel::makeCluster(parallel::detectCores() - 1)
parallel::clusterSetRNGStream(cl)
parallel::invisible(clusterEvalQ(cl, {
  library('tidyr')
}))
clusterExport(cl, c(
    "pivot_responses_wide",
    "pivot_responses_long",
    "generate_splits",
    "random_split",
    "simulate_null_cor",
    "generate_response_profiles",
    "correlate_binary_response_profiles"
))

# Run in parallel (this will take several minutes)
null_repsim_cor <- parSapply(
    cl,
    1:1000,
    FUN = function(i, data) simulate_null_cor(data),
    data = lapply(associations, pivot_responses_wide)
)

summary(null_repsim_cor)
print(mean(null_repsim_cor < repsim_cor))

save(null_repsim_cor, file = 'null-repsim-cor-adult-child.Rdata')
save(repsim_cor, file = 'repsim-cor-adult-child.Rdata')
