tokenize_phrase <- function(.data, phrase, phrase_token = NULL, sort_tokens = TRUE) {
    # Sort by transcript, utterance, and token order
    if (sort_tokens) {
        ix <- with(.data, order(transcript_id, utterance_id, token_order))
        .data <- .data[ix, ]
    }

    # Find instances of the phrase as a sequence of tokens.
    n <- length(phrase)
    ix <- which(tolower(.data$gloss) == tolower(phrase[1]))
    z <- !logical(length(ix))
    for (i in 2:n) {
        z <- z & tolower(.data$gloss[ix + (i - 1)]) == tolower(phrase[i])
    }
    ix <- ix[z]
    ixx <- matrix(ix, nrow = n, ncol = length(ix), byrow = TRUE) + seq(0, n - 1)

    # Check that instances of the phrase do not span multiple utterances.
    uid <- matrix(.data$utterance_id[c(ixx)], nrow = n, ncol = ncol(ixx))
    z <- apply(uid, 2, function(x) all(x[1] == x[-1]))
    ix <- ix[z]
    ixx <- ixx[ , z]
    if (is.null(phrase_token)) {
        phrase_token <- paste(phrase, collapse = '+')
    }
    cat("Found", sum(z), "instances of phrase", phrase_token, "\n")
    if (sum(z) > 0) {
        stem_token <- character(ncol(ixx))
        pos_token <- character(ncol(ixx))
        for (i in 1:ncol(ixx)) {
            stem_token[i] <- paste(.data$stem[ixx[, i]], collapse = '+')
            pos_token[i] <- paste(.data$part_of_speech[ixx[, i]], collapse = '+')
        }
        .data$gloss[ix] <- phrase_token
        .data$stem[ix] <- stem_token
        .data$part_of_speech[ix] <- pos_token
        ixd <- c(ixx[-1, ])
        .data <- .data[-ixd, ]
    }
    return(.data)
}
