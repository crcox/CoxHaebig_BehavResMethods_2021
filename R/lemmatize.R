lemmatize <- function(x, dict) {
    dict[[1]] <- as.character(dict[[1]])
    dict[[2]] <- as.character(dict[[2]])
    x <- sub("'s$", "", tolower(trimws(as.character(x))))
    x <- gsub("([a-z])\\1{2,}", replacement = "\\1\\1", x = x)
    lemma <- lemmatize_words(x, dictionary = dict)
    return(lemma)
}

# lemmatize_words and check_dictionary are taken from the textstem package by Tyler Rinker.
# Rinker, T. W. (2018). textstem: Tools for stemming and lemmatizing text
#     version 0.1.4. Buffalo, New York. http://github.com/trinker/textstem
lemmatize_words <- function(x, dictionary) {
    locs <- match(tolower(x), dictionary[[1]])
    x[!is.na(locs)] <- dictionary[locs[!is.na(locs)], ][[2]]
    return(x)
}
check_dictionary <- function(x) {
    if (anyDuplicated(x[[1]]) > 0) {
        stop("Duplicate tokens found in column one of the lemma dictionary.")
    }
}
