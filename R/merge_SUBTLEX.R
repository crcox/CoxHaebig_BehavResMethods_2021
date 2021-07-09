#' Merge SUBTLEX frequency norms into dataframe of cued responses
#'
#' @param R A dataframe of unique responses in a column named RESPONSE
#'   associated with each cue in a column named CUE.
#' @param SUBTLEX The table of frequency norms (and other
#'   psycholinguistic variables) published by Brysbaert et al 2012.
#' @return A dataframe with the same number of rows as R, with additional
#'   columns including psycholinguistic norms.
#' @details Table of frequency norms can be found at
#'   https://www.ugent.be/pp/experimentele-psychologie/en/research/documents/subtlexus/subtlexus4.zip
#'   and the Behavior Research Methods article can be found
#'   https://www.ugent.be/pp/experimentele-psychologie/en/research/documents/subtlexus/brysbaertnew.pdf
#'   
#'   SUBTLEX word frequency data are derived from television and film subtitles.
merge_SUBTLEX <- function(R, SUBTLEX) {
  R <- left_join(
    x = R,
    y = select(SUBTLEX, RESPONSE = Word, SUBTLWF, Lg10WF, SUBTLCD, Lg10CD),
    by = 'RESPONSE'
  )
  return(R)
}