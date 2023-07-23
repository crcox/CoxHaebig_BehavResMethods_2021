#' Merge Kuperman AoA norms into dataframe of cued responses
#'
#' @param R A dataframe of unique responses in a column named RESPONSE
#'   associated with each cue in a column named CUE.
#' @param KupAoA The table of age of acquisition (AoA) ratings (and other
#'   psycholinguistic variables) published by Kuperman et al 2012.
#' @return A dataframe with the same number of rows as R, with additional
#'   columns including psycholinguistic norms.
#' @details Table of AoA ratings can be found at
#'   http://crr.ugent.be/papers/AoA_ratings_Kuperman_et_al_BRM.zip and the
#'   Behavior Research Methods article can be found
#'   http://crr.ugent.be/papers/Kuperman%20et%20al%20AoA%20ratings.pdf
#'   
#'   Taken from page 9 of linked manuscript by Kuperman et al. (2012):
#'   "We used the same instructions as for the collection of the Bristol norms
#'   (StadthagenGonzalez & Davis, 2006). Participants were asked for each word
#'   to enter the age (in years) at which they thought they had learned the
#'   word. It was specified that by learning a word, "we mean the age at which
#'   you would have understood that word if somebody had used it in front of
#'   you, EVEN IF YOU DID NOT use, read or write it at the time". Unlike many
#'   other studies, we did not ask participants to use a 7-point Likert rating
#'   scale, because this artificially restricts the response range and is also
#'   more difficult for participants to use (see Ghyselinck, De Moor, &
#'   Brysbaert, 2000, for a comparison of both methods; also see Figure 3
#'   below). When participants did not know a word, they were asked to enter the
#'   letter x. This prevented us from collecting wrong AoA ratings and also
#'   provided us with an estimate of how familiar responders were with the
#'   words. A complete list of 362 words (300 test words, 10 calibrator words,
#'   and 52 control words) took some 20 minutes to complete. Participants were
#'   paid half a dollar cent per rated word (i.e., $1.81 for a validly completed
#'   list)."
merge_kuperman <- function(R, KupAoA) {
  z <- KupAoA$Word != KupAoA$Alternative.spelling
  R <- left_join(R, select(KupAoA, RESPONSE = Word, Freq_pm, Dom_PoS_SUBTLEX, Nletters, Nphon, Nsyll, Lemma_highest_PoS, AoA_Kup, AoA_Kup_lem), by = "RESPONSE")
  R <- left_join(R, select(subset(KupAoA, z), RESPONSE = Alternative.spelling, Freq_pm, Dom_PoS_SUBTLEX, Nletters, Nphon, Nsyll, Lemma_highest_PoS, AoA_Kup, AoA_Kup_lem), by = "RESPONSE")
  # The second join will result in duplicated columns. The following
  # consolidates the data
  R$Freq_pm.x <- ifelse(is.na(R$Freq_pm.x),R$Freq_pm.y,R$Freq_pm.x)
  R$Dom_PoS_SUBTLEX.x <- ifelse(is.na(R$Dom_PoS_SUBTLEX.x),R$Dom_PoS_SUBTLEX.y,R$Dom_PoS_SUBTLEX.x)
  R$Nletters.x <- ifelse(is.na(R$Nletters.x),R$Nletters.y,R$Nletters.x)
  R$Nphon.x <- ifelse(is.na(R$Nphon.x),R$Nphon.y,R$Nphon.x)
  R$Nsyll.x <- ifelse(is.na(R$Nsyll.x),R$Nsyll.y,R$Nsyll.x)
  R$Lemma_highest_PoS.x <- ifelse(is.na(R$Lemma_highest_PoS.x),R$Lemma_highest_PoS.y,R$Lemma_highest_PoS.x)
  R$AoA_Kup.x <- ifelse(is.na(R$AoA_Kup.x),R$AoA_Kup.y,R$AoA_Kup.x)
  R$AoA_Kup_lem.x <- ifelse(is.na(R$AoA_Kup_lem.x),R$AoA_Kup_lem.y,R$AoA_Kup_lem.x)
  # Finally, select and rename columns
  R <- R %>%
    select(-ends_with(".y")) %>%
    rename(
      Freq_pm = Freq_pm.x,
      Dom_PoS_SUBTLEX = Dom_PoS_SUBTLEX.x,
      Nletters = Nletters.x,
      Nphon = Nphon.x,
      Nsyll = Nsyll.x,
      Lemma_highest_PoS = Lemma_highest_PoS.x,
      AoA_Kup = AoA_Kup.x,
      AoA_Kup_lem = AoA_Kup_lem.x
    )
  return(R)
}
