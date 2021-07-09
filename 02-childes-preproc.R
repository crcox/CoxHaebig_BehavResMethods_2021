library('dplyr')
source('./R/lemmatize.R')
source('./R/recode_nonwords.R')
source('./R/tokenize_phrase.R')

# Get CHILDES tokens ----
# Local version acquired 13 June 2020 using:
#    roles <- c("Adult", "Father", "Mother", "Aunt", "Uncle", "Grandmother",
#              "Grandfather", "Teacher", "Babysitter", "Nurse", "Doctor",
#              "Clinician", "Therapist")
#    childes_tokens_preproc <- childesr::get_tokens(collection = "Eng-NA", role = roles, token='*')
#    save(childes_tokens_preproc, file = './data/childes_tokens_preproc_Eng-NA_2020June13_raw.Rdata')

# Preprocessing ...
load('./data/cdi-metadata-preproc.Rdata')
load('./data/childes-tokens.Rdata')
lemma_dict_childes <- read.csv('./data/lemma-dict-childes.csv')

childes_tokens_preproc <- subset(childes_tokens, childes_tokens[['target_child_age']] <= 60)
childes_tokens_preproc <- tokenize_phrase(childes_tokens_preproc, c("so", "big"),       sort_tokens = TRUE)
childes_tokens_preproc <- tokenize_phrase(childes_tokens_preproc, c("go", "potty"),     sort_tokens = FALSE)
childes_tokens_preproc <- tokenize_phrase(childes_tokens_preproc, c("night", "night"),  sort_tokens = FALSE)
childes_tokens_preproc <- tokenize_phrase(childes_tokens_preproc, c("thank", "you"),    sort_tokens = FALSE)
childes_tokens_preproc <- tokenize_phrase(childes_tokens_preproc, c("all", "gone"),     sort_tokens = FALSE)
childes_tokens_preproc <- tokenize_phrase(childes_tokens_preproc, c("next", "to"),      sort_tokens = FALSE)
childes_tokens_preproc <- tokenize_phrase(childes_tokens_preproc, c("on", "top", "of"), sort_tokens = FALSE)
childes_tokens_preproc <- tokenize_phrase(childes_tokens_preproc, c("a", "lot"),        sort_tokens = FALSE)
childes_tokens_preproc <- tokenize_phrase(childes_tokens_preproc, c("going", "to"),     sort_tokens = FALSE)
childes_tokens_preproc <- tokenize_phrase(childes_tokens_preproc, c("got", "to"),       sort_tokens = FALSE)
childes_tokens_preproc <- tokenize_phrase(childes_tokens_preproc, c("have", "to"),      sort_tokens = FALSE)
childes_tokens_preproc <- tokenize_phrase(childes_tokens_preproc, c("let", "me"),       sort_tokens = FALSE)
childes_tokens_preproc <- tokenize_phrase(childes_tokens_preproc, c("want", "to"),      sort_tokens = FALSE)
childes_tokens_preproc <- tokenize_phrase(childes_tokens_preproc, c("did", "ya"), phrase_token = "did+you", sort_tokens = FALSE)
childes_tokens_preproc <- tokenize_phrase(childes_tokens_preproc, c("did", "you"),      sort_tokens = FALSE)
childes_tokens_preproc <- tokenize_phrase(childes_tokens_preproc, c("need", "to"),      sort_tokens = FALSE)
childes_tokens_preproc <- tokenize_phrase(childes_tokens_preproc, c("try", "to"),       sort_tokens = FALSE)
childes_tokens_preproc <- tokenize_phrase(childes_tokens_preproc, c("play", "pen"),     sort_tokens = FALSE)
childes_tokens_preproc$lemma <- lemmatize(childes_tokens_preproc$gloss, dict = lemma_dict_childes)
childes_tokens_preproc$lemma <- recode_nonwords(childes_tokens_preproc$lemma, childes_tokens_preproc$part_of_speech)
childes_tokens_preproc <- subset(childes_tokens_preproc, !is.na(childes_tokens_preproc$lemma))

# Extract information about target children ----
child_id_fields <- c(
  "target_child_id",
  "target_child_age",
  "target_child_name",
  "target_child_sex"
)
childes_children <- subset(
  childes_tokens_preproc,
  !duplicated(childes_tokens_preproc[child_id_fields]),
  select = child_id_fields
)
cat("Summary of target child age in CHILDES selection\n")
cat("------------------------------------------------\n")
print(summary(childes_children$target_child_age))
stem(childes_children$target_child_age)

# Confirm that all CDI word lemmas appear among CHILDES tokens ----
all_lemmas_appear_in_childes <- all(cdi_metadata_preproc$lemma %in% childes_tokens_preproc$lemma)
stopifnot(all_lemmas_appear_in_childes)
cat("All study lemmas occur in CHILDES.")

# Arrange tokens by utterances within transcripts ----
childes_tokens_preproc <- dplyr::arrange(
  childes_tokens_preproc,
  transcript_id,
  utterance_id,
  token_order
)

# Save processed CHILDES tokens ----
save(childes_tokens_preproc, file = './data/childes-tokens-preproc.Rdata')
save(childes_children, file = './data/childes-children.Rdata')
