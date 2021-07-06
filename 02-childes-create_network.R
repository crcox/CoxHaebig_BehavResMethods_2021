library('dplyr')
library('netbuildr')

# Load CHILDES tokens and CDI metadata ----
load('./data/childes-tokens-preproc.Rdata')
load('./data/cdi-metadata-preproc.Rdata')

# Compute cooccurrence matrix with window size 5 ----
# N.B Riordan and Jones (2007) Cog. Sci. used windows size 10
k <- 5L # Hills et al. (2010) J. Mem. Lang.

# Ensure tokens are ordered and blocked by ordered utterances within
# transcripts
childes_tokens_preproc <- dplyr::arrange(
  childes_tokens_preproc,
  transcript_id,
  utterance_id,
  token_order
)

transcripts <- split(childes_tokens_preproc$lemma, childes_tokens_preproc$transcript_id)
cooccurrences <- netbuildr::create_cooccurrence_matrix(
  tokens = transcripts,
  window_size = k,
  types = cdi_metadata_preproc$lemma
)

# Save network variants ----
# UNWEIGHTED
assocnet_childes_preproc <- cooccurrences > 0

if (!dir.exists("./network")) dir.create("./network")
save(assocnet_childes_preproc, file = "./network/assocnet-childes-preproc.Rdata")
