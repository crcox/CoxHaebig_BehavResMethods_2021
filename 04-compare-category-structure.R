library('igraph')
library('boot')

trim_prefix <- function(x, prefix) {
    return(trimws(x, which = "left", whitespace = prefix))
}

trim_suffix <- function(x, suffix) {
    return(trimws(x, which = "right", whitespace = suffix))
}

combn_and_paste <- function(x, m, sep = "-") {
  return(apply(combn(x, m), MARGIN = 2, paste, collapse = sep))
}

load_to_list <- function(files) {
    X <- new.env()
    lapply(files, load, envir = X)
    return(as.list(X))
}

get_bca_ci_values <- function(x, data) {
    boot::boot.ci(data, index = x, type = "bca")[["bca"]][4:5]
}

summary.boot <- function(b, index = 1:length(b$t0)) {
    out <- cbind(
        b$t0[index],
        apply(b$t[, index], 2, sd),
        t(vapply(index, get_bca_ci_values, numeric(2), data = b))
    )
    colnames(out) <- c("mean", "SE", "ci.l", "ci.u")
    cat('\n')
    print(out)
    cat('\n')
    cat("Confidence intervals are 95% computed via", b$sim, "bootstrap and BCa.\n")
    cat("Based on", b$R, "bootstrap replicates.", '\n')
    cat("Calculations and intervals on original scale.", '\n')
    cat('\n')
    invisible(out)
}

# Load processed word associations ----
load('./data/cdi-metadata-preproc.Rdata')
assocnet <- load_to_list(c(
  "./network/assocnet-adult-preproc.Rdata",
  "./network/assocnet-child-preproc.Rdata",
  "./network/assocnet-childes-preproc.Rdata"
))
names(assocnet) <- trim_suffix(trim_prefix(names(assocnet), "assocnet_"), "_preproc")

# Set diagonal to zero (no loops in network)
assocnet <- lapply(assocnet, function(m) {diag(m) <- FALSE; return(m)})
assocnet_igraph <- lapply(assocnet, igraph::graph_from_adjacency_matrix, mode = 'directed')
assocnet_distances <- lapply(assocnet_igraph, igraph::distances)
replace_inf_with_max <- function(d) {
  m <- max(is.finite(d))
  d[is.infinite(d)] <- m
  return(d)
}
assocnet_distances <- lapply(assocnet_distances, replace_inf_with_max)
assocnet_distances <- lapply(assocnet_distances, as.dist)

x <- merge(
  data.frame(lemma = labels(assocnet_distances[[1]])),
  cdi_metadata_preproc[c('lemma', 'num_item_id', 'category')],
  sort = FALSE
)

stopifnot(all.equal(x$lemma, labels(assocnet_distances[[1]])))
stopifnot(all.equal(x$lemma, labels(assocnet_distances[[2]])))
stopifnot(all.equal(x$lemma, labels(assocnet_distances[[3]])))

within_cat <- outer(x$category, x$category, FUN = "==")
within_cat_vec <- within_cat[lower.tri(within_cat)]

bootfun <- function(distances, ix, within_cat) {
  # ix is ignored
  d <- lapply(distances, as.vector)
  between_cat <- !within_cat
  w <- vapply(d, function(x, z) return(mean(sample(x[z], replace = TRUE))), numeric(1), z = within_cat)
  b <- vapply(d, function(x, z) return(mean(sample(x[z], replace = TRUE))), numeric(1), z = between_cat)
  x <- w / b
  out <- c(x, x[1] - x[2], x[1] - x[3], x[2] - x[3])
  names(out) <- c(names(d), combn_and_paste(names(d), 2))
  return(out)
}

BB <- boot::boot(assocnet_distances[c('adult', 'child', 'childes')], bootfun, 10000, within_cat = within_cat_vec, parallel = "multicore", ncpus = 30)

colnames(BB$t) <- names(BB$t0)
boot::boot.ci(BB, index = "adult", type = "bca")
boot::boot.ci(BB, index = "child", type = "bca")
boot::boot.ci(BB, index = "childes", type = "bca")
boot::boot.ci(BB, index = "adult-child", type = "bca")
boot::boot.ci(BB, index = "adult-childes", type = "bca")

boot::boot.ci(BB, index = "child-childes", type = "bca")

comparisons <- c("adult-child", "adult-childes", "child-childes")
summary(BB, index = comparisons)
