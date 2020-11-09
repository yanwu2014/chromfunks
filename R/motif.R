#' Computes motif enrichment with a hypergeometric distribution
#'
#' @param meta.features Dataframe with gc content and peak names
#' @param features Vector peak names to test for enrichment. All peaks must be in meta.features
#' @param motif.all Matrix of peaks by motifs
#' @param n.background Number of background peaks to sample
#' @param verbose Print progress
#'
#' @return Dataframe with motif enrichment results
#' @import Signac
#' @export
#'
MotifEnrich <- function(meta.features, features, motif.all,
                        n.background = 40000, verbose = TRUE) {
  rownames(meta.features) <- meta.features$Peak
  background <- Signac::MatchRegionStats(meta.feature = meta.features, regions = features,
                                         n = n.background, verbose = verbose)
  if (verbose) {
    message("Testing motif enrichment in ", length(x = features),
            " regions")
  }
  query.motifs <- motif.all[features, ]
  background.motifs <- motif.all[background, ]

  query.counts <- colSums(x = query.motifs)
  background.counts <- colSums(x = background.motifs)
  percent.observed <- query.counts/length(x = features) * 100
  percent.background <- background.counts/length(x = background) * 100
  fold.enrichment <- percent.observed/percent.background
  p.list <- vector(mode = "numeric")
  for (i in seq_along(along.with = query.counts)) {
    p.list[[i]] <- phyper(q = query.counts[[i]] - 1, m = background.counts[[i]],
                          n = nrow(x = background.motifs) - background.counts[[i]],
                          k = length(x = features), lower.tail = FALSE)
  }
  # qv <- qvalue::qvalue(p = p.list)
  qv <- p.adjust(p.list, method = "BH")
  results <- data.frame(motif = names(x = query.counts), observed = query.counts,
                        background = background.counts, percent.observed = percent.observed,
                        percent.background = percent.background, fold.enrichment = fold.enrichment,
                        pvalue = p.list, FDR = qv,
                        stringsAsFactors = FALSE)
  if (nrow(x = results) == 0) {
    return(results)
  }
  else {
    return(results[order(results[, 7], -results[, 6]), ])
  }
}
