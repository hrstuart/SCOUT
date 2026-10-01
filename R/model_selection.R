# =============================================================================
# Hybrid selection: AIC/AICc for BM1 vs OU1, penalized-fit LRT for regime number
# =============================================================================
# Why not test everything against BM1?
#  * SCOUT puts a Gamma(2,1) prior on alpha (weight lambda2*N) for every OU fit and none on BM1,
#    so alpha-hat never approaches 0 and OU1 is NOT nested in BM1 in practice: on true-BM1 genes
#    ll(OU1) < ll(BM1) in 84-100% of fits. A chi^2 LRT between them is invalid -> AIC/AICc only.
#  * OUM-type fits (any partition) vs OU1 share the alpha prior and thetas are unpenalized, so that
#    pair stays nested; its chi^2 LRT held 4-5% size at alpha=3 in the 260922 regime-perturb sims.
#  * The multi-regime model tested is the IC-best of n_multi candidates, so its p-value is
#    Bonferroni-corrected by n_multi (uncorrected p<0.05 left 11-18% of true OU1 called multi).
#
# Procedure per gene:
#   1. winner = argmin IC over all candidates (ic = 'AICc' or 'AIC', recomputed from ll/k/n).
#   2. winner is BM1 or OU1 -> keep it.
#   3. winner is multi-regime -> LRT winner vs OU1 (mvMORPH convention: larger model is the alt,
#      df = k difference, negative LR clipped to 0). p * n_multi < alpha -> keep winner;
#      otherwise on_fail = 'fallback' -> IC-best of BM1/OU1, or 'ambiguous' -> labelled ambiguous.
# lrt_status: 'not tested' (BM1 best), 'not eligible' (OU1 best; no valid LRT vs BM1),
#             'tested: significant' / 'tested: not significant' (multi-regime winner vs OU1).
# The statistics use SCOUT's unpenalized ll_total at penalized estimates: report as non-standard.

#' Hybrid IC + LRT model selection for one gene
#'
#' AIC/AICc decides between BM1 and OU1 (no valid LRT between them under SCOUT's alpha prior); a
#' multi-regime winner is confirmed by an LRT against OU1, Bonferroni-corrected by the number of
#' multi-regime candidates. See the header of \code{R/model_selection.R} for the rationale.
#'
#' @param gene_df one gene's candidate fits, one row per regime (last EM iteration only).
#' @param ic information criterion used to pick the winner: \code{'AICc'} or \code{'AIC'}.
#' @param on_fail when the multi-regime winner is not significant: \code{'fallback'} calls the
#'   IC-best of BM1/OU1, \code{'ambiguous'} labels the gene ambiguous.
#' @param alpha significance level.
#' @param correction \code{'bonferroni'} (by n_multi) or \code{'none'}.
#' @param bm_regime,ou_regime names of the BM1 and OU1 rows.
#' @param regime_col,loglik_col,npar_col,ntips_col column names.
#' @return A one-row tibble: ic_winner, best_simple, delta_ic_simple, n_multi, LRT fields,
#'   \code{call}, \code{lrt_status} and \code{flag}.
#' @export
select_model_hybrid <- function(gene_df,
                                ic          = c("AICc", "AIC"),
                                on_fail     = c("fallback", "ambiguous"),
                                alpha       = 0.05,
                                correction  = c("bonferroni", "none"),
                                bm_regime   = "BM1",
                                ou_regime   = "OU1",
                                regime_col  = "regime",
                                loglik_col  = "ll_total",
                                npar_col    = "param.count",
                                ntips_col   = "ntips") {
  ic         <- match.arg(ic)
  on_fail    <- match.arg(on_fail)
  correction <- match.arg(correction)

  reg <- as.character(gene_df[[regime_col]])
  if (anyDuplicated(reg)) {
    stop("More than one row per regime (e.g. several EM iterations); keep the last iter per regime first.")
  }
  for (r in c(bm_regime, ou_regime)) {
    if (!r %in% reg) stop(sprintf("No '%s' row for this gene.", r))
  }

  ll <- gene_df[[loglik_col]]
  k  <- gene_df[[npar_col]]
  n  <- gene_df[[ntips_col]][1]
  crit <- if (ic == "AIC") -2 * ll + 2 * k else -2 * ll + 2 * k * n / (n - k - 1)
  names(crit) <- names(ll) <- names(k) <- reg

  is_multi <- !reg %in% c(bm_regime, ou_regime)
  n_multi  <- sum(is_multi)
  winner   <- reg[which.min(crit)]
  simple   <- if (crit[bm_regime] <= crit[ou_regime]) bm_regime else ou_regime

  out <- dplyr::tibble(
    ic                = ic,
    ic_winner         = winner,
    best_simple       = simple,
    delta_ic_simple   = unname(crit[simple] - min(crit)),   # >0: winner beats best of BM1/OU1
    n_multi           = n_multi,
    lrt_alt           = NA_character_,
    lrt_null          = NA_character_,
    lrt_statistic     = NA_real_,
    lrt_statistic_raw = NA_real_,
    lrt_df            = NA_integer_,
    lrt_p_value       = NA_real_,
    lrt_p_corrected   = NA_real_,
    call              = winner,
    lrt_status        = NA_character_,
    flag              = ""
  )

  # BM1 / OU1 winners: the IC decision stands, with the reason no LRT was run
  if (winner == bm_regime) {
    out$lrt_status <- "not tested"
    out$flag       <- sprintf("%s best by %s; simplest candidate, no null to test against", bm_regime, ic)
    return(out)
  }
  if (winner == ou_regime) {
    out$lrt_status <- "not eligible"
    out$flag       <- sprintf("%s best by %s; LRT vs %s not valid (not nested under SCOUT's alpha penalty)",
                              ou_regime, ic, bm_regime)
    return(out)
  }

  raw  <- 2 * (ll[[winner]] - ll[[ou_regime]])
  stat <- max(raw, 0)
  df   <- as.integer(k[[winner]] - k[[ou_regime]])
  if (df <= 0) stop(sprintf("'%s' does not have more parameters than '%s'.", winner, ou_regime))
  p    <- stats::pchisq(stat, df = df, lower.tail = FALSE)
  pc   <- if (correction == "bonferroni") min(1, p * n_multi) else p

  out$lrt_alt           <- winner
  out$lrt_null          <- ou_regime
  out$lrt_statistic     <- stat
  out$lrt_statistic_raw <- raw
  out$lrt_df            <- df
  out$lrt_p_value       <- p
  out$lrt_p_corrected   <- pc
  out$lrt_status        <- if (pc < alpha) "tested: significant" else "tested: not significant"
  flags <- character(0)
  if (raw < 0) flags <- c(flags, "negative LR clipped to 0")
  if (pc >= alpha) {
    out$call <- if (on_fail == "fallback") simple else "ambiguous"
    flags    <- c(flags, sprintf("%s not significant vs %s", winner, ou_regime))
  }
  out$flag <- paste(flags, collapse = "; ")
  out
}

#' Hybrid IC + LRT model selection across genes
#'
#' Applies \code{select_model_hybrid()} to every gene of a SCOUT
#' \code{*_all_genes_full_history.csv} (or the \code{annotated} history inside \code{SCOUT()}).
#' Rows are reduced to the last \code{iter_col} per (group_vars, regime) first, as
#' \code{annotate_history()} does. \code{fdr_method} adds a genome-wide adjustment of
#' \code{lrt_p_corrected} (NA for BM1/OU1 winners) over all rows passed in -- run once per dataset
#' if FDR should be per dataset. With \code{call_on_fdr = TRUE}, multi-regime calls are re-decided
#' on that adjusted p instead of the per-gene one.
#'
#' @param df a full-history data.frame, or one or more paths to full-history CSVs (row-bound).
#' @param group_vars columns identifying one gene's candidate set.
#' @param iter_col EM iteration column; only the last iteration per regime is kept.
#' @param regime_col regime column.
#' @param fdr_method \code{p.adjust} method, or NULL for none.
#' @param call_on_fdr re-decide multi-regime calls on the FDR-adjusted p.
#' @param ... passed to \code{select_model_hybrid()} (e.g. \code{ic}, \code{alpha}, \code{on_fail}).
#' @return One row per gene with the hybrid call, \code{decided_on} and (with FDR) \code{lrt_p_fdr}.
#' @importFrom rlang .data
#' @export
run_hybrid_pipeline <- function(df, group_vars = c("dataset", "gene_name"), iter_col = "iter",
                                regime_col = "regime", fdr_method = "BH", call_on_fdr = FALSE, ...) {
  if (is.character(df)) df <- dplyr::bind_rows(lapply(df, utils::read.csv))
  df <- df %>% dplyr::ungroup() %>% dplyr::select(-dplyr::any_of(c("X", "...1")))  # write.csv row-name column
  missing <- setdiff(c(group_vars, regime_col), names(df))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
  if (iter_col %in% names(df)) {
    df <- df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c(group_vars, regime_col)))) %>%
      dplyr::slice_max(.data[[iter_col]], n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
  }

  args <- list(...)
  res  <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::reframe(select_model_hybrid(dplyr::pick(dplyr::everything()), regime_col = regime_col, ...))

  res$decided_on <- ifelse(is.na(res$lrt_p_corrected), "IC", "per-gene p")
  if (is.null(fdr_method)) return(res)
  res$lrt_p_fdr <- stats::p.adjust(res$lrt_p_corrected, method = fdr_method)

  if (call_on_fdr) {
    # flag keeps the per-gene verdict; decided_on marks calls re-made on the FDR-adjusted p
    a       <- if (is.null(args$alpha)) 0.05 else args$alpha
    on_fail <- if (is.null(args$on_fail)) "fallback" else match.arg(args$on_fail, c("fallback", "ambiguous"))
    tested  <- !is.na(res$lrt_p_fdr)
    res$call[tested] <- ifelse(res$lrt_p_fdr[tested] < a, res$ic_winner[tested],
                               if (on_fail == "fallback") res$best_simple[tested] else "ambiguous")
    res$lrt_status[tested] <- ifelse(res$lrt_p_fdr[tested] < a, "tested: significant", "tested: not significant")
    res$decided_on[tested] <- paste0("FDR p (", fdr_method, ")")
  }
  res
}
