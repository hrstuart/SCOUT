# lrt_utils.R -- likelihood-ratio-test model selection on SCOUT fit histories: AIC-best model with
# a confirmatory LRT vs BM1 (run_lrt_pipeline), or hybrid AIC/AICc selection with a
# Bonferroni-corrected multi-regime-vs-OU1 LRT (run_hybrid_pipeline).
#
# Usage:
#   source('lrt_utils.R')
#   sel <- run_hybrid_pipeline(full_history_df)      # one row per (dataset, gene_name)

library(dplyr)
library(tibble)

# --- LRT helper (same convention as before: expects TRUE log-likelihood) ---
lrt_stat <- function(ll_null, ll_alt, k_null, k_alt) {
  df_diff <- k_alt - k_null
  if (df_diff <= 0) {
    stop("Alt model must have strictly more free parameters than the null model.")
  }
  stat <- 2 * (ll_alt - ll_null)
  if (stat < 0) {
    warning("Negative LRT statistic; clipping to 0 (check optimizer convergence).")
    stat <- 0
  }
  tibble(
    lrt_statistic = stat,
    lrt_df        = df_diff,
    lrt_p_value   = pchisq(stat, df = df_diff, lower.tail = FALSE)
  )
}

# --- AICc-first, LRT-confirmatory pipeline for one gene's candidate set ---
select_model_pipeline <- function(gene_df,
                                   model_col  = "model",
                                   regime_col = "regime",
                                   loglik_col = "ll_total",
                                   npar_col   = "param.count",
                                   aic_col    = "AIC",
                                   null_model = "BM1",
                                   alpha      = 0.05) {

  null_row <- gene_df[gene_df[[model_col]] == null_model, ]
  if (nrow(null_row) != 1) {
    stop(sprintf("Expected exactly 1 '%s' row, found %d.", null_model, nrow(null_row)))
  }

  # Recompute the AIC-best row rather than trusting best_fit blindly
  best_row <- gene_df[which.min(gene_df[[aic_col]]), ]

  out <- tibble(
    best_regime    = best_row[[regime_col]],
    best_model     = best_row[[model_col]],
    best_AIC       = best_row[[aic_col]],
    delta_AIC_null = null_row[[aic_col]] - best_row[[aic_col]]  # >0 means best beats BM1
  )

  # Case 1: BM1 itself is AIC-best -> no LRT needed, no OU support
  if (identical(best_row[[model_col]], null_model)) {
    return(out %>% mutate(
      lrt_statistic = NA_real_,
      lrt_df        = NA_integer_,
      lrt_p_value   = NA_real_,
      decision      = "BM1 (no OU support)"
    ))
  }

  # Case 2: single confirmatory LRT of the AIC-winning model vs BM1
  lrt_res <- lrt_stat(
    ll_null = null_row[[loglik_col]], ll_alt = best_row[[loglik_col]],
    k_null  = null_row[[npar_col]],   k_alt  = best_row[[npar_col]]
  )

  out %>%
    bind_cols(lrt_res) %>%
    mutate(decision = if_else(
      lrt_p_value < alpha,
      paste0(best_model, " supported (p < ", alpha, ")"),
      "AIC favors OU but not LRT-significant vs BM1"
    ))
}

# --- apply across genes (and any extra grouping vars) + genome-wide FDR ---
run_lrt_pipeline <- function(df, group_vars = "gene_name", fdr_method = "BH", ...) {
  results <- df %>%
    group_by(across(all_of(group_vars))) %>%
    reframe(select_model_pipeline(pick(everything()), ...))

  results$lrt_p_adj <- p.adjust(results$lrt_p_value, method = fdr_method)
  results
}

# =============================================================================
# Hybrid selection: AIC/AICc for BM1 vs OU1, penalized-fit LRT for regime number
# =============================================================================
# Why not test everything against BM1 (select_model_pipeline above)?
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

  out <- tibble(
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
  p    <- pchisq(stat, df = df, lower.tail = FALSE)
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

# --- apply across genes of a SCOUT *_all_genes_full_history.csv --------------------------------
# `df` is a full_history data frame, or one or more paths to full_history CSVs (row-bound). Rows are
# reduced to the last `iter_col` per (group_vars, regime) before selection, as annotate_history does.
# fdr_method adds a genome-wide adjustment of lrt_p_corrected (NA for BM1/OU1 winners), computed over
# all rows passed in -- run once per dataset if FDR should be per dataset. With call_on_fdr = TRUE,
# multi-regime calls are re-decided on that adjusted p instead of the per-gene one.
run_hybrid_pipeline <- function(df, group_vars = c("dataset", "gene_name"), iter_col = "iter",
                                regime_col = "regime", fdr_method = "BH", call_on_fdr = FALSE, ...) {
  if (is.character(df)) df <- bind_rows(lapply(df, read.csv))
  df <- df %>% dplyr::select(-any_of(c("X", "...1")))            # write.csv row-name column
  missing <- setdiff(c(group_vars, regime_col), names(df))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
  if (iter_col %in% names(df)) {
    df <- df %>%
      group_by(across(all_of(c(group_vars, regime_col)))) %>%
      slice_max(.data[[iter_col]], n = 1, with_ties = FALSE) %>%
      ungroup()
  }

  args <- list(...)
  res  <- df %>%
    group_by(across(all_of(group_vars))) %>%
    reframe(select_model_hybrid(pick(everything()), regime_col = regime_col, ...))

  res$decided_on <- ifelse(is.na(res$lrt_p_corrected), "IC", "per-gene p")
  if (is.null(fdr_method)) return(res)
  res$lrt_p_fdr <- p.adjust(res$lrt_p_corrected, method = fdr_method)

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
