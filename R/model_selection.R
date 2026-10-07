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
#
# Speed: the per-gene decision (hybrid_decide) works on plain vectors and returns a plain list.
# Building a one-row tibble per gene cost ~4 ms and dominated whole-genome runs; the bulk paths
# (run_hybrid_pipeline, scout_model_selection) assemble each output column once instead.

#' Resolve select_model_hybrid() options (match.arg and defaults) once per call.
#' @noRd
hybrid_args <- function(ic          = c("AICc", "AIC"),
                        on_fail     = c("fallback", "ambiguous"),
                        alpha       = 0.05,
                        correction  = c("bonferroni", "none"),
                        bm_regime   = "BM1",
                        ou_regime   = "OU1",
                        loglik_col  = "ll_total",
                        npar_col    = "param.count",
                        ntips_col   = "ntips") {
  list(ic = match.arg(ic), on_fail = match.arg(on_fail), alpha = alpha,
       correction = match.arg(correction), bm_regime = bm_regime, ou_regime = ou_regime,
       loglik_col = loglik_col, npar_col = npar_col, ntips_col = ntips_col)
}

#' The hybrid decision for one gene's candidate set, as a plain list (one element per output
#' column of select_model_hybrid()). `reg`, `ll`, `k` are per-regime vectors, `n` the tip count.
#' @noRd
hybrid_decide <- function(reg, ll, k, n, a) {
  if (anyDuplicated(reg)) {
    stop("More than one row per regime (e.g. several EM iterations); keep the last iter per regime first.")
  }
  for (r in c(a$bm_regime, a$ou_regime)) {
    if (!r %in% reg) stop(sprintf("No '%s' row for this gene.", r))
  }
  bm <- a$bm_regime; ou <- a$ou_regime; ic <- a$ic

  crit <- if (ic == "AIC") -2 * ll + 2 * k else -2 * ll + 2 * k * n / (n - k - 1)
  names(crit) <- names(ll) <- names(k) <- reg

  n_multi <- sum(!reg %in% c(bm, ou))
  winner  <- reg[which.min(crit)]
  simple  <- if (crit[bm] <= crit[ou]) bm else ou

  out <- list(
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
  if (winner == bm) {
    out$lrt_status <- "not tested"
    out$flag       <- sprintf("%s best by %s; simplest candidate, no null to test against", bm, ic)
    return(out)
  }
  if (winner == ou) {
    out$lrt_status <- "not eligible"
    out$flag       <- sprintf("%s best by %s; LRT vs %s not valid (not nested under SCOUT's alpha penalty)",
                              ou, ic, bm)
    return(out)
  }

  raw  <- 2 * (ll[[winner]] - ll[[ou]])
  stat <- max(raw, 0)
  df   <- as.integer(k[[winner]] - k[[ou]])
  if (df <= 0) stop(sprintf("'%s' does not have more parameters than '%s'.", winner, ou))
  p    <- stats::pchisq(stat, df = df, lower.tail = FALSE)
  pc   <- if (a$correction == "bonferroni") min(1, p * n_multi) else p

  out$lrt_alt           <- winner
  out$lrt_null          <- ou
  out$lrt_statistic     <- stat
  out$lrt_statistic_raw <- raw
  out$lrt_df            <- df
  out$lrt_p_value       <- p
  out$lrt_p_corrected   <- pc
  out$lrt_status        <- if (pc < a$alpha) "tested: significant" else "tested: not significant"
  flags <- character(0)
  if (raw < 0) flags <- c(flags, "negative LR clipped to 0")
  if (pc >= a$alpha) {
    out$call <- if (a$on_fail == "fallback") simple else "ambiguous"
    flags    <- c(flags, sprintf("%s not significant vs %s", winner, ou))
  }
  out$flag <- paste(flags, collapse = "; ")
  out
}

#' Keep the last `iter_col` row per `cols` group, as
#' group_by(cols) |> slice_max(iter, n = 1, with_ties = FALSE) |> ungroup() does -- same rows, same
#' order (groups sorted, strings in C-locale order; ties keep the earlier row) -- in one radix sort
#' instead of a per-group evaluation, which took ~30 s on 150k groups.
#' @noRd
last_iter_rows <- function(df, cols, iter_col) {
  o <- do.call(order, c(unname(as.list(df[cols])), list(df[[iter_col]]),
                        list(decreasing = c(rep(FALSE, length(cols)), TRUE), method = "radix")))
  keep <- o[!duplicated(df[o, cols, drop = FALSE])]
  out <- df[keep, , drop = FALSE]
  if (!inherits(out, "tbl_df")) rownames(out) <- NULL
  out
}

#' Column template for hybrid_decide() results (names and types), for vapply.
#' @noRd
HYBRID_TEMPLATE <- list(ic = "", ic_winner = "", best_simple = "", delta_ic_simple = 0,
                        n_multi = 0L, lrt_alt = "", lrt_null = "", lrt_statistic = 0,
                        lrt_statistic_raw = 0, lrt_df = 0L, lrt_p_value = 0,
                        lrt_p_corrected = 0, call = "", lrt_status = "", flag = "")

#' Turn a list of per-row named lists into a list of columns, one vapply per column.
#' @noRd
rows_to_columns <- function(rows, template) {
  lapply(stats::setNames(names(template), names(template)), function(cn)
    vapply(rows, function(r) r[[cn]], template[[cn]], USE.NAMES = FALSE))
}

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
  a <- hybrid_args(ic, on_fail, alpha, correction, bm_regime, ou_regime,
                   loglik_col, npar_col, ntips_col)
  dplyr::as_tibble(hybrid_decide(as.character(gene_df[[regime_col]]), gene_df[[loglik_col]],
                                 gene_df[[npar_col]], gene_df[[ntips_col]][1], a))
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
  if (iter_col %in% names(df)) df <- last_iter_rows(df, c(group_vars, regime_col), iter_col)

  a  <- do.call(hybrid_args, list(...))
  gd <- dplyr::group_by(df, dplyr::across(dplyr::all_of(group_vars)))
  idx <- dplyr::group_rows(gd)                                   # dplyr's group order
  reg <- as.character(df[[regime_col]]); ll <- df[[a$loglik_col]]
  k   <- df[[a$npar_col]];               n  <- df[[a$ntips_col]]
  rows <- lapply(idx, function(i) hybrid_decide(reg[i], ll[i], k[i], n[i][1], a))
  res  <- dplyr::bind_cols(dplyr::group_keys(gd),
                           dplyr::as_tibble(rows_to_columns(rows, HYBRID_TEMPLATE)))

  res$decided_on <- ifelse(is.na(res$lrt_p_corrected), "IC", "per-gene p")
  if (is.null(fdr_method)) return(res)
  res$lrt_p_fdr <- stats::p.adjust(res$lrt_p_corrected, method = fdr_method)

  if (call_on_fdr) {
    # flag keeps the per-gene verdict; decided_on marks calls re-made on the FDR-adjusted p
    tested <- !is.na(res$lrt_p_fdr)
    res$call[tested] <- ifelse(res$lrt_p_fdr[tested] < a$alpha, res$ic_winner[tested],
                               if (a$on_fail == "fallback") res$best_simple[tested] else "ambiguous")
    res$lrt_status[tested] <- ifelse(res$lrt_p_fdr[tested] < a$alpha, "tested: significant", "tested: not significant")
    res$decided_on[tested] <- paste0("FDR p (", fdr_method, ")")
  }
  res
}

# =============================================================================
# Per-gene model selection table: IC or hybrid call, weights, parameters, QC flags
# =============================================================================

#' Read SCOUT parameter tables into a named list of per-regime data.frames.
#' Accepts SCOUT()$SCOUT_params (a named list) or paths to *_all_genes_<regime>_parameters.csv.
#' @noRd
read_scout_params <- function(params) {
  if (is.null(params)) return(NULL)
  if (is.data.frame(params)) stop("`params` must be a named list of per-regime data.frames or file paths.")
  if (is.character(params)) {
    regs <- names(params)
    if (is.null(regs)) regs <- sub("^.*_all_genes_(.+)_parameters\\.csv$", "\\1", basename(params))
    params <- stats::setNames(lapply(params, function(f) utils::read.table(f, header = TRUE)), regs)
  }
  if (is.null(names(params)) || any(names(params) == "")) stop("`params` must be named by regime.")
  params
}

#' Per-gene model selection with QC flags
#'
#' Builds one row per gene from SCOUT's full fit history: the selected model (by AIC/AICc alone, or
#' by the hybrid IC + LRT procedure of \code{select_model_hybrid()}), its information-criterion
#' weight and margin over the next-best candidate, its fitted parameters, and one logical column
#' per quality flag. \code{final_set} is TRUE for genes that raise none of the flags named in
#' \code{exclude}.
#'
#' Convergence comes from the EM \code{converge} column: \code{'converged'},
#' \code{'early_stop_oscillating'}, or NA (stopped without converging -- the iteration cap, or a
#' failed iteration). \code{runSCOUT()} runs at most 100 EM iterations, which is not recorded in the
#' output, so a fit is taken to have hit the cap when \code{iter >= max_iter} without converging.
#' Single-pass fits (\code{method = 'MTF'} or \code{'SM'}) never converge in this sense; use
#' \code{check_convergence = FALSE} for them.
#'
#' @param history SCOUT full history: a data.frame (e.g. \code{*_all_genes_full_history.csv}, or
#'   the \code{annotated} table inside \code{SCOUT()}), or one or more paths to such CSVs.
#' @param params optional fitted parameters, for the theta columns: \code{SCOUT()$SCOUT_params}
#'   or paths to \code{*_all_genes_<regime>_parameters.csv} (named by regime, or named as SCOUT
#'   writes them). Needs a single dataset in \code{history}.
#' @param ic \code{'AICc'} or \code{'AIC'}, recomputed from ll_total, param.count and ntips.
#' @param hybrid if TRUE, select with \code{select_model_hybrid()} (needs BM1 and OU1 candidates);
#'   if FALSE, the selected model is simply the \code{ic}-best.
#' @param hybrid_alpha,hybrid_on_fail,hybrid_correction passed to \code{select_model_hybrid()}
#'   as \code{alpha}, \code{on_fail}, \code{correction}.
#' @param fdr_method \code{p.adjust} method for \code{lrt_p_fdr} across genes (hybrid only).
#' @param max_iter EM iteration cap used in the fit (100 in \code{runSCOUT()}).
#' @param check_convergence if FALSE, the convergence flags are all FALSE (single-pass fits).
#' @param flag_scope \code{'selected'}: convergence flags describe the selected model only;
#'   \code{'all'}: they are raised if ANY candidate model of the gene failed, since a
#'   non-converged competitor's likelihood also distorts the comparison.
#' @param exclude flags that remove a gene from \code{final_set}: any of \code{'not_converged'},
#'   \code{'max_iter'}, \code{'oscillating'}, \code{'ambiguous'}, \code{'missing_fit'}.
#' @param group_vars,regime_col,iter_col columns identifying a gene, its candidate regime, and the
#'   EM iteration (only the last iteration per regime is used).
#' @return A data.frame, one row per gene: selected_regime, selected_model, next_best_regime;
#'   \code{<ic>}, \code{<ic>_weight}, \code{delta_<ic>_next} (next-best minus selected; negative
#'   when the hybrid kept a simpler model than the IC winner); ll_total, param.count, ntips,
#'   alpha (NA for BM1), sigma, tau, theta_* (with \code{params}); iter, converge,
#'   n_candidates_not_converged; hybrid columns (ic_winner, lrt_*); the flag_* columns and
#'   final_set.
#' @examples
#' \dontrun{
#' sel <- scout_model_selection('out/run_all_genes_full_history.csv',
#'          params = Sys.glob('out/run_all_genes_*_parameters.csv'), ic = 'AICc', hybrid = TRUE)
#' table(sel$selected_regime[sel$final_set])
#' }
#' @export
scout_model_selection <- function(history,
                                  params            = NULL,
                                  ic                = c("AICc", "AIC"),
                                  hybrid            = TRUE,
                                  hybrid_alpha      = 0.05,
                                  hybrid_on_fail    = c("fallback", "ambiguous"),
                                  hybrid_correction = c("bonferroni", "none"),
                                  fdr_method        = "BH",
                                  max_iter          = 100,
                                  check_convergence = TRUE,
                                  flag_scope        = c("selected", "all"),
                                  exclude           = c("not_converged", "max_iter", "oscillating",
                                                        "ambiguous", "missing_fit"),
                                  group_vars        = c("dataset", "gene_name"),
                                  regime_col        = "regime",
                                  iter_col          = "iter") {
  ic                <- match.arg(ic)
  hybrid_on_fail    <- match.arg(hybrid_on_fail)
  hybrid_correction <- match.arg(hybrid_correction)
  flag_scope        <- match.arg(flag_scope)
  all_flags <- c("not_converged", "max_iter", "oscillating", "ambiguous", "missing_fit")
  bad <- setdiff(exclude, all_flags)
  if (length(bad)) stop("Unknown flag(s) in `exclude`: ", paste(bad, collapse = ", "))

  df <- if (is.character(history)) dplyr::bind_rows(lapply(history, utils::read.csv)) else history
  df <- as.data.frame(dplyr::ungroup(df))
  df <- df[, setdiff(names(df), c("X", "...1")), drop = FALSE]       # write.csv row-name column
  group_vars <- intersect(group_vars, names(df))                       # 'dataset' is optional
  need <- c("gene_name", regime_col, "ll_total", "param.count", "ntips")
  miss <- setdiff(need, names(df))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
  if (!"converge" %in% names(df)) df$converge <- NA_character_
  if (!"model" %in% names(df)) df$model <- df[[regime_col]]
  if (iter_col %in% names(df)) {
    df <- last_iter_rows(df, c(group_vars, regime_col), iter_col)
  } else {
    df[[iter_col]] <- NA_integer_
  }

  candidates <- sort(unique(as.character(df[[regime_col]])))
  if (hybrid && !all(c("BM1", "OU1") %in% candidates)) {
    stop("hybrid = TRUE needs BM1 and OU1 among the fitted regimes; use hybrid = FALSE.")
  }
  a <- hybrid_args(ic = ic, on_fail = hybrid_on_fail, alpha = hybrid_alpha,
                   correction = hybrid_correction)

  par_list <- read_scout_params(params)
  if (!is.null(par_list) && "dataset" %in% names(df) && length(unique(df$dataset)) > 1) {
    stop("`params` can only be matched to a history with a single dataset.")
  }
  theta_cols <- unique(unlist(lapply(par_list, function(p) grep("^theta_", names(p), value = TRUE))))
  # theta lookup: regime -> theta column -> value by gene (first row per gene, as before)
  theta_lut <- lapply(par_list, function(p) {
    first <- !duplicated(as.character(p$gene_name))
    lapply(stats::setNames(intersect(theta_cols, names(p)), intersect(theta_cols, names(p))),
           function(tc) stats::setNames(p[[tc]][first], as.character(p$gene_name)[first]))
  })

  # Column vectors, indexed per gene below -- no per-gene data.frame subsetting.
  reg_all   <- as.character(df[[regime_col]])
  ll_all    <- df$ll_total
  k_all     <- df$param.count
  n_all     <- df$ntips
  conv_all  <- as.character(df$converge)
  it_all    <- df[[iter_col]]
  model_all <- as.character(df$model)
  gene_all  <- as.character(df$gene_name)
  alpha_all <- if ("alpha" %in% names(df)) df$alpha else NULL
  sigma_all <- if ("sigma" %in% names(df)) df$sigma else NULL
  tau_all   <- if ("tau" %in% names(df)) df$tau else NULL

  not_conv <- function(cv) is.na(cv) | cv != "converged"
  hit_cap  <- function(cv, it) is.na(cv) & !is.na(it) & it >= max_iter
  oscill   <- function(cv) !is.na(cv) & cv == "early_stop_oscillating"

  ic_col <- ic; w_col <- paste0(ic, "_weight"); d_col <- paste0("delta_", ic, "_next")
  hyb_cols <- c("ic_winner", "lrt_status", "lrt_alt", "lrt_statistic", "lrt_df",
                "lrt_p_value", "lrt_p_corrected")
  # Output columns after the group variables, in order, with their NA of the right type.
  template <- c(
    list(selected_regime = NA_character_, selected_model = NA_character_,
         next_best_regime = NA_character_),
    stats::setNames(list(NA_real_, NA_real_, NA_real_), c(ic_col, w_col, d_col)),
    list(ll_total = ll_all[NA_integer_][1], param.count = k_all[NA_integer_][1],
         ntips = n_all[NA_integer_][1], alpha = NA_real_, sigma = NA_real_, tau = NA_real_),
    stats::setNames(rep(list(NA_real_), length(theta_cols)), theta_cols),
    list(iter = it_all[NA_integer_][1], converge = NA_character_,
         n_candidates_not_converged = NA_integer_),
    if (hybrid) list(ic_winner = NA_character_, lrt_status = NA_character_,
                     lrt_alt = NA_character_, lrt_statistic = NA_real_, lrt_df = NA_integer_,
                     lrt_p_value = NA_real_, lrt_p_corrected = NA_real_),
    stats::setNames(rep(list(FALSE), length(all_flags)), paste0("flag_", all_flags)))

  one_gene <- function(ix) {
    row <- template
    reg <- reg_all[ix]; ll <- ll_all[ix]; k <- k_all[ix]; n <- n_all[ix][1]
    crit <- if (ic == "AIC") -2 * ll + 2 * k else -2 * ll + 2 * k * n / (n - k - 1)
    names(crit) <- reg
    ok <- is.finite(crit)
    row$flag_missing_fit <- !all(candidates %in% reg[ok])
    if (!any(ok)) return(row)                                # no usable fit: all NA, flagged
    w <- rep(NA_real_, length(crit))
    d <- crit[ok] - min(crit[ok])
    w[ok] <- exp(-0.5 * d) / sum(exp(-0.5 * d))
    ic_winner <- reg[ok][which.min(crit[ok])]

    hyb <- NULL
    selected <- ic_winner
    if (hybrid) {
      hyb <- tryCatch(hybrid_decide(reg[ok], ll[ok], k[ok], n_all[ix[ok]][1], a),
                      error = function(e) NULL)
      if (!is.null(hyb)) selected <- hyb$call
    }
    ambiguous <- identical(selected, "ambiguous")
    shown <- if (ambiguous) ic_winner else selected          # whose parameters are reported
    i <- which(reg == shown)
    others <- ok & reg != shown
    nb <- if (any(others)) reg[others][which.min(crit[others])] else NA_character_
    j <- ix[i]                                               # row of the reported fit
    mdl <- model_all[j]

    row$selected_regime  <- selected
    row$selected_model   <- if (ambiguous) NA_character_ else mdl
    row$next_best_regime <- nb
    row[[ic_col]] <- unname(crit[i])
    row[[w_col]]  <- w[i]
    row[[d_col]]  <- if (is.na(nb)) NA_real_ else unname(crit[nb] - crit[i])
    row$ll_total    <- ll[i]
    row$param.count <- k[i]
    row$ntips       <- n
    row$alpha <- if (mdl == "BM1" || is.null(alpha_all)) NA_real_ else alpha_all[j]
    row$sigma <- if (is.null(sigma_all)) NA_real_ else sigma_all[j]
    row$tau   <- if (is.null(tau_all)) NA_real_ else tau_all[j]
    if (length(theta_cols)) {
      lut <- theta_lut[[shown]]
      for (tc in names(lut)) {
        v <- lut[[tc]][gene_all[ix[1]]]
        if (!is.na(names(v))) row[[tc]] <- unname(v)
      }
    }
    row$iter     <- it_all[j]
    row$converge <- conv_all[j]
    row$n_candidates_not_converged <- sum(not_conv(conv_all[ix]))

    if (!is.null(hyb)) {
      for (cc in hyb_cols) row[[cc]] <- hyb[[cc]]
    } else if (hybrid) {
      row$ic_winner  <- ic_winner
      row$lrt_status <- "hybrid failed"
    }

    sc <- if (flag_scope == "all") ix else j
    cv <- conv_all[sc]; it <- it_all[sc]
    row$flag_not_converged <- check_convergence && any(not_conv(cv))
    row$flag_max_iter      <- check_convergence && any(hit_cap(cv, it))
    row$flag_oscillating   <- check_convergence && any(oscill(cv))
    row$flag_ambiguous     <- ambiguous
    row
  }

  # Same grouping and group order as before: one split of row indices, not of the data.frame.
  keys <- interaction(df[group_vars], drop = TRUE, lex.order = TRUE)
  idx  <- split(seq_len(nrow(df)), keys)
  rows <- lapply(idx, one_gene)
  first <- vapply(idx, `[`, 1L, 1L, USE.NAMES = FALSE)
  out <- df[first, group_vars, drop = FALSE]
  cols <- lapply(stats::setNames(names(template), names(template)), function(cn)
    unlist(lapply(rows, `[[`, cn), use.names = FALSE))
  for (cn in names(cols)) out[[cn]] <- cols[[cn]]

  if (hybrid && !is.null(fdr_method)) {
    out$lrt_p_fdr <- stats::p.adjust(out$lrt_p_corrected, method = fdr_method)
  }
  flag_mat <- as.matrix(out[, paste0("flag_", exclude), drop = FALSE])
  out$final_set <- if (length(exclude)) rowSums(flag_mat) == 0 else TRUE
  out$ic <- ic
  rownames(out) <- NULL
  as.data.frame(out)
}
