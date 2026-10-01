## Tests for the hybrid IC + LRT model selection (select_model_hybrid / run_hybrid_pipeline).

## One gene's candidate set. k: BM1 2, OU1 3, OU2 4, OU3 5 (thetas add one each).
gene_rows <- function(gene, ll, ntips = 200, dataset = 'd1') {
    data.frame(dataset = dataset, gene_name = gene,
               regime = c('BM1', 'OU1', 'OU2', 'OU3'),
               ll_total = ll, param.count = c(2, 3, 4, 5), ntips = ntips, iter = 10,
               stringsAsFactors = FALSE)
}

test_that('BM1 and OU1 winners are decided by IC alone', {
    bm <- select_model_hybrid(gene_rows('g', c(-100, -100.2, -100.1, -100)))
    expect_equal(bm$call, 'BM1')
    expect_equal(bm$lrt_status, 'not tested')
    ou <- select_model_hybrid(gene_rows('g', c(-110, -100, -99.9, -99.8)))
    expect_equal(ou$call, 'OU1')
    expect_equal(ou$lrt_status, 'not eligible')
    expect_true(is.na(ou$lrt_p_value))
})

test_that('a multi-regime winner is kept only if its Bonferroni LRT vs OU1 is significant', {
    strong <- select_model_hybrid(gene_rows('g', c(-120, -110, -95, -95)))
    expect_equal(strong$ic_winner, 'OU2')
    expect_equal(strong$call, 'OU2')
    expect_equal(strong$lrt_df, 1L)
    expect_equal(strong$lrt_statistic, 30)
    expect_equal(strong$lrt_p_corrected, min(1, 2 * pchisq(30, 1, lower.tail = FALSE)))
    expect_equal(strong$lrt_status, 'tested: significant')

    # LR = 2 * 2.2 = 4.4: p ~ 0.036 uncorrected (wins by AIC), 0.072 after x2 -> not significant
    weak <- gene_rows('g', c(-120, -110, -107.8, -110))
    w <- select_model_hybrid(weak, ic = 'AIC')
    expect_equal(w$ic_winner, 'OU2')
    expect_lt(w$lrt_p_value, 0.05)
    expect_gt(w$lrt_p_corrected, 0.05)
    expect_equal(w$call, 'OU1')                                     # fallback to best of BM1/OU1
    expect_equal(select_model_hybrid(weak, ic = 'AIC', on_fail = 'ambiguous')$call, 'ambiguous')
    expect_equal(select_model_hybrid(weak, ic = 'AIC', correction = 'none')$call, 'OU2')
})

test_that('AICc and AIC are recomputed from ll, k and ntips', {
    r <- gene_rows('g', c(-100, -99, -98, -97), ntips = 12)
    a <- select_model_hybrid(r, ic = 'AICc')
    crit <- -2 * r$ll_total + 2 * r$param.count * 12 / (12 - r$param.count - 1)
    expect_equal(a$ic_winner, r$regime[which.min(crit)])
})

test_that('bad candidate sets are errors', {
    r <- gene_rows('g', c(-100, -99, -98, -97))
    expect_error(select_model_hybrid(rbind(r, r)), 'More than one row')
    expect_error(select_model_hybrid(r[r$regime != 'OU1', ]), "No 'OU1' row")
})

test_that('run_hybrid_pipeline keeps the last iteration and adds FDR', {
    early <- gene_rows('g1', c(-200, -200, -200, -200)); early$iter <- 1
    df <- rbind(early, gene_rows('g1', c(-120, -110, -95, -95)),
                gene_rows('g2', c(-100, -100.2, -100.1, -100)),
                gene_rows('g3', c(-120, -110, -107.8, -110)))
    df$X <- seq_len(nrow(df))                                       # write.csv row-name column
    res <- run_hybrid_pipeline(df, ic = 'AIC')
    expect_equal(nrow(res), 3)
    expect_equal(res$call, c('OU2', 'BM1', 'OU1'))
    expect_equal(res$decided_on, c('per-gene p', 'IC', 'per-gene p'))
    expect_equal(res$lrt_p_fdr, p.adjust(res$lrt_p_corrected, 'BH'))

    fdr <- run_hybrid_pipeline(df, ic = 'AIC', call_on_fdr = TRUE)
    expect_equal(fdr$decided_on[!is.na(fdr$lrt_p_fdr)], rep('FDR p (BH)', 2))
    expect_false('lrt_p_fdr' %in% names(run_hybrid_pipeline(df, ic = 'AIC', fdr_method = NULL)))

    tf <- tempfile(fileext = '.csv'); on.exit(unlink(tf))
    write.csv(df[, names(df) != 'X'], tf)
    expect_equal(run_hybrid_pipeline(tf, ic = 'AIC')$call, res$call)
})

## ---- scout_model_selection ------------------------------------------------------------

## full-history rows with fit diagnostics; converge defaults to 'converged'
fit_rows <- function(gene, ll, ntips = 200, converge = 'converged', iter = 8) {
    r <- gene_rows(gene, ll, ntips)
    r$model <- c('BM1', 'OU1', 'OUM', 'OUM')
    r$converge <- converge; r$iter <- iter
    r$alpha <- c(1e-10, 0.5, 0.6, 0.7); r$sigma <- 1:4; r$tau <- 0.2
    r
}

test_that('AIC and AICc select differently on a small tree; weights and deltas are consistent', {
    r <- fit_rows('g', c(-100, -98, -96.5, -96), ntips = 12)
    a  <- scout_model_selection(r, ic = 'AIC',  hybrid = FALSE)
    ac <- scout_model_selection(r, ic = 'AICc', hybrid = FALSE)
    expect_equal(a$selected_regime, 'OU2')                          # AIC 204 202 201 202
    expect_equal(ac$selected_regime, 'OU1')                         # AICc 205.3 205 206.7 212
    expect_equal(a$AIC, 201)
    expect_equal(a$delta_AIC_next, 1)
    expect_equal(a$next_best_regime, 'OU1')                         # first of the tie at 202
    w <- exp(-0.5 * (c(204, 202, 201, 202) - 201)); w <- w / sum(w)
    expect_equal(a$AIC_weight, w[3])
    expect_true(all(c('AICc', 'AICc_weight', 'delta_AICc_next') %in% names(ac)))
    expect_false(any(grepl('^lrt_', names(a))))
    expect_equal(a$alpha, 0.6); expect_equal(a$sigma, 3)
    expect_true(a$final_set)
})

test_that('hybrid selection matches select_model_hybrid and can sit below the IC winner', {
    r <- fit_rows('g', c(-120, -110, -107.8, -110))
    s <- scout_model_selection(r, ic = 'AIC', hybrid = TRUE)
    h <- select_model_hybrid(r, ic = 'AIC')
    expect_equal(s$selected_regime, h$call)                         # OU1 (fallback)
    expect_equal(s$ic_winner, 'OU2')
    expect_equal(s$lrt_p_corrected, h$lrt_p_corrected)
    expect_equal(s$delta_AIC_next, 223.6 - 226)                     # negative: winner was OU2
    amb <- scout_model_selection(r, ic = 'AIC', hybrid_on_fail = 'ambiguous')
    expect_equal(amb$selected_regime, 'ambiguous')
    expect_true(amb$flag_ambiguous)
    expect_false(amb$final_set)
    expect_equal(amb$alpha, 0.6)                                    # parameters of the IC winner
    expect_error(scout_model_selection(r[r$regime != 'OU1', ]), 'needs BM1 and OU1')
})

test_that('BM1 alpha is reported as NA', {
    s <- scout_model_selection(fit_rows('g', c(-100, -100.2, -100.1, -100)), hybrid = FALSE)
    expect_equal(s$selected_model, 'BM1')
    expect_true(is.na(s$alpha))
})

test_that('convergence flags, flag_scope, check_convergence and exclude', {
    strong <- c(-120, -110, -95, -95)                               # OU2 selected
    capped <- fit_rows('capped', strong, converge = NA, iter = 100)
    osc    <- fit_rows('osc', strong, converge = 'early_stop_oscillating', iter = 40)
    failed <- fit_rows('failed', strong, converge = NA, iter = 3)
    side   <- fit_rows('side', strong); side$converge[1] <- NA; side$iter[1] <- 100  # BM1 only
    s <- scout_model_selection(rbind(capped, osc, failed, side), ic = 'AIC', hybrid = FALSE)
    rownames(s) <- s$gene_name
    expect_equal(s['capped', c('flag_not_converged', 'flag_max_iter', 'flag_oscillating')],
                 data.frame(flag_not_converged = TRUE, flag_max_iter = TRUE, flag_oscillating = FALSE,
                            row.names = 'capped'))
    expect_true(s['osc', 'flag_oscillating'] && s['osc', 'flag_not_converged'])
    expect_false(s['osc', 'flag_max_iter'])
    expect_true(s['failed', 'flag_not_converged']); expect_false(s['failed', 'flag_max_iter'])
    expect_false(s['side', 'flag_not_converged'])                   # selected OU2 converged
    expect_equal(s['side', 'n_candidates_not_converged'], 1)
    expect_equal(s$final_set, c(FALSE, FALSE, FALSE, TRUE))

    all_scope <- scout_model_selection(side, ic = 'AIC', hybrid = FALSE, flag_scope = 'all')
    expect_true(all_scope$flag_max_iter); expect_false(all_scope$final_set)
    off <- scout_model_selection(capped, ic = 'AIC', hybrid = FALSE, check_convergence = FALSE)
    expect_false(off$flag_not_converged); expect_true(off$final_set)
    keep_osc <- scout_model_selection(osc, ic = 'AIC', hybrid = FALSE, exclude = 'max_iter')
    expect_true(keep_osc$flag_oscillating); expect_true(keep_osc$final_set)
    expect_error(scout_model_selection(osc, exclude = 'nope'), 'Unknown flag')
})

test_that('a missing or non-finite candidate fit is flagged', {
    a <- fit_rows('a', c(-120, -110, -95, -95))
    b <- fit_rows('b', c(-120, -110, -95, NA))
    s <- scout_model_selection(rbind(a, b[b$regime != 'OU2', ], fit_rows('c', c(-120, -110, -95, -95))[-3, ]),
                               ic = 'AIC', hybrid = FALSE)
    expect_equal(s$flag_missing_fit, c(FALSE, TRUE, TRUE))
    expect_equal(s$final_set, c(TRUE, FALSE, FALSE))
})

test_that('theta columns come from params, as a list or as SCOUT parameter files', {
    r <- rbind(fit_rows('g1', c(-120, -110, -95, -95)), fit_rows('g2', c(-100, -100.2, -100.1, -100)))
    pars <- list(BM1 = data.frame(gene_name = c('g1', 'g2'), theta_1 = c(1, 2)),
                 OU1 = data.frame(gene_name = c('g1', 'g2'), theta_1 = c(3, 4)),
                 OU2 = data.frame(gene_name = c('g1', 'g2'), theta_a = c(5, 6), theta_b = c(7, 8)),
                 OU3 = data.frame(gene_name = c('g1', 'g2'), theta_a = 9, theta_b = 9, theta_c = 9))
    s <- scout_model_selection(r, params = pars, ic = 'AIC', hybrid = FALSE)
    expect_equal(s$selected_regime, c('OU2', 'BM1'))
    expect_equal(s$theta_a, c(5, NA)); expect_equal(s$theta_b, c(7, NA))
    expect_equal(s$theta_1, c(NA, 2)); expect_true(all(is.na(s$theta_c)))

    d <- tempfile(); dir.create(d); on.exit(unlink(d, recursive = TRUE))
    files <- vapply(names(pars), function(rg) {
        f <- file.path(d, sprintf('run_all_genes_%s_parameters.csv', rg))
        write.table(pars[[rg]], f); f }, '')
    s2 <- scout_model_selection(r, params = unname(files), ic = 'AIC', hybrid = FALSE)
    expect_equal(s2[, c('theta_1', 'theta_a', 'theta_b')], s[, c('theta_1', 'theta_a', 'theta_b')])
    r$dataset[r$gene_name == 'g2'] <- 'd2'
    expect_error(scout_model_selection(r, params = pars, hybrid = FALSE), 'single dataset')
})
