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
