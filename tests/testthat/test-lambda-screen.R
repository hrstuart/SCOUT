## Tests for the Pagel's lambda pre-fit screen: the fast eigen engine against phytools, the
## non-ultrametric fallback, input handling that must match SCOUT/formatSCOUT, and the helpers.

suppressPackageStartupMessages(library(ape))

lambda_fixture <- function(n = 60, seed = 3, ultrametric = TRUE) {
    set.seed(seed)
    phy <- if (ultrametric) rcoal(n) else rtree(n)
    phy$tip.label <- paste0('c', seq_len(n))
    bm <- rTraitCont(phy, sigma = 1)[phy$tip.label]
    counts <- data.frame(cellBC = phy$tip.label,
                         OU2 = rep(c('a', 'b'), length.out = n),
                         bm_gene = bm,
                         bm_gene2 = 2 * bm + rnorm(n, sd = 0.05),
                         noise_gene = rnorm(n),
                         flat_gene = 5,
                         stringsAsFactors = FALSE)
    list(phy = phy, counts = counts)
}

test_that('fast engine reproduces phytools::phylosig on an ultrametric tree', {
    skip_if_not_installed('phytools')
    f <- lambda_fixture()
    a <- lambda_screen(f$counts, f$phy, species_key = 'cellBC', regimes = 'OU2',
                       normalize = FALSE, engine = 'fast', verbose = FALSE)
    b <- lambda_screen(f$counts, f$phy, species_key = 'cellBC', regimes = 'OU2',
                       normalize = FALSE, engine = 'phylosig', verbose = FALSE)
    expect_equal(attr(a, 'engine'), 'fast')
    expect_equal(attr(b, 'engine'), 'phylosig')
    for (col in c('lambda', 'lambda_logL', 'lambda_logL0', 'lambda_P')) {
        expect_equal(a[[col]], b[[col]], tolerance = 1e-6, info = col)
    }
})

test_that('BM genes score high lambda, noise scores near zero, constant genes are NA', {
    f <- lambda_fixture()
    s <- lambda_screen(f$counts, f$phy, species_key = 'cellBC', regimes = 'OU2',
                       normalize = FALSE, verbose = FALSE)
    rownames(s) <- s$gene_name
    expect_setequal(s$gene_name, c('bm_gene', 'bm_gene2', 'noise_gene', 'flat_gene'))
    expect_gt(s['bm_gene', 'lambda'], 0.8)
    expect_lt(s['bm_gene', 'lambda_P'], 1e-3)
    expect_lt(s['noise_gene', 'lambda'], 0.2)
    expect_true(is.na(s['flat_gene', 'lambda']))
    expect_true(attr(s, 'ultrametric'))
})

test_that('auto engine reports a non-ultrametric tree and uses the slow engine', {
    skip_if_not_installed('phytools')
    f <- lambda_fixture(n = 30, ultrametric = FALSE)
    expect_message(s <- lambda_screen(f$counts, f$phy, species_key = 'cellBC', regimes = 'OU2',
                                      normalize = FALSE, verbose = FALSE),
                   'not ultrametric')
    expect_equal(attr(s, 'engine'), 'phylosig')
    expect_false(attr(s, 'ultrametric'))
    expect_equal(attr(s, 'max_lambda'), 1)
    expect_error(lambda_screen(f$counts, f$phy, species_key = 'cellBC', regimes = 'OU2',
                               normalize = FALSE, engine = 'fast', verbose = FALSE), 'not ultrametric')
})

test_that('inputs are read as SCOUT reads them: file, name mangling, tip subset, log1p', {
    f <- lambda_fixture()
    cn <- f$counts
    names(cn)[names(cn) == 'bm_gene2'] <- 'HLA-DRA'
    extra <- cn[1:5, ]; extra$cellBC <- paste0('extra', 1:5)      # cells not in the tree
    cn <- rbind(cn, extra)
    cn$noise_gene <- abs(cn$noise_gene)
    cn$`HLA-DRA` <- abs(cn$`HLA-DRA`)
    tf <- tempfile(fileext = '.csv'); on.exit(unlink(tf))
    write.csv(cn, tf)                                              # leading index column
    s <- lambda_screen(tf, f$phy, species_key = 'cellBC', regimes = 'OU2',
                       genes = c('HLA-DRA', 'noise_gene'), normalize = TRUE, verbose = FALSE)
    expect_equal(s$gene_name, c('HLA.DRA', 'noise_gene'))
    expect_equal(unique(s$ntips), Ntip(f$phy))
    expect_equal(s$mean_expr[2], mean(log1p(abs(f$counts$noise_gene))), tolerance = 1e-12)
})

test_that('normalize = TRUE on negative (already log-scale) values warns', {
    f <- lambda_fixture()
    f$counts$bm_gene <- f$counts$bm_gene - 5                        # values below -1: log1p NaN
    expect_warning(s <- lambda_screen(f$counts, f$phy, species_key = 'cellBC', regimes = 'OU2',
                                      genes = 'bm_gene', verbose = FALSE), 'normalize = FALSE')
    expect_true(is.na(s$lambda))
})

test_that('missing tips are an error, not a silent subset', {
    f <- lambda_fixture()
    expect_error(lambda_screen(f$counts[-1, ], f$phy, species_key = 'cellBC', regimes = 'OU2',
                               verbose = FALSE), 'Missing counts')
})

test_that('cutoffs, select_lambda_genes and plot_lambda_screen agree', {
    f <- lambda_fixture()
    s <- lambda_screen(f$counts, f$phy, species_key = 'cellBC', regimes = 'OU2',
                       normalize = FALSE, lambda_min = 0.5, verbose = FALSE)
    expect_true('pass' %in% names(s))
    keep <- select_lambda_genes(s, lambda_min = 0.5)
    expect_setequal(keep, s$gene_name[s$pass])
    expect_true(all(c('bm_gene', 'bm_gene2') %in% keep))
    expect_false(any(c('noise_gene', 'flat_gene') %in% keep))      # NA never passes
    expect_error(select_lambda_genes(s), 'lambda_min')
    pdf(NULL); on.exit(dev.off())
    p <- plot_lambda_screen(s, lambda_min = 0.5)
    expect_equal(p$n_kept, length(keep))
    expect_equal(p$n_na, 1)
})

test_that('formatSCOUT accepts a single screened gene', {
    f <- lambda_fixture()
    out <- formatSCOUT(f$phy, f$counts, outpath = tempdir(), species_key = 'cellBC',
                       quant_traits = 'bm_gene', regimes = 'BM1', anc_infer = 'skip')
    expect_equal(out$gene_cols, 'bm_gene')
})
