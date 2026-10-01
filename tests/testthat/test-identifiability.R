## Tests for the model-class identifiability sweep. These cover the tree-variant construction and
## the analytic reference; the empirical sweep itself is too slow for a test suite.

suppressPackageStartupMessages({
    library(ape)
    library(castor)
})

id_fixture <- function(n = 40, seed = 11) {
    set.seed(seed)
    P <- matrix(c(0.90, 0.08, 0.02,
                  0.10, 0.80, 0.10,
                  0.02, 0.08, 0.90), 3, 3, byrow = TRUE)
    phy <- generate_tree_hbd_reverse(n, rho = 1, lambda = 1, mu = 0)$trees[[1]]
    phy$tip.label <- paste0('t', seq_len(n))
    st <- simulate_mk_model(phy, Q = expm::logm(P), include_nodes = FALSE)
    phy$states <- paste0('state', st$tip_states)
    SCOUT:::prep_dropout_tree(phy, verbose = FALSE)
}

phy <- id_fixture()

## ---- tree variants --------------------------------------------------------------------

test_that('tree_variant rescales to the requested height and preserves the topology', {
    for (br in c('keep', 'unit')) {
        for (dp in c('keep', 'equalise', 'vary')) {
            v <- tree_variant(phy, br, dp, height = 7, randseed = 1)
            expect_equal(max(node.depth.edgelength(v)[seq_len(Ntip(v))]), 7, tolerance = 1e-8,
                         info = paste(br, dp))
            expect_identical(v$edge, phy$edge)          # topology untouched
            expect_identical(v$tip.label, phy$tip.label)
            expect_true(all(v$edge.length > 0))
        }
    }
})

test_that('the factorial arms have the intended branch-length structure', {
    v <- tree_variant_set(phy, height = 7, randseed = 1)
    expect_named(v, c('yule_ultra', 'unit_ultra', 'unit_vary', 'yule_vary'))

    internal_cv <- function(p) {
        e <- p$edge.length[p$edge[, 2] > Ntip(p)]
        stats::sd(e) / mean(e)
    }
    # unit arms have uniform internal branches; yule arms do not
    expect_equal(internal_cv(v$unit_ultra), 0, tolerance = 1e-10)
    expect_equal(internal_cv(v$unit_vary), 0, tolerance = 1e-10)
    expect_gt(internal_cv(v$yule_ultra), 0.3)
    expect_gt(internal_cv(v$yule_vary), 0.3)

    # ultra arms are ultrametric; vary arms are not
    expect_true(is.ultrametric(v$yule_ultra))
    expect_true(is.ultrametric(v$unit_ultra))
    expect_false(is.ultrametric(v$unit_vary))
    expect_false(is.ultrametric(v$yule_vary))

    # every arm shares the imposed height
    for (nm in names(v)) {
        expect_equal(max(node.depth.edgelength(v[[nm]])[seq_len(Ntip(phy))]), 7, tolerance = 1e-8)
    }
})

test_that('a tree with no branch lengths is treated as unit-length, matching formatSCOUT', {
    bare <- phy
    bare$edge.length <- NULL
    v <- tree_variant(bare, 'keep', 'keep', height = NULL)
    expect_true(all(v$edge.length == 1))
})

## ---- analytic identifiability ---------------------------------------------------------

test_that('BM and OU become indistinguishable in shape as alpha*H approaches zero', {
    v <- tree_variant_set(phy, height = 7, randseed = 1)
    r <- identifiability_reference(v$yule_ultra, alphaH = c(0.1, 0.5, 2, 10, 40))

    # correlation structures coincide in the small-alpha limit ...
    expect_gt(r$cor_shape[r$alphaH == 0.1], 0.99)
    # ... and separate as alpha*H grows
    expect_lt(r$cor_shape[r$alphaH == 40], 0.4)
    expect_true(all(diff(r$cor_shape) < 0))          # monotone decreasing
    expect_true(all(is.finite(r$resid_sd)))
    expect_true(all(r$resid_sd >= 0))
})

test_that('the tip-variance discriminator exists only when tip depths vary', {
    v <- tree_variant_set(phy, height = 7, randseed = 1)
    ah <- c(1, 5, 20)

    for (nm in c('yule_ultra', 'unit_ultra')) {
        r <- identifiability_reference(v[[nm]], alphaH = ah)
        # BM tip variance is sigma^2 * depth, so on an ultrametric tree it is constant
        expect_equal(unique(r$cv_var_bm), 0, tolerance = 1e-8, info = nm)
    }
    for (nm in c('unit_vary', 'yule_vary')) {
        r <- identifiability_reference(v[[nm]], alphaH = ah)
        expect_gt(unique(r$cv_var_bm), 0.02, label = nm)
    }
    # the OU stationary variance never depends on depth
    for (nm in names(v)) {
        r <- identifiability_reference(v[[nm]], alphaH = ah)
        expect_equal(unique(r$cv_var_ou), 0, tolerance = 1e-8, info = nm)
    }
})

test_that('alpha is derived from alpha*H and the reported height', {
    v <- tree_variant_set(phy, height = 7, randseed = 1)
    r <- identifiability_reference(v$unit_vary, alphaH = c(1, 14))
    expect_equal(r$height, c(7, 7), tolerance = 1e-8)
    expect_equal(r$alpha, c(1 / 7, 2), tolerance = 1e-8)
})

test_that('vcv_inputs reproduces what preprocessTree assembles', {
    n <- Ntip(phy)
    info <- SCOUT:::vcv_inputs(phy)
    rd <- node.depth.edgelength(phy)
    expect_equal(info$n_leaves, n)
    expect_equal(info$leaf_dists, rd[seq_len(n)])
    # shared_lengths[i,j] is the root-to-MRCA distance
    m <- mrca(phy)
    expect_equal(info$shared_lengths[3, 7], rd[m[3, 7]])
    expect_equal(info$shared_lengths[5, 5], rd[5])   # a tip's MRCA with itself is itself
})
