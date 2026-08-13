## Regression tests for the held-out dropout validation machinery.
## These exercise predict_masked_latent against e_step, against brute-force linear algebra,
## and against data drawn from SCOUT's own generative model.

suppressPackageStartupMessages({
    library(ape)
    library(castor)
})

## ---- shared fixture -------------------------------------------------------------------

make_fixture <- function(n = 40, seed = 11) {
    set.seed(seed)
    P <- matrix(c(0.90, 0.08, 0.02,
                  0.10, 0.80, 0.10,
                  0.02, 0.08, 0.90), 3, 3, byrow = TRUE)
    phy <- generate_tree_hbd_reverse(n, rho = 1, lambda = 1, mu = 0)$trees[[1]]
    phy$tip.label <- paste0('t', seq_len(n))
    st <- simulate_mk_model(phy, Q = expm::logm(P), include_nodes = FALSE)
    states <- setNames(paste0('state', st$tip_states), phy$tip.label)
    phy <- SCOUT:::prep_dropout_tree(phy, states, verbose = FALSE)

    meta <- data.frame(species = phy$tip.label, OUM = phy$states,
                       g1 = abs(rnorm(n, 8, 1)), g2 = abs(rnorm(n, 8, 1)),
                       stringsAsFactors = FALSE)
    idata <- suppressMessages(
        formatSCOUT(phy, meta, tempdir(), species_key = 'species',
                    regimes = c('BM1', 'OU1', 'OUM'), anc_infer = 'ape', normalize = FALSE))
    list(phy = phy, states = setNames(phy$states, phy$tip.label), idata = idata)
}

setup_regime <- function(fx, reg) {
    ti <- SCOUT:::preprocessTree(fx$idata, reg = reg, scaleHeight = FALSE,
                                 root.fixed = FALSE, skipTau = FALSE)
    gi <- SCOUT:::preprocessGene(fx$idata, 'g1', ti$tree_const, ti$defaults$root.state,
                                 ti$defaults$add.root, ti$defaults$model, ti$defaults$skipTau)
    gi$param_init <- SCOUT:::verify_initial_params(gi$param_init, ti$defaults, FALSE)
    tips <- ti$tree_const$tree$tip.label
    list(ti = ti, gi = gi, tips = tips, X = gi$X[tips])
}

fx <- make_fixture()

## ---- the estimator reduces to the E-step ----------------------------------------------

test_that('predict_masked_latent reduces to e_step when conditioning on every tip', {
    for (reg in c('BM1', 'OU1', 'OUM')) {
        s <- setup_regime(fx, reg)
        d <- SCOUT:::align_regimes(s$ti$defaults, s$gi$param_init$theta)
        d$verbose <- FALSE

        er <- SCOUT:::e_step(s$X, s$ti$tree_const$tree, s$ti$tree_const$edges,
                             s$gi$param_init, d, diagnose = FALSE)
        pr <- SCOUT:::predict_masked_latent(s$gi$param_init, s$ti$defaults,
                                            held_idx = seq_along(s$tips), X_full = s$X,
                                            kept_idx = seq_along(s$tips))

        expect_equal(as.vector(er$E_Z), unname(pr$mean), tolerance = 1e-8,
                     info = paste('regime', reg))
        expect_equal(diag(er$Var_Z_given_X), unname(pr$var), tolerance = 1e-8,
                     info = paste('regime', reg))
    }
})

test_that('a single held-out tip matches a brute-force conditional', {
    s <- setup_regime(fx, 'OUM')
    d <- SCOUT:::align_regimes(s$ti$defaults, s$gi$param_init$theta)
    p <- s$gi$param_init

    V <- SCOUT:::compute_VCV(d$parsed_alt_tree, p$alpha, p$sigma, add.root = d$add.root)
    W <- SCOUT:::compute_W_matrix(d$parsed_alt_tree, p$alpha, add.root = d$add.root)
    mu <- as.vector(W %*% p$theta)

    h <- 3L
    o <- setdiff(seq_along(s$tips), h)
    Voo <- V[o, o] + diag(p$tau^2, length(o)) + 1e-8 * diag(length(o))
    brute <- mu[h] + V[h, o, drop = FALSE] %*% solve(Voo, s$X[o] - mu[o])

    pr <- SCOUT:::predict_masked_latent(p, s$ti$defaults, held_idx = h, X_full = s$X)
    expect_equal(as.numeric(brute), unname(pr$mean), tolerance = 1e-8)
})

## ---- held-out means held out ----------------------------------------------------------

test_that('predictions do not depend on the masked observations', {
    s <- setup_regime(fx, 'OUM')
    cl <- SCOUT:::sample_clade(fx$phy, 4, 10, states = fx$states)
    h <- match(cl$tips, s$tips)

    X2 <- s$X
    X2[h] <- X2[h] + 1000

    a <- SCOUT:::predict_masked_latent(s$gi$param_init, s$ti$defaults, h, s$X)
    b <- SCOUT:::predict_masked_latent(s$gi$param_init, s$ti$defaults, h, X2)
    expect_equal(a$mean, b$mean, tolerance = 1e-10)
})

## ---- ordering: this is the bug that motivated the test file ---------------------------

test_that('output rows follow the order of held_idx, not sorted order', {
    s <- setup_regime(fx, 'OUM')
    h <- match(SCOUT:::sample_random_tips(fx$phy, 8, states = fx$states), s$tips)
    perm <- rev(h)

    a <- SCOUT:::predict_masked_latent(s$gi$param_init, s$ti$defaults, h, s$X)
    b <- SCOUT:::predict_masked_latent(s$gi$param_init, s$ti$defaults, perm, s$X)

    expect_identical(names(a$mean), s$tips[h])
    expect_identical(names(b$mean), s$tips[perm])
    expect_equal(unname(rev(b$mean)), unname(a$mean), tolerance = 1e-12)
})

test_that('duplicate held-out tips are rejected', {
    s <- setup_regime(fx, 'OUM')
    expect_error(SCOUT:::predict_masked_latent(s$gi$param_init, s$ti$defaults,
                                               c(1L, 2L, 2L), s$X), 'duplicate')
})

## ---- statistical correctness under the true model --------------------------------------

test_that('with true parameters the conditional beats the mean baseline and is calibrated', {
    s <- setup_regime(fx, 'OUM')
    d <- SCOUT:::align_regimes(s$ti$defaults, s$gi$param_init$theta)
    n <- length(s$tips)

    p <- list(tau = 0.1, alpha = 0.3, sigma = 1, theta = s$gi$param_init$theta)
    V <- SCOUT:::compute_VCV(d$parsed_alt_tree, p$alpha, p$sigma, add.root = d$add.root)
    W <- SCOUT:::compute_W_matrix(d$parsed_alt_tree, p$alpha, add.root = d$add.root)
    mu <- as.vector(W %*% p$theta)
    L <- chol(V + 1e-8 * diag(n))

    set.seed(5)
    acc <- t(replicate(40, {
        Z <- as.vector(mu + t(L) %*% rnorm(n))
        X <- setNames(Z + rnorm(n, 0, p$tau), s$tips)
        h <- match(SCOUT:::sample_random_tips(fx$phy, 8, states = fx$states), s$tips)
        pr <- SCOUT:::predict_masked_latent(p, s$ti$defaults, h, X)
        c(masked = sqrt(mean((pr$mean - Z[h])^2)),
          global = sqrt(mean((mean(X[pr$kept_idx]) - Z[h])^2)),
          cov = mean(abs(pr$mean - Z[h]) <= 1.96 * sqrt(pr$var)))
    }))

    expect_lt(mean(acc[, 'masked']), 0.85 * mean(acc[, 'global']))
    expect_equal(mean(acc[, 'cov']), 0.95, tolerance = 0.07)
})

## ---- masking arms ----------------------------------------------------------------------

test_that('a masked clade gives a rank-one cross-covariance but scattered tips do not', {
    s <- setup_regime(fx, 'OUM')
    d <- SCOUT:::align_regimes(s$ti$defaults, s$gi$param_init$theta)
    V <- SCOUT:::compute_VCV(d$parsed_alt_tree, s$gi$param_init$alpha,
                             s$gi$param_init$sigma, add.root = d$add.root)

    cl <- SCOUT:::sample_clade(fx$phy, 5, 12, states = fx$states)
    h <- match(cl$tips, s$tips)
    expect_equal(qr(V[h, setdiff(seq_along(s$tips), h), drop = FALSE])$rank, 1L)

    hr <- match(SCOUT:::sample_random_tips(fx$phy, length(cl$tips), states = fx$states), s$tips)
    expect_gt(qr(V[hr, setdiff(seq_along(s$tips), hr), drop = FALSE])$rank, 1L)
})

test_that('sample_clade respects size bounds and keeps every regime represented', {
    for (i in 1:15) {
        cl <- SCOUT:::sample_clade(fx$phy, 4, 10, states = fx$states)
        if (is.null(cl)) next
        expect_gte(length(cl$tips), 4)
        expect_lte(length(cl$tips), 10)
        kept <- setdiff(fx$phy$tip.label, cl$tips)
        expect_setequal(unique(fx$states[kept]), unique(fx$states))
        # never one of the root's children, which would collapse the root
        root <- length(fx$phy$tip.label) + 1
        expect_false(cl$node %in% fx$phy$edge[fx$phy$edge[, 1] == root, 2])
    }
})

test_that('descendant_tips agrees with ape::extract.clade', {
    desc <- SCOUT:::descendant_tips(fx$phy)
    n <- length(fx$phy$tip.label)
    for (nd in sample((n + 2):(n + fx$phy$Nnode), 8)) {
        expect_setequal(fx$phy$tip.label[desc[[nd]]], extract.clade(fx$phy, nd)$tip.label)
    }
})

## ---- guards ----------------------------------------------------------------------------

test_that('compute_W_matrix errors on a regime with no matching column', {
    s <- setup_regime(fx, 'OUM')
    bad <- s$ti$defaults$parsed_alt_tree
    bad$unique_regimes <- c('state1', 'not_a_regime')
    expect_error(SCOUT:::compute_W_matrix(bad, 0.5, add.root = FALSE), 'matches 0 columns')
})

test_that('check_regime_coverage flags a regime the fit does not cover', {
    s <- setup_regime(fx, 'OUM')
    d <- SCOUT:::align_regimes(s$ti$defaults, s$gi$param_init$theta)
    th <- s$gi$param_init$theta
    expect_length(SCOUT:::check_regime_coverage(d, th), 0)
    expect_identical(SCOUT:::check_regime_coverage(d, th[-1]), names(th)[1])
})

## ---- the simulator needs a preordered tree ---------------------------------------------

test_that('OUwie.sim.edited propagates the root value on a non-preordered tree', {
    # castor builds trees whose edge matrix is not in preorder; before the reorder inside
    # OUwie.sim.edited that silently discarded theta0 and the phylogenetic covariance.
    preordered <- function(p) {
        seen <- rep(FALSE, max(p$edge))
        seen[length(p$tip.label) + 1] <- TRUE
        for (i in seq_len(nrow(p$edge))) {
            if (!seen[p$edge[i, 1]]) return(FALSE)
            seen[p$edge[i, 2]] <- TRUE
        }
        TRUE
    }
    expect_false(preordered(fx$phy))   # fixture really is the awkward case

    meta <- data.frame(cellID = fx$phy$tip.label, cluster = as.factor(fx$phy$states))
    set.seed(1)
    v <- replicate(25, SCOUT:::OUwie.sim.edited(fx$phy, meta, alpha = rep(1e-10, 3),
                                                sigma.sq = rep(1, 3), theta0 = 10,
                                                theta = rep(10, 3), simmap.tree = FALSE)$X)
    # a BM starting at theta0 = 10 must stay centred there, not collapse to 0
    expect_equal(mean(v), 10, tolerance = 1)
})

## ---- caller-controlled mask selection --------------------------------------------------

test_that('eligible_clades enumerates exactly what sample_clade draws from', {
    cand <- eligible_clades(fx$phy, 4, 10, states = fx$states)
    expect_gt(nrow(cand), 0)
    expect_true(all(cand$size >= 4 & cand$size <= 10))
    expect_identical(cand$size, lengths(attr(cand, 'tips')))

    # root children are excluded, so masking never collapses the root
    root <- length(fx$phy$tip.label) + 1
    expect_false(any(cand$node %in% fx$phy$edge[fx$phy$edge[, 1] == root, 2]))

    # every regime survives removal of any candidate
    all_reg <- unique(as.character(fx$states))
    for (i in seq_len(nrow(cand))) {
        kept <- setdiff(fx$phy$tip.label, attr(cand, 'tips')[[i]])
        expect_setequal(unique(as.character(fx$states[kept])), all_reg)
    }

    for (i in 1:20) {
        cl <- SCOUT:::sample_clade(fx$phy, 4, 10, states = fx$states)
        expect_true(cl$node %in% cand$node)
        expect_setequal(cl$tips, attr(cand, 'tips')[[match(cl$node, cand$node)]])
    }
})

test_that('clade_tips agrees with extract.clade and rejects bad nodes', {
    cand <- eligible_clades(fx$phy, 4, 10, states = fx$states)
    nd <- cand$node[1]
    expect_setequal(SCOUT:::clade_tips(fx$phy, nd), extract.clade(fx$phy, nd)$tip.label)
    expect_error(SCOUT:::clade_tips(fx$phy, 1e6), 'not a node')
})

test_that('dropout_states builds a composite that is stricter than either column', {
    meta <- data.frame(species = c('a', 'b', 'c'), r1 = c('X', 'X', 'Y'),
                       r2 = c('P', 'Q', 'Q'), stringsAsFactors = FALSE)
    st <- dropout_states(meta, c('r1', 'r2'))
    expect_identical(unname(st), c('X.P', 'X.Q', 'Y.Q'))
    expect_identical(names(st), c('a', 'b', 'c'))
    # 3 composite cells vs 2 per column: keeping every composite cell keeps every column value
    expect_gt(length(unique(st)), length(unique(meta$r1)))
    expect_error(dropout_states(meta, c('r1', 'nope')), 'not found')
})

test_that('prep_dropout_tree can skip ancestral state inference', {
    set.seed(3)
    phy <- generate_tree_hbd_reverse(20, rho = 1, lambda = 1, mu = 0)$trees[[1]]
    phy$tip.label <- paste0('t', seq_len(20))
    phy$node.label <- NULL
    st <- setNames(rep(c('a', 'b'), each = 10), phy$tip.label)
    expect_null(SCOUT:::prep_dropout_tree(phy, st, verbose = FALSE, infer_nodes = FALSE)$node.label)
    expect_false(is.null(SCOUT:::prep_dropout_tree(phy, st, verbose = FALSE)$node.label))
})

test_that('score_recovery reports wider coverage on the observed scale', {
    set.seed(5)
    truth <- rnorm(200)
    pred <- truth + rnorm(200, sd = 0.5)
    sd_lat <- rep(0.2, 200)
    sc <- SCOUT:::score_recovery(truth, pred, pred_sd = sd_lat,
                                 pred_sd_obs = sqrt(sd_lat^2 + 0.5^2))
    expect_gt(sc$coverage_obs, sc$coverage)
    expect_true(is.na(SCOUT:::score_recovery(truth, pred)$coverage_obs))
})

test_that('runSCOUT.dropout validates random_frac before doing any work', {
    fx <- make_fixture(n = 40, seed = 11)
    d <- tempfile()

    # Entries must be NA (size-matched) or a proper fraction.
    expect_error(runSCOUT.dropout(fx$phy, d, random_frac = 1.5), 'must be NA')
    expect_error(runSCOUT.dropout(fx$phy, d, random_frac = 0), 'must be NA')
    expect_error(runSCOUT.dropout(fx$phy, d, random_frac = -0.1), 'must be NA')
    # Duplicates would collide as arm names.
    expect_error(runSCOUT.dropout(fx$phy, d, random_frac = c(0.1, 0.1)), 'duplicate')
    expect_error(runSCOUT.dropout(fx$phy, d, random_frac = numeric(0)), 'at least one entry')

    # The largest fraction must leave enough tips to fit on: 0.9 * 40 = 36 masked of 40.
    expect_error(runSCOUT.dropout(fx$phy, d, clade_min = 4, clade_max = 10, random_frac = 0.9),
                 'too few')
    # NULL is the documented default and must stay a legal value.
    expect_null(formals(runSCOUT.dropout)$random_frac)
})
