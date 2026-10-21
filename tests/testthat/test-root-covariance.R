root_cov_fixture <- function(unequal = FALSE, fixed = FALSE) {
    phy <- ape::read.tree(text = if (unequal)
        '((a:0.5,b:1):1,(c:1,d:2):0.5);' else
        '((a:1,b:1):1,(c:1,d:1):1);')
    depths <- ape::node.depth.edgelength(phy)
    states <- setNames(rep('1', 7), 1:7)
    info <- list(n_leaves = 4L, leaf_dists = depths[1:4],
        shared_lengths = matrix(depths[ape::mrca(phy)], 4, 4),
        regime_paths = build_regime_paths(phy, states, depths),
        unique_regimes = if (fixed) c('root', '1') else '1',
        root.state = if (fixed) 'root' else '1')
    list(phy = phy, info = info,
         map = getMapFromNode(phy, factor(rep(1, 4)), factor(rep(1, 3)), 0.5))
}

root_cov_loglik <- function(x, mu, covariance) {
    L <- chol(covariance)
    residual <- as.numeric(x - mu)
    -0.5 * (length(x) * log(2 * pi) + 2 * sum(log(diag(L))) +
        sum(backsolve(L, residual, transpose = TRUE)^2))
}

test_that('three-point covariance and likelihood match dense root models', {
    for (unequal in c(FALSE, TRUE)) for (fixed in c(FALSE, TRUE)) {
        f <- root_cov_fixture(unequal, fixed)
        for (alpha in c(0.1, 0.7, 2)) for (tau in c(0, 0.3)) {
            sigma <- 1.3
            # Input root edges must not leak into the transformed root model.
            f$phy$root.edge <- 7
            paths <- lapply(1:4, function(i) getPathToRoot(f$phy, i))
            tr <- transformPhy.new(f$phy, f$map, alpha, sigma, paths, add.root = fixed)
            expect_equal(tr, transformPhy.new(f$phy, f$map, alpha, sigma, add.root = fixed))
            tips <- tr$tree$edge[, 2] <= 4
            tr$tree$edge.length[tips] <- tr$tree$edge.length[tips] +
                tau^2 / tr$diagWeight[tr$tree$edge[tips, 2]]^2
            D <- diag(tr$diagWeight)
            actual <- D %*% (ape::vcv(tr$tree) + tr$tree$root.edge) %*% D
            expected <- compute_VCV(f$info, alpha, sigma, add.root = fixed)
            diag(expected) <- diag(expected) + tau^2
            expect_equal(unname(actual), expected, tolerance = 1e-8)
            x <- setNames(c(1, 2, 4, 3), f$phy$tip.label)
            mu <- matrix(c(2, 2.5, 3, 2), 4, dimnames = list(f$phy$tip.label, NULL))
            expect_equal(loglik_OU_threepoint(x, list(tr, mu, 4)),
                         root_cov_loglik(x, mu, expected), tolerance = 1e-8)
        }
    }
})

test_that('fixed-root Brownian transformation retains its zero-alpha limit', {
    f <- root_cov_fixture()
    tr <- transformPhy.new(f$phy, f$map, 0, 1.3)
    expect_equal(tr$tree$root.edge, 0)
    expect_equal(unname(tr$diagWeight), rep(1, 4))
    expect_equal(ape::vcv(tr$tree), 1.3 * ape::vcv(f$phy))
    for (alpha in c(0, -1, Inf, NA_real_)) {
        expect_error(transformPhy.new(f$phy, f$map, alpha, 1.3, add.root = FALSE),
                     'Stationary root covariance requires positive finite alpha', fixed = TRUE)
    }
})

test_that('M-step objective and final score use the requested root covariance', {
    skip_if_not_installed('nloptr')
    # Use a local function environment so this test also works when sources are loaded directly.
    fit_step <- m_step_tipfog
    environment(fit_step) <- list2env(list(nloptr = nloptr::nloptr),
                                    parent = environment(m_step_tipfog))
    for (unequal in c(FALSE, TRUE)) for (fixed in c(FALSE, TRUE)) {
        for (skip_tau in c(FALSE, TRUE)) {
            f <- root_cov_fixture(unequal, fixed)
            defaults <- list(model = 'OU1', assume.station = !fixed, add.root = fixed,
                skipTau = skip_tau, weight_type = 'alt', algorithm = 'three.point',
                transf = 'ouwie', parsed_alt_tree = f$info,
                tp_const = list(map = f$map, tippath = NULL), verbose = FALSE,
                alpha_shape = 2, alpha_rate = 1, tau_prior_mean = 0.2, tau_prior_sd = 0.1)
            params <- list(alpha = 0.7, sigma = 1.3,
                           theta = if (fixed) c(2, 2) else 2)
            if (!skip_tau) params$tau <- 0.3
            x <- setNames(c(1, 2, 4, 3), f$phy$tip.label)
            fit <- fit_step(x, x, params, f$phy, NULL, defaults, 0.2, 0.2)
            V <- compute_VCV(f$info, fit$alpha, fit$sigma, add.root = fixed)
            diag(V) <- diag(V) + fit$tau^2
            mu <- compute_W_matrix(f$info, fit$alpha) %*% fit$theta
            ll <- root_cov_loglik(x, mu, V)
            expect_equal(fit$ll_total, ll, tolerance = 1e-8)
            expect_equal(fit$ll_ou, ll, tolerance = 1e-8)
            prior <- 4 * 0.2 * (dgamma(fit$alpha, 2, rate = 1, log = TRUE) +
                dnorm(fit$tau, sd(x) * 0.2, 0.1, log = TRUE))
            expect_equal(fit$objective, ll + prior, tolerance = 1e-8)
        }
    }
})
