test_that('root weights use regime labels regardless of factor codes or column order', {
    # All internal nodes are B, so its factor code is 1, while the first column is A.
    phy <- ape::read.tree(text = '((a:1,b:1):1,(c:1,d:1):1);')
    root_dists <- ape::node.depth.edgelength(phy)
    node_regimes <- c('A', 'A', 'B', 'B', 'B', 'B', 'B')
    names(node_regimes) <- seq_along(node_regimes)
    info <- list(
        n_leaves = 4L,
        unique_regimes = c('A', 'B'),
        leaf_dists = root_dists[1:4],
        regime_paths = build_regime_paths(phy, node_regimes, root_dists),
        root.state = factor('B')
    )

    weights <- compute_W_matrix(info, alpha = 1)
    # A occupies the terminal branches (length 1) for tips a and b.
    expect_equal(unname(weights[, 'A']), c(1 - exp(-1), 1 - exp(-1), 0, 0))
    expect_equal(unname(weights[, 'B']), c(exp(-1), exp(-1), 1, 1))
    expect_equal(unname(rowSums(weights)), rep(1, 4))

    info$root.state <- 'B'
    expect_equal(compute_W_matrix(info, alpha = 1), weights)
    info$root.state <- factor('B')
    info$unique_regimes <- c('B', 'A')
    expect_equal(compute_W_matrix(info, alpha = 1), weights[, c('B', 'A')])
})

test_that('root regime must identify exactly one weights-matrix column', {
    info <- list(n_leaves = 1L, unique_regimes = c('A', 'B'))
    for (root in list(NULL, NA_character_, 'missing', c('A', 'B'))) {
        info$root.state <- root
        expect_error(compute_W_matrix(info, alpha = 1),
                     'Root regime must match exactly one weights-matrix column', fixed = TRUE)
    }
    info$root.state <- 'B'
    info$unique_regimes <- c('B', 'B')
    expect_error(compute_W_matrix(info, alpha = 1),
                 'Root regime must match exactly one weights-matrix column', fixed = TRUE)
})
