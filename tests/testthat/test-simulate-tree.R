## Tests for simulate_tree_states(): castor hbd tree + Mk regime states for OU counts simulation.

suppressPackageStartupMessages({
    library(ape)
    library(castor)
})

P3 <- matrix(c(0.90, 0.08, 0.02,
               0.10, 0.80, 0.10,
               0.02, 0.08, 0.90), 3, 3, byrow = TRUE)

## build_tree() from scripts/simulated_masking.R (cfg$tree_seed -> seed)
build_tree_ref <- function(n, seed) {
    set.seed(seed)
    Q <- expm::logm(P3)
    obj <- generate_tree_hbd_reverse(n, rho = 1, lambda = 1, mu = 0)
    phy <- obj$trees[[1]]
    phy$tip.label <- paste0('t', phy$tip.label)
    st  <- simulate_mk_model(phy, Q = Q, include_nodes = FALSE)
    phy$states <- as.factor(paste0('state', st$tip_states))
    phy <- infer_anc(phy)
    prep_dropout_tree(phy, states = setNames(as.character(phy$states), phy$tip.label), verbose = FALSE)
}

test_that('reproduces build_tree() from the dropout scripts for the same seed', {
    for (sd in 1:2) {
        a <- build_tree_ref(96, sd)
        b <- simulate_tree_states(P3, ncells = 96, seed = sd, verbose = FALSE)
        expect_identical(b$edge, a$edge)
        expect_equal(b$edge.length, a$edge.length)
        expect_identical(b$tip.label, a$tip.label)
        expect_identical(b$states, as.character(a$states))
        expect_identical(b$node.label, a$node.label)
    }
})

test_that('output is an ultrametric binary tree with valid regime labels', {
    phy <- simulate_tree_states(P3, ncells = 80, seed = 4, verbose = FALSE)
    expect_s3_class(phy, 'phylo')
    expect_equal(Ntip(phy), 80)
    expect_true(is.ultrametric(phy)); expect_true(is.binary(phy))
    expect_length(phy$states, 80)
    expect_length(phy$node.label, phy$Nnode)
    expect_setequal(unique(phy$states), paste0('state', 1:3))       # require_all_states
    node_lv <- sort(unique(phy$node.label)); tip_lv <- sort(unique(phy$states))
    expect_identical(node_lv, tip_lv[seq_along(node_lv)])           # prefix: no code shift
    expect_equal(rowSums(attr(phy, 'Q')), rep(0, 3), tolerance = 1e-10)
})

test_that('state names come from rownames(P) or state_names; true node states are available', {
    Pn <- P3; rownames(Pn) <- colnames(Pn) <- c('A', 'B', 'C')
    phy <- simulate_tree_states(Pn, ncells = 60, seed = 2, verbose = FALSE)
    expect_true(all(phy$states %in% c('A', 'B', 'C')))
    phy2 <- simulate_tree_states(P3, ncells = 60, seed = 2, state_names = c('x', 'y', 'z'),
                                 node_states = 'true', verbose = FALSE)
    expect_true(all(c(phy2$states, phy2$node.label) %in% c('x', 'y', 'z')))
    expect_equal(attr(phy2, 'node_states'), 'true')
    expect_error(simulate_tree_states(P3, 20, state_names = c('a', 'a', 'b')), 'distinct')
})

test_that('invalid transition matrices are rejected', {
    expect_error(simulate_tree_states(matrix(c(0.5, 0.6, 0.5, 0.4), 2, byrow = TRUE), 20), 'rows summing to 1')
    flip <- matrix(c(0.1, 0.9, 0.9, 0.1), 2, 2)                    # not embeddable
    expect_error(simulate_tree_states(flip, 20, verbose = FALSE), 'embeddable')
    expect_error(simulate_tree_states(P3[1:2, ], 20), 'square')
})

test_that('an impossible state requirement fails after max_tries', {
    expect_error(suppressWarnings(simulate_tree_states(P3, ncells = 2, max_tries = 3, seed = 1, verbose = FALSE)),
                 'No acceptable state draw')
    phy <- suppressWarnings(  # ape::ace on 2 tips cannot estimate rate SEs
        simulate_tree_states(P3, ncells = 2, require_all_states = FALSE, seed = 1, verbose = FALSE))
    expect_equal(Ntip(phy), 2)
})

test_that('the tree feeds simulate_test_data(sim_lin = FALSE)', {
    phy <- simulate_tree_states(P3, ncells = 64, seed = 3, verbose = FALSE)
    od <- tempfile(); dir.create(od); on.exit(unlink(od, recursive = TRUE))
    sim <- suppressWarnings(simulate_test_data(ngenes = 3, ncells = 64, tree = phy, outdir = od,
                                               out_prefix = 't', a = 3, s = 1, t0 = 5, theta_step = 2,
                                               nevfs = 3, sim_lin = FALSE))
    cnt <- sim$simulated_OU$counts[[1]]
    expect_equal(nrow(cnt), 64)
    expect_identical(as.character(cnt$OUM), phy$states[match(cnt$species, phy$tip.label)])
})
