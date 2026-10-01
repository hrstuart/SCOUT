## Tree + regime-state simulator that feeds simulate_test_data(sim_lin = FALSE).
##
## A birth-death tree from castor::generate_tree_hbd_reverse, tip regimes from a continuous-time
## Markov chain whose rate matrix is Q = logm(P) for a one-time-unit transition matrix P, and
## internal-node regime labels by ancestral state reconstruction (or the simulated truth).
## Ported from build_tree() in scripts/simulated_masking.R (the 260817 counts LOO sweep and the
## 260811 v2 dropout run): with the same seed and the default arguments the tree and tip states are
## identical to that function's.

#' Rate matrix Q = logm(P) from a transition matrix, with the checks logm() does not do.
#' @noRd
rate_matrix <- function(P, tol = 1e-8) {
    if (!is.matrix(P) || nrow(P) != ncol(P) || nrow(P) < 2) stop('`P` must be a square matrix (k >= 2).')
    if (any(!is.finite(P)) || any(P < 0) || any(abs(rowSums(P) - 1) > 1e-6)) {
        stop('`P` must be a transition matrix: non-negative, rows summing to 1.')
    }
    # logm() warns (NaN) or errors when no real logarithm exists; both mean "not embeddable".
    Q <- tryCatch(suppressWarnings(expm::logm(P)), error = function(e) matrix(NaN, nrow(P), ncol(P)))
    if (is.complex(Q)) Q <- Re(Q)
    off <- Q[row(Q) != col(Q)]
    if (any(!is.finite(Q)) || any(off < -tol) || any(abs(rowSums(Q)) > 1e-6)) {
        stop('`P` has no valid rate matrix: logm(P) has negative off-diagonal rates (P is not ',
             'embeddable in a continuous-time chain). Increase the diagonal (stickier states).')
    }
    Q[row(Q) != col(Q)] <- pmax(off, 0)
    diag(Q) <- 0
    diag(Q) <- -rowSums(Q)
    Q
}

#' Simulate a lineage tree with regime states for OU counts simulation
#'
#' Takes a regime transition matrix and a number of cells. Generates an ultrametric birth-death
#' tree of \code{ncells} tips with \code{castor::generate_tree_hbd_reverse}, tip regime states from
#' a continuous-time Markov chain (\code{castor::simulate_mk_model}) with rate matrix
#' \code{Q = expm::logm(P)}, and internal-node regime labels. The result can be passed straight to
#' \code{simulate_test_data(tree = ..., ncells = ncells, sim_lin = FALSE)} to simulate BM1 / OU1 /
#' OUM genes and counts on it.
#'
#' With the 3-state matrix in the example and the same \code{seed}, the tree and tip states are
#' identical to \code{build_tree()} in \code{scripts/simulated_masking.R} (the v2 dropout and
#' counts LOO runs).
#'
#' REGIME CODES. The OU simulator (\code{OUwie.sim.edited}) codes node and tip regimes with two
#' independent \code{factor()} calls and assigns thetas by position, so a regime that appears at the
#' tips but on no internal node silently shifts every later regime's theta -- unless it sorts last.
#' Each draw is therefore accepted only if the internal-node state set is a prefix of the sorted
#' tip state set; otherwise the states (not the tree) are redrawn, up to \code{max_tries}. With
#' \code{require_all_states = TRUE} every state must also appear at the tips. Redraws condition the
#' state distribution on these events; the number taken is in \code{attr(tree, 'state_draws')}.
#'
#' @param P k x k regime transition matrix over one time unit (rows sum to 1); state i is row i.
#'   Must be embeddable (\code{logm(P)} a valid rate matrix), which holds for diagonally dominant
#'   (sticky) matrices.
#' @param ncells number of cells (tips).
#' @param birth_rate,death_rate,rho speciation rate, extinction rate and sampling fraction for
#'   \code{generate_tree_hbd_reverse} (\code{lambda}, \code{mu}, \code{rho}).
#' @param state_names names of the k states; default \code{rownames(P)} if set, else
#'   \code{state1..statek}.
#' @param tip_prefix prefix for tip labels (\code{t1, t2, ...}).
#' @param node_states \code{'inferred'} (default): internal nodes labelled by ancestral state
#'   reconstruction from the tips (\code{infer_anc()}), as when fitting real data;
#'   \code{'true'}: the states the Markov chain actually visited at internal nodes.
#' @param anc_infer method for \code{infer_anc()}: \code{'ape'} or \code{'castor'}.
#' @param require_all_states require every state to appear at the tips.
#' @param max_tries maximum number of state draws.
#' @param seed optional RNG seed (\code{set.seed}) for the whole simulation.
#' @param verbose log tree and state summaries.
#' @return A \code{phylo} with \code{$states} (tip states, tip.label order) and
#'   \code{$node.label} (internal-node states), and attributes \code{Q}, \code{state_draws},
#'   \code{node_states}.
#' @examples
#' \dontrun{
#' P <- matrix(c(0.90, 0.08, 0.02,
#'               0.10, 0.80, 0.10,
#'               0.02, 0.08, 0.90), 3, 3, byrow = TRUE)
#' phy <- simulate_tree_states(P, ncells = 256, seed = 1)
#' table(phy$states)
#' sim <- simulate_test_data(ngenes = 30, ncells = 256, tree = phy, outdir = 'sim',
#'                           out_prefix = 'hbd256', a = 3, s = 1, t0 = 5, theta_step = 2,
#'                           sim_lin = FALSE)
#' }
#' @export
simulate_tree_states <- function(P,
                                 ncells,
                                 birth_rate = 1,
                                 death_rate = 0,
                                 rho = 1,
                                 state_names = NULL,
                                 tip_prefix = 't',
                                 node_states = c('inferred', 'true'),
                                 anc_infer = 'ape',
                                 require_all_states = TRUE,
                                 max_tries = 100,
                                 seed = NULL,
                                 verbose = TRUE) {
    node_states <- match.arg(node_states)
    if (!is.null(seed)) set.seed(seed)
    Q <- rate_matrix(P)
    k <- nrow(Q)
    if (is.null(state_names)) state_names <- rownames(P) %||% paste0('state', seq_len(k))
    if (length(state_names) != k || anyDuplicated(state_names)) {
        stop(sprintf('`state_names` must be %d distinct names.', k))
    }

    obj <- castor::generate_tree_hbd_reverse(ncells, rho = rho, lambda = birth_rate, mu = death_rate)
    if (!isTRUE(obj$success) || !length(obj$trees)) {
        stop('generate_tree_hbd_reverse failed: ', obj$error %||% 'no tree returned')
    }
    phy <- obj$trees[[1]]
    phy$tip.label <- paste0(tip_prefix, phy$tip.label)

    for (draw in seq_len(max_tries)) {
        st <- castor::simulate_mk_model(phy, Q = Q, include_nodes = node_states == 'true')
        tips <- state_names[st$tip_states]
        cand <- phy
        cand$states <- tips
        if (node_states == 'true') {
            cand$node.label <- state_names[st$node_states]
        } else {
            cand <- infer_anc(cand, anc_infer)
        }
        tip_lv  <- sort(unique(tips))
        node_lv <- sort(unique(cand$node.label))
        all_ok    <- !require_all_states || length(tip_lv) == k
        prefix_ok <- length(node_lv) <= length(tip_lv) &&
                     identical(node_lv, tip_lv[seq_along(node_lv)])
        if (all_ok && prefix_ok) break
        if (draw == max_tries) {
            stop(sprintf(paste0('No acceptable state draw in %d tries (tip states: %s; node states: %s). ',
                                'Rare states are vanishing: raise ncells, adjust P, set ',
                                'require_all_states = FALSE, or name the rarest state to sort last.'),
                         max_tries, paste(tip_lv, collapse = ' '), paste(node_lv, collapse = ' ')))
        }
    }

    # Same normalisation runSCOUT.dropout applies, so the simulated tree is the tree that is fitted.
    phy <- prep_dropout_tree(cand, states = stats::setNames(cand$states, cand$tip.label),
                             verbose = FALSE, infer_nodes = FALSE)
    attr(phy, 'Q') <- Q
    attr(phy, 'state_draws') <- draw
    attr(phy, 'node_states') <- node_states
    if (verbose) {
        tab <- table(factor(phy$states, levels = state_names))
        log_message(sprintf('simulate_tree_states: %d tips, height %.3f, %s; tip states %s (%d draw%s)',
                            ncells, max(ape::node.depth.edgelength(phy)), node_states,
                            paste(names(tab), tab, sep = '=', collapse = ' '), draw,
                            if (draw == 1) '' else 's'), verbose = TRUE)
    }
    phy
}
