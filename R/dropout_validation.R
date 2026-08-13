
####### HELD-OUT DROPOUT VALIDATION #######
#
# SCOUT models tip expression as X = Z + eps, where Z ~ N(W theta, V) is an OU process on the
# lineage tree and eps ~ N(0, tau^2 I) is tip fog. These functions mask part of the tree, fit
# SCOUT-EM without it, and then predict the masked tips' latent expression from the retained
# tips alone. If that prediction beats a no-phylogeny baseline, the OU-on-tree covariance is
# carrying real information across the tree.
#
# e_step() computes the posterior over ALL tips jointly:
#     E_Z = W theta + V (V + tau^2 I)^-1 (X - W theta)
# Partitioning tips into retained `o` and held-out `h`, and using Cov(Z_h, X_o) = V_ho and
# Var(X_o) = V_oo + tau^2 I, the same Gaussian algebra gives
#     E[Z_h | X_o]   = W_h theta + V_ho (V_oo + tau^2 I)^-1 (X_o - W_o theta)
#     Var[Z_h | X_o] = V_hh      - V_ho (V_oo + tau^2 I)^-1 V_oh
# e_step is the special case h = o = all tips, which is what the reduction test exercises.


`%||%` <- function(x, y) if (is.null(x)) y else x


######## TREE MASKING ########

#' Descendant tip indices for every node of a tree
#' @param phy Tree of class phylo.
#' @return List indexed by node number. Entries 1:Ntip are the tips themselves.
descendant_tips <- function(phy) {
    n_tips <- length(phy$tip.label)
    n_all <- n_tips + phy$Nnode

    desc <- vector('list', n_all)
    for (i in seq_len(n_tips)) desc[[i]] <- i

    # postorder guarantees every child is visited before its parent.
    ed <- reorder.phylo(phy, 'postorder')$edge
    for (e in seq_len(nrow(ed))) {
        desc[[ed[e, 1]]] <- c(desc[[ed[e, 1]]], desc[[ed[e, 2]]])
    }

    return(desc)
}

#' Every clade of a tree that is eligible to be masked
#'
#' Enumerates the candidate set that `sample_clade` draws from. Exposed separately so a caller
#' can inspect what the tree actually offers -- and pick a deliberate subset via
#' `runSCOUT.dropout(mask_nodes = )` -- before committing to a long run.
#'
#' @param phy Tree of class phylo.
#' @param min_size Minimum number of tips in the clade.
#' @param max_size Maximum number of tips in the clade.
#' @param states Named vector of tip regimes. Clades whose removal would eliminate a regime
#'   from the retained tips are rejected, because a regime with no retained tips has no
#'   estimable theta and would silently corrupt the weights matrix.
#' @param exclude_root_children Whether to reject the two children of the root. Dropping one
#'   of them collapses the root and shifts every root-to-tip distance, so the masked-tree
#'   parameters would no longer be on the same time scale as the full tree.
#' @return data.frame with `node` and `size`, in increasing node order, carrying the tip labels
#'   as the `tips` attribute (a list aligned to the rows).
#' @export
eligible_clades <- function(phy, min_size, max_size, states = NULL, exclude_root_children = TRUE) {
    n_tips <- length(phy$tip.label)
    root <- n_tips + 1
    desc <- descendant_tips(phy)

    internal <- setdiff((n_tips + 1):(n_tips + phy$Nnode), root)
    if (exclude_root_children) {
        internal <- setdiff(internal, phy$edge[phy$edge[, 1] == root, 2])
    }

    all_regimes <- if (is.null(states)) NULL else unique(as.character(states[phy$tip.label]))

    nodes <- integer(0)
    tip_sets <- list()
    for (nd in internal) {
        tips <- phy$tip.label[desc[[nd]]]
        if (length(tips) < min_size || length(tips) > max_size) next
        if (!is.null(states)) {
            kept <- setdiff(phy$tip.label, tips)
            if (!setequal(unique(as.character(states[kept])), all_regimes)) next
        }
        nodes <- c(nodes, nd)
        tip_sets[[length(tip_sets) + 1]] <- tips
    }

    out <- data.frame(node = nodes, size = lengths(tip_sets))
    attr(out, 'tips') <- tip_sets
    return(out)
}

#' Tip labels below a given node
#' @param phy Tree of class phylo.
#' @param node Node number.
#' @return Character vector of tip labels.
clade_tips <- function(phy, node) {
    n_all <- length(phy$tip.label) + phy$Nnode
    if (length(node) != 1 || is.na(node) || node < 1 || node > n_all) {
        stop(sprintf('Node %s is not a node of this tree (1..%d).', as.character(node), n_all))
    }
    return(phy$tip.label[descendant_tips(phy)[[node]]])
}

#' Sample a random clade to mask
#' @inheritParams eligible_clades
#' @return List with `node` and `tips`, or NULL when no clade satisfies the constraints.
sample_clade <- function(phy, min_size, max_size, states = NULL, exclude_root_children = TRUE) {
    cand <- eligible_clades(phy, min_size, max_size, states, exclude_root_children)
    if (nrow(cand) == 0) {
        return(NULL)
    }

    i <- sample(nrow(cand), 1)
    return(list(node = cand$node[i], tips = attr(cand, 'tips')[[i]]))
}

#' Composite tip states for mask eligibility
#'
#' Pastes several regime columns into one state per tip. Requiring every *composite* cell to keep
#' at least one retained tip is stricter than requiring it of each column separately, which is the
#' safe direction: no individual regime can be wiped out by a mask that the composite accepts.
#'
#' @param meta Metadata data.frame.
#' @param cols Column names to combine.
#' @param species_key Column holding the tip labels.
#' @return Character vector of states named by tip label.
#' @export
dropout_states <- function(meta, cols, species_key = 'species') {
    missing <- setdiff(c(cols, species_key), colnames(meta))
    if (length(missing) > 0) {
        stop(sprintf('Columns not found in the metadata: %s', paste(missing, collapse = ', ')))
    }
    # paste, not interaction(): interaction() returns a factor whose level set is the full cross
    # product, and every consumer here works on values rather than levels.
    s <- do.call(paste, c(lapply(cols, function(cc) as.character(meta[[cc]])), sep = '.'))
    return(setNames(s, as.character(meta[[species_key]])))
}

#' Sample k tips uniformly at random, keeping every regime represented
#' @param phy Tree of class phylo.
#' @param k Number of tips to mask.
#' @param states Named vector of tip regimes, used as in sample_clade.
#' @param max_tries How many times to resample before giving up.
#' @return Character vector of tip labels, or NULL.
sample_random_tips <- function(phy, k, states = NULL, max_tries = 50) {
    all_regimes <- if (is.null(states)) NULL else unique(as.character(states[phy$tip.label]))

    for (i in seq_len(max_tries)) {
        tips <- sample(phy$tip.label, k)
        if (is.null(states)) {
            return(tips)
        }
        kept <- setdiff(phy$tip.label, tips)
        if (setequal(unique(as.character(states[kept])), all_regimes)) {
            return(tips)
        }
    }

    return(NULL)
}


######## PREDICTION ########

#' Align the regime bookkeeping in `defaults` to a fitted theta vector
#' Mirrors the patch runSCOUT applies before calling run_em, which is what keeps the columns
#' of the weights matrix lined up with theta.
#' @param defaults Defaults list from preprocessTree.
#' @param theta Named theta vector from a fit.
#' @return The modified defaults list.
align_regimes <- function(defaults, theta) {
    unq <- defaults$parsed_alt_tree$unique_regimes
    if (length(theta) != length(unq) || any(names(theta) != as.character(unq))) {
        defaults$parsed_alt_tree$unique_regimes <- names(theta)
    }
    if (!'root' %in% defaults$parsed_alt_tree$unique_regimes & defaults$add.root) {
        defaults$parsed_alt_tree$unique_regimes <- c('root', defaults$parsed_alt_tree$unique_regimes)
    }

    # compute_W_matrix does `W[, root.state]`, and preprocessTree leaves root.state as an element
    # of the int.states FACTOR, so that subscript resolves by level code rather than by name. It
    # lands on the right column only when the root regime sits at the same index in the factor's
    # levels as in the (re-ordered) regime names. That holds for the datasets in use here, but it
    # is a coincidence, so make a violation loud instead of silently building a wrong W.
    rs <- defaults$parsed_alt_tree$root.state
    if (is.factor(rs)) {
        regs <- as.character(defaults$parsed_alt_tree$unique_regimes)
        by_name <- which(regs == as.character(rs))
        if (length(by_name) != 1 || by_name != as.integer(rs)) {
            stop(sprintf(paste0('Root regime `%s` is column %s by name but %d by factor code. ',
                                'compute_W_matrix would add the root contribution to the wrong ',
                                'column. Regimes: %s'),
                         as.character(rs),
                         if (length(by_name) == 1) as.character(by_name) else 'unmatched',
                         as.integer(rs), paste(regs, collapse = ', ')))
        }
    }

    return(defaults)
}

#' Recover a fitted parameter set with theta's regime names restored
#' nloptr returns an unnamed solution vector, so m_step_tipfog strips the names off theta and
#' every downstream consumer has to put them back from the fit's own regime list. This mirrors
#' what extract_parameters does, including dropping the leading 'root' entry for BM1, where
#' runSCOUT prepends a root regime but theta still carries a single value.
#' @param res One element of a runSCOUT result list.
#' @return The `paras` list with a named theta, or NULL if the fit is unusable.
fit_paras <- function(res) {
    if (is.null(res) || is.null(res$paras) || is.null(res$paras$theta)) {
        return(NULL)
    }
    paras <- res$paras
    theta <- unlist(paras$theta)
    regs <- as.character(res$defaults$parsed_alt_tree$unique_regimes)
    if (identical(res$defaults$model, 'BM1')) {
        regs <- regs[-1]
    }
    if (length(regs) != length(theta)) {
        return(NULL)
    }
    names(theta) <- regs
    paras$theta <- theta
    return(paras)
}

#' Regimes painted on a tree that a fitted theta does not cover
#' A regime with no retained tips has no estimable theta. Transferring such a fit onto the
#' full tree would leave part of the weights matrix undefined, so callers should treat a
#' non-empty result as a failed replicate rather than proceeding.
#' @param defaults Defaults list from preprocessTree for the target (full) tree.
#' @param theta Named theta vector from the masked fit.
#' @return Character vector of uncovered regimes.
check_regime_coverage <- function(defaults, theta) {
    tree_regs <- setdiff(as.character(defaults$parsed_alt_tree$unique_regimes), 'root')
    fit_regs <- setdiff(names(theta), 'root')
    return(setdiff(tree_regs, fit_regs))
}

#' Posterior over held-out tips given only the retained tips
#' @param paras List with tau, alpha, sigma, theta -- typically the `paras` element of a fit
#'   produced on a masked tree.
#' @param defaults Defaults list from preprocessTree for the FULL tree. The prediction is made
#'   on the full tree's geometry at the masked tree's parameters.
#' @param held_idx Indices of the held-out tips, relative to the tip order of the tree inside
#'   `defaults` (i.e. the cladewise-reordered tree returned by preprocessTree).
#' @param X_full Named numeric vector of observed values for every tip of the full tree, in
#'   that same tip order.
#' @param kept_idx Indices to condition on. Defaults to the complement of `held_idx`. Passing
#'   the full index set for both reduces this function to e_step.
#' @param verbose Whether to print numerical failures.
#' @return List with `mean`, `var`, `prior_mean` (the W theta term alone) and the indices
#'   used, or NULL if the covariance could not be inverted.
predict_masked_latent <- function(paras, defaults, held_idx, X_full, kept_idx = NULL, verbose = FALSE) {
    n <- defaults$parsed_alt_tree$n_leaves
    if (length(X_full) != n) {
        stop(sprintf('X has %d entries but the tree has %d tips.', length(X_full), n))
    }

    # The returned rows follow the caller's order, so do NOT sort -- callers compare against
    # X_full[held_idx] and a reordering here would silently misalign them.
    held_idx <- as.integer(held_idx)
    if (length(held_idx) == 0 || any(held_idx < 1) || any(held_idx > n)) {
        stop('held_idx must be a non-empty set of tip indices within the tree.')
    }
    if (anyDuplicated(held_idx)) {
        stop('held_idx contains duplicate tips.')
    }
    if (is.null(kept_idx)) {
        kept_idx <- setdiff(seq_len(n), held_idx)
    } else {
        kept_idx <- as.integer(kept_idx)
    }
    if (length(kept_idx) == 0) {
        stop('No tips left to condition on.')
    }
    if (paras$tau < 0) {
        stop('Error, tau is < 0.')
    }

    defaults <- align_regimes(defaults, paras$theta)

    V <- compute_VCV(defaults$parsed_alt_tree, paras$alpha, paras$sigma, add.root = defaults$add.root)
    W_ <- compute_W_matrix(defaults$parsed_alt_tree, paras$alpha, add.root = defaults$add.root)

    # BM1 with a fixed root carries theta twice, as in e_step.
    if (defaults$add.root == TRUE & defaults$model == 'BM1') {
        theta <- c(paras$theta, paras$theta)
    } else {
        theta <- paras$theta
    }
    if (dim(W_)[2] != length(theta)) {
        stop('Warning: Dimensions of weights matrix does not match inputted theta.\n')
    }

    mu <- as.vector(W_ %*% theta)

    n_kept <- length(kept_idx)
    ridge <- 1e-8 * diag(n_kept)
    V_oo <- V[kept_idx, kept_idx, drop = FALSE] + diag(paras$tau^2, n_kept) + ridge
    V_ho <- V[held_idx, kept_idx, drop = FALSE]

    # Solve once for both the mean correction and the variance reduction.
    rhs <- cbind(X_full[kept_idx] - mu[kept_idx], t(V_ho))

    L <- tryCatch(chol(V_oo), error = function(err) {
        if (verbose) print('Warning: tree likely bad. Cannot calculate inverse with cholesky decomp.')
        return(NULL)
    })
    if (is.null(L)) {
        A <- tryCatch(pseudoinverse(V_oo) %*% rhs, error = function(err) {
            if (verbose) print(err$message)
            return(NULL)
        })
    } else {
        A <- backsolve(L, backsolve(L, rhs, transpose = TRUE))
    }
    if (is.null(A)) {
        return(NULL)
    }

    post_mean <- as.vector(mu[held_idx] + V_ho %*% A[, 1])
    post_var <- diag(V[held_idx, held_idx, drop = FALSE] - V_ho %*% A[, -1, drop = FALSE])
    post_var <- pmax(as.vector(post_var), 0)

    names(post_mean) <- names(X_full)[held_idx]
    names(post_var) <- names(X_full)[held_idx]

    return(list(mean = post_mean,
                var = post_var,
                prior_mean = mu[held_idx],
                held_idx = held_idx,
                kept_idx = kept_idx))
}


######## SCORING ########

#' Score a set of predictions against the simulated truth
#' @param truth Numeric vector of true values at the held-out tips.
#' @param pred Numeric vector of predicted values.
#' @param pred_sd Optional posterior SD of the latent process, used for `coverage`.
#' @param pred_sd_obs Optional predictive SD on the OBSERVED scale, i.e. sqrt(Var[Z_h|X_o] + tau^2),
#'   used for `coverage_obs`. On real data the target is X_h = Z_h + eps, so `coverage` is
#'   structurally low for reasons unrelated to the model and `coverage_obs` is the honest number.
#' @param level Coverage level.
#' @return One-row data.frame of metrics.
score_recovery <- function(truth, pred, pred_sd = NULL, pred_sd_obs = NULL, level = 0.95) {
    ok <- is.finite(truth) & is.finite(pred)
    n <- sum(ok)

    out <- data.frame(n = n, rmse = NA_real_, mae = NA_real_, bias = NA_real_,
                      pearson = NA_real_, spearman = NA_real_,
                      mean_error = NA_real_, coverage = NA_real_, coverage_obs = NA_real_)
    if (n < 1) {
        return(out)
    }

    tr <- truth[ok]
    pr <- pred[ok]
    err <- pr - tr

    out$rmse <- sqrt(mean(err^2))
    out$mae <- mean(abs(err))
    out$bias <- mean(err)
    # For a masked clade this is the informative one: within-clade structure is not
    # recoverable, but the clade's overall level is.
    out$mean_error <- mean(pr) - mean(tr)

    # cor() warns and returns NA on a constant vector, which is the expected case for the
    # clade arm on an ultrametric tree.
    if (n >= 3 && sd(tr) > 0 && sd(pr) > 0) {
        out$pearson <- suppressWarnings(cor(tr, pr))
        out$spearman <- suppressWarnings(cor(tr, pr, method = 'spearman'))
    }

    z <- qnorm(1 - (1 - level) / 2)
    if (!is.null(pred_sd)) {
        psd <- pred_sd[ok]
        valid <- is.finite(psd)
        if (any(valid)) {
            out$coverage <- mean(abs(err[valid]) <= z * psd[valid])
        }
    }
    if (!is.null(pred_sd_obs)) {
        psd <- pred_sd_obs[ok]
        valid <- is.finite(psd)
        if (any(valid)) {
            out$coverage_obs <- mean(abs(err[valid]) <= z * psd[valid])
        }
    }

    return(out)
}

#' AIC for a single fit, using the real tip count
#' annotate_history() hardcodes ntips = 256, which is wrong for any masked tree, so model
#' selection here is computed from scratch.
#' @param res One element of a runSCOUT result list.
#' @param ntips Number of tips the fit was actually run on.
#' @return List with ll, k, aic and aicc.
fit_aic <- function(res, ntips) {
    if (is.null(res) || is.null(res$paras)) {
        return(list(ll = NA_real_, k = NA_integer_, aic = NA_real_, aicc = NA_real_))
    }
    ll <- suppressWarnings(as.numeric(res$final$ll_total))
    k <- length(unlist(res$paras))
    if (identical(res$defaults$model, 'BM1')) {
        k <- k - 1
    }
    if (length(ll) != 1 || !is.finite(ll)) {
        return(list(ll = NA_real_, k = k, aic = NA_real_, aicc = NA_real_))
    }
    aic <- -2 * ll + 2 * k
    aicc <- if (ntips - k - 1 > 0) -2 * ll + (2 * k * (ntips / (ntips - k - 1))) else NA_real_
    return(list(ll = ll, k = k, aic = aic, aicc = aicc))
}


######## ORCHESTRATION ########

#' Prepare a tree for the dropout workflow
#' Enforces the two contracts simulate_test_data requires of an externally generated tree:
#' labelled internal nodes and tip states.
#' @param tree Path to a newick file or an object of class phylo.
#' @param states Optional named vector of tip regimes.
#' @param logfile Optional log file.
#' @param verbose Whether to log.
#' @param infer_nodes Whether to reconstruct ancestral states when the tree has no node labels.
#'   The simulator needs them; formatSCOUT does not, because it re-derives node labels per regime
#'   and resets `node.label` to NULL on every iteration.
#' @return A phylo object with `states` and `node.label` set.
#' @export
prep_dropout_tree <- function(tree, states = NULL, logfile = NULL, verbose = TRUE,
                              infer_nodes = TRUE) {
    if (inherits(tree, 'phylo')) {
        phy <- tree
    } else {
        phy <- ape::read.tree(tree)
    }
    if (is.null(phy) || !inherits(phy, 'phylo')) {
        stop('Could not read a tree. Provide a newick file path or a single phylo object.')
    }

    # Resolve states against the tip labels before any topology edits.
    if (is.null(states)) {
        if (!'states' %in% names(phy)) {
            stop('No tip states found. Supply `states` or attach them to the tree as `tree$states`.')
        }
        states <- phy$states
    }
    if (is.null(names(states))) {
        if (length(states) != length(phy$tip.label)) {
            stop('`states` is unnamed, so it must have one entry per tip in tip.label order.')
        }
        names(states) <- phy$tip.label
    }
    missing_states <- setdiff(phy$tip.label, names(states))
    if (length(missing_states) > 0) {
        stop(sprintf('Missing states for %d tips, e.g. %s', length(missing_states),
                     paste(head(missing_states, 5), collapse = ', ')))
    }
    states <- setNames(as.character(states), names(states))

    # ape's ancestral state reconstruction requires a fully dichotomous tree.
    if (phy$Nnode != length(phy$tip.label) - 1) {
        log_message('Tree is not fully dichotomous. Resolving polytomies with multi2di.', logfile, verbose = verbose)
        phy <- multi2di(phy)
    }

    if (is.null(phy$edge.length)) {
        log_message('No edge lengths, replacing with 1s.', logfile, verbose = verbose)
        phy$edge.length <- rep(1, Nedge(phy))
    } else if (sum(phy$edge.length == 0) > 0) {
        log_message('Found 0 length edges. Replacing 0 edge lengths with 1e-7.', logfile, verbose = verbose)
        phy$edge.length[phy$edge.length == 0] <- 1e-7
    }

    phy$states <- as.character(states[phy$tip.label])

    if (infer_nodes && is.null(phy$node.label)) {
        log_message('Internal nodes are unlabelled. Inferring ancestral states with infer_anc.', logfile, verbose = verbose)
        phy <- infer_anc(phy)
    }

    return(phy)
}

#' Index a runSCOUT result list by gene and regime
#' @param fit A runSCOUT result list.
#' @return Named list keyed by "<gene>||<regime>".
index_fits <- function(fit) {
    out <- list()
    for (res in fit) {
        out[[paste(res$settings[[1]], res$settings[[2]], sep = '||')]] <- res
    }
    return(out)
}

#' Held-out validation of SCOUT-EM by masking part of the tree
#'
#' Repeatedly masks part of a tree, refits SCOUT-EM on the reduced tree, and predicts the masked
#' tips' latent expression from the retained tips alone. Predictions are scored against the truth
#' and against baselines.
#'
#' Runs in one of two modes. With `data = NULL` expression is **simulated** on the tree and scored
#' against the noiseless latent EVFs. With `data` supplied the **observed** matrix is fit and
#' scored against the observed values at the masked tips -- which are X_h = Z_h + eps, so RMSE
#' carries an irreducible floor near tau and the result must be read from the comparison between
#' predictors rather than from the absolute number.
#'
#' The four predictors are: `masked` (masked-fit parameters conditioned on the retained tips --
#' the thing under test), `prior` (the W_h theta term alone, so `masked - prior` isolates what the
#' phylogenetic covariance V_ho contributes), `oracle` (full-tree parameters, still conditioned on
#' retained tips only -- it never sees X_h, and bounds how much of the gap is parameter-estimation
#' error rather than an information limit), and `global` (the retained-tip mean; no model).
#'
#' Two masking arms are run per replicate. In the `clade` arm every tip below a sampled node is
#' masked; because a monophyletic clade joins the rest of the tree through a single branch, V_ho
#' is rank one. On an ultrametric tree that makes within-clade predictions constant, so `pearson`
#' is undefined and `mean_error` -- the clade's overall level -- is the informative metric. On a
#' NON-ultrametric tree the predictions instead form a one-parameter family indexed by tip depth,
#' so `pearson` will be near +/-1 by construction and reflects depth rather than clade-internal
#' structure; do not read it as recovery. The `random` arm masks the same number of tips scattered
#' across the tree, where V_ho has full rank and per-tip recovery is testable.
#'
#' The random arm can be repeated at several fixed mask sizes in one run via `random_frac`, giving
#' a dose-response over how much of the tree is removed. That needs its own argument because clade
#' sizes cannot be dialled to arbitrary values -- a tree only contains the clade sizes its topology
#' happens to offer, and on a 256-tip birth tree the eligible sizes typically leave gaps wide
#' enough that whole percentage bands are unreachable.
#'
#' @param tree Path to a newick file, or a phylo object.
#' @param results_dir Directory for output files.
#' @param states Named vector of tip regimes. Optional if the tree carries `tree$states`, or in
#'   data mode where it is derived from `mask_states`.
#' @param data Real expression data: a data.frame or a path to a CSV with `row.names = 1`. When
#'   supplied the simulation is skipped and this matrix is both fit and scored against.
#' @param species_key Column of `data` holding the tip labels.
#' @param blacklist Non-gene columns of `data` to exclude, passed to formatSCOUT. Omitting the
#'   character columns here makes them "genes" and poisons `gene_cols` with NAs.
#' @param normalize Whether formatSCOUT should log1p the matrix. Leave FALSE for data that is
#'   already normalised.
#' @param genes Gene subset, passed to formatSCOUT as `quant_traits`. Cost is linear in this, so
#'   on a large matrix it is effectively required -- see the timing note in the vignette.
#' @param mask_states Column names of `data` combined by `dropout_states` into the composite tip
#'   state used for mask eligibility. Defaults to the OUM-style entries of `regimes`.
#' @param mask_nodes Explicit internal node numbers to mask, one per replicate, instead of
#'   sampling. Overrides `n_replicates`, `clade_min` and `clade_max`. Use `eligible_clades` to
#'   choose them. This also decouples mask choice from the RNG stream that runSCOUT consumes in
#'   proportion to the gene count, so masks stay fixed when the gene set changes.
#' @param mask_draw When are clades drawn: `'upfront'` draws all of them before any fitting, so
#'   no clade repeats and the draw does not depend on how much RNG runSCOUT consumes (which
#'   scales with the gene count). `'inline'` draws one per replicate as the loop runs, which is
#'   the historical behaviour and is kept so older seeded runs stay reproducible. Defaults to
#'   `'inline'` in sim mode and `'upfront'` in data mode.
#' @param full_fit Cached reference pass to reuse: a list with `fit`, or a path to one saved by
#'   `save_full_fit`. Skips refitting the full tree.
#' @param save_full_fit Whether to write the full-tree reference fit for later reuse.
#' @param checkpoint Whether to write a rolling `<testid>_checkpoint.rds` after every replicate.
#' @param sim_res Result of simulate_test_data. When NULL, data are simulated internally.
#'   Supplying one reuses a single simulation across every replicate, which isolates the effect
#'   of clade choice from simulation noise. Ignored in data mode.
#' @param sim_key Which parameter combination of `sim_res` to use. Defaults to the first.
#' @param regimes Regimes to fit, passed to formatSCOUT.
#' @param n_replicates Number of masks to draw.
#' @param clade_min Minimum masked clade size.
#' @param clade_max Maximum masked clade size.
#' @param arms Which masking arms to run: any of 'clade' and 'random'.
#' @param random_frac Mask sizes for the random arm, as fractions of the tip count. `NULL` (the
#'   default) keeps the historical behaviour of one random arm size-matched to the replicate's
#'   clade. A vector runs one random arm per entry, named `random_05pct`, `random_10pct` and so
#'   on; an `NA` entry means the size-matched arm, so `c(NA, 0.05, 0.2)` runs the matched control
#'   alongside the sweep. Each entry must be NA or in (0, 1), and the largest must still leave
#'   more than 10 tips. Ignored when `arms` excludes 'random'. Note a clade is still drawn every
#'   replicate even for a random-only run, since the matched arm and the replicate's bookkeeping
#'   depend on it.
#' @param keep_predictions Which per-tip predictions to retain: 'ref_model' keeps only each gene's
#'   reference regime, 'all' keeps every regime, 'none' keeps metrics only. 'true_model' is
#'   accepted as an alias for 'ref_model'.
#' @param ngenes,nevfs,a,s,t0,theta_step Simulation settings, passed to simulate_test_data.
#'   Only used when `sim_res` is NULL and `data` is NULL.
#' @param lambda1,lambda2,fixed_root,tau_prior_mean,tau_prior_sd Fitting settings, passed to
#'   runSCOUT.
#' @param randseed Seed set once at the start of the run.
#' @param testid Run ID. Generated when NULL.
#' @param cores Cores for runSCOUT.
#' @param logfile Optional log file.
#' @param verbose Whether to log progress.
#' @return List with `clades`, `predictions`, `metrics`, `params`, `model_select` and
#'   `settings`. The same tables are written to `results_dir`. In `metrics` and `predictions`,
#'   `is_selected` marks the regime the masked fit itself chose by AIC -- filter on it for the
#'   headline result, since that selection never saw the held-out tips.
#' @import ape
#' @import stringr
#' @import dplyr
#' @export
runSCOUT.dropout <- function(tree,
    results_dir,
    states = NULL,
    data = NULL,
    species_key = 'species',
    blacklist = NULL,
    normalize = FALSE,
    genes = NULL,
    mask_states = NULL,
    mask_nodes = NULL,
    mask_draw = NULL,
    full_fit = NULL,
    save_full_fit = TRUE,
    checkpoint = TRUE,
    sim_res = NULL,
    sim_key = NULL,
    regimes = c('BM1', 'OU1', 'OUM'),
    n_replicates = 20,
    clade_min = 5,
    clade_max = 30,
    arms = c('clade', 'random'),
    random_frac = NULL,
    keep_predictions = 'true_model',
    ngenes = 20,
    nevfs = 20,
    a = 3,
    s = 1,
    t0 = 3,
    theta_step = 2,
    lambda1 = 0.2,
    lambda2 = 0.2,
    fixed_root = FALSE,
    tau_prior_mean = 0.2,
    tau_prior_sd = 0.1,
    randseed = NULL,
    testid = NULL,
    cores = 1,
    logfile = NULL,
    verbose = TRUE) {

    if (!is.null(randseed)) {
        set.seed(randseed)
    }
    mode <- if (is.null(data)) 'sim' else 'data'
    if (mode == 'data' && !is.null(sim_res)) {
        stop('Supply either `data` or `sim_res`, not both.')
    }
    if (is.null(mask_draw)) {
        mask_draw <- if (mode == 'data') 'upfront' else 'inline'
    }
    if (!mask_draw %in% c('upfront', 'inline')) {
        stop("mask_draw must be 'upfront' or 'inline'.")
    }
    if (identical(keep_predictions, 'true_model')) {
        keep_predictions <- 'ref_model'
    }
    if (!keep_predictions %in% c('ref_model', 'all', 'none')) {
        stop("keep_predictions must be one of 'ref_model', 'all' or 'none'.")
    }
    arms <- intersect(arms, c('clade', 'random'))
    if (length(arms) == 0) {
        stop("arms must contain at least one of 'clade' and 'random'.")
    }
    if (!is.null(random_frac)) {
        random_frac <- as.numeric(random_frac)
        if (length(random_frac) == 0) {
            stop('random_frac must have at least one entry, or be NULL.')
        }
        bad <- !is.na(random_frac) & (random_frac <= 0 | random_frac >= 1)
        if (any(bad)) {
            stop('random_frac entries must be NA (size-matched to the clade) or in (0, 1).')
        }
        if (anyDuplicated(random_frac)) {
            stop('random_frac has duplicate entries; each one becomes an arm and arm names must be unique.')
        }
    }

    create_directory_if_not_exists(results_dir)
    tid <- if (is.null(testid)) {
        paste0('SCOUTDROP_', paste(sample(c(0:9, letters, LETTERS), 8, replace = TRUE), collapse = ""))
    } else {
        testid
    }
    log_message(sprintf('Dropout validation run ID = %s', tid), logfile, verbose = TRUE)
    log_message(sprintf('Data will be saved to --> %s', results_dir), logfile, verbose = verbose)
    log_message('=========================================================\n', logfile, verbose = verbose)

    ############### Data (data mode) ###############
    # Read before the tree, because the composite mask-eligibility states come from the metadata.
    if (mode == 'data') {
        # as.data.frame is not cosmetic: formatSCOUT tests `class(x) == 'data.frame'`, which is a
        # length-3 condition for a tibble and errors under R >= 4.2.
        meta <- as.data.frame(if (is.data.frame(data)) data else read.csv(data, row.names = 1))
        if (!species_key %in% colnames(meta)) {
            stop(sprintf('species_key `%s` is not a column of `data`.', species_key))
        }
        if (anyDuplicated(meta[[species_key]])) {
            stop('`data` has duplicate species. formatSCOUT uses them as row names.')
        }
        if (is.null(states)) {
            mcols <- if (is.null(mask_states)) setdiff(regimes, c('BM1', 'OU1')) else mask_states
            if (length(mcols) == 0) {
                stop('No regime columns to build mask states from. Set `mask_states` or `states`.')
            }
            states <- dropout_states(meta, mcols, species_key)
            log_message(sprintf('Mask eligibility states from: %s', paste(mcols, collapse = ', ')),
                        logfile, verbose = verbose)
        }
    }

    ############### Tree ###############
    # Ancestral states are only needed by the simulator; formatSCOUT re-derives node labels per
    # regime and resets node.label to NULL each iteration, so in data mode inferring them here is
    # wasted work and one more failure surface.
    phy <- prep_dropout_tree(tree, states, logfile, verbose, infer_nodes = (mode == 'sim'))
    n_tips <- length(phy$tip.label)
    tip_states <- setNames(as.character(phy$states), phy$tip.label)
    log_message(sprintf('Tree has %d tips across %d regimes.', n_tips, length(unique(tip_states))),
                logfile, verbose = verbose)
    small <- names(which(table(tip_states) < 3))
    if (length(small) > 0) {
        log_message(sprintf('Regimes with fewer than 3 tips (fragile under masking): %s',
                            paste(small, collapse = ', ')), logfile, verbose = TRUE)
    }

    if (is.null(mask_nodes) && clade_max >= n_tips - 10) {
        stop(sprintf('clade_max (%d) leaves too few tips to fit on out of %d.', clade_max, n_tips))
    }
    if (!is.null(random_frac)) {
        k_rand <- pmax(1L, round(random_frac[!is.na(random_frac)] * n_tips))
        if (any(k_rand >= n_tips - 10)) {
            stop(sprintf('random_frac up to %.2f masks %d of %d tips, leaving too few to fit on.',
                         max(random_frac, na.rm = TRUE), max(k_rand), n_tips))
        }
    }

    ############### Simulate (sim mode) ###############
    if (mode == 'sim') {
        if (is.null(sim_res)) {
            sim_dir <- file.path(results_dir, 'simulated_data')
            # simulate_test_data does not create outdir when sim_lin = FALSE, but still writes to it.
            create_directory_if_not_exists(sim_dir)
            log_message('Simulating expression on the complete tree.', logfile, verbose = verbose)
            sim_res <- simulate_test_data(ngenes = ngenes, ncells = n_tips, tree = phy,
                                          outdir = sim_dir, out_prefix = tid,
                                          a = a, s = s, t0 = t0, theta_step = theta_step,
                                          nevfs = nevfs, sim_lin = FALSE)
        }

        if (is.null(sim_key)) {
            sim_key <- names(sim_res$simulated_OU$evfs)[1]
        }
        if (!sim_key %in% names(sim_res$simulated_OU$evfs)) {
            stop(sprintf('sim_key `%s` not found. Available: %s', sim_key,
                         paste(names(sim_res$simulated_OU$evfs), collapse = ', ')))
        }
        # The EVF matrix is the raw OUwie.sim output, i.e. Z itself with no observation noise, so
        # errors against it are on the same scale as the latent process. The counts matrix is a
        # beta-Poisson draw from a mixture of EVFs and has no 1:1 gene-to-latent correspondence.
        meta <- sim_res$simulated_OU$evfs[[sim_key]]
        rownames(meta) <- meta$species
        if (!setequal(meta$species, phy$tip.label)) {
            stop('Simulated data and tree tip labels do not match. Was sim_res generated on this tree?')
        }

        sim_gene_cols <- setdiff(colnames(meta), c('OUM', 'species'))
        min_trait <- min(meta[, sim_gene_cols], na.rm = TRUE)
        if (min_trait <= 0) {
            warning(sprintf('Simulated traits reach %.4f. m_step_tipfog bounds theta below at 1e-10, so non-positive traits cannot be fit -- raise t0.', min_trait))
        }
    } else {
        if (!setequal(meta[[species_key]], phy$tip.label)) {
            stop(sprintf('Data and tree tip labels do not match: %d tips, %d species, %d shared.',
                         n_tips, length(meta[[species_key]]),
                         length(intersect(meta[[species_key]], phy$tip.label))))
        }
    }

    ############### Full-tree reference ###############
    log_message('Preparing the complete tree.', logfile, verbose = verbose)
    idata_full <- formatSCOUT(tree_path = phy, metadata_path = meta, outpath = results_dir,
                              species_key = species_key, quant_traits = genes, regimes = regimes,
                              anc_infer = 'ape', normalize = normalize, blacklist = blacklist,
                              logfile = logfile)
    gene_cols <- idata_full$gene_cols
    # A character column that escaped `blacklist` becomes NA here rather than failing loudly:
    # formatSCOUT coerces the gene block through as.matrix, and the zero-sum filter then injects
    # NA into gene_cols.
    if (length(gene_cols) == 0 || anyNA(gene_cols)) {
        stop(sprintf('gene_cols is unusable (%d entries, %d NA). Check `blacklist` and `genes`.',
                     length(gene_cols), sum(is.na(gene_cols))))
    }
    log_message(sprintf('Fitting %d genes x %d regimes.', length(gene_cols), length(regimes)),
                logfile, verbose = verbose)

    if (mode == 'data') {
        gm <- as.matrix(idata_full$meta_data[, gene_cols, drop = FALSE])
        log_message(sprintf('Observed data: range [%.3f, %.3f], %.1f%% zeros, %d genes whose tree-wide mean is below 1e-7 (theta sits on its lower bound there).',
                            min(gm, na.rm = TRUE), max(gm, na.rm = TRUE), 100 * mean(gm == 0, na.rm = TRUE),
                            sum(colMeans(gm, na.rm = TRUE) < 1e-7)),
                    logfile, verbose = verbose)
    }

    # scaleHeight is deliberately FALSE: it divides edge lengths by Tmax, so parameters fit on
    # a masked tree of a different height would not transfer onto the full tree's covariance.
    tree_info_full <- list()
    for (reg in regimes) {
        tree_info_full[[reg]] <- preprocessTree(inputs = idata_full, reg = reg, scaleHeight = FALSE,
                                                root.fixed = fixed_root, root.age = NULL, skipTau = FALSE)
    }
    # Every quantity in parsed_alt_tree is indexed by the cladewise tip order preprocessTree
    # produces, so held-out indices are only meaningful if all regimes share that order.
    full_tip_order <- tree_info_full[[regimes[1]]]$tree_const$tree$tip.label
    for (reg in regimes) {
        if (!identical(tree_info_full[[reg]]$tree_const$tree$tip.label, full_tip_order)) {
            stop(sprintf('Regime %s reorders the tips relative to %s. Cannot index held-out tips consistently.', reg, regimes[1]))
        }
    }

    if (is.null(full_fit)) {
        log_message('Fitting SCOUT-EM on the complete tree (reference).', logfile, verbose = verbose)
        fit_full <- runSCOUT(idata_full, lambda1 = lambda1, lambda2 = lambda2, fixed.root = fixed_root,
                             tau_prior_mean = tau_prior_mean, tau_prior_sd = tau_prior_sd,
                             scaleHeight = FALSE, runid = paste0(tid, '_full'),
                             cores = cores, logfile = logfile, verbose = FALSE)
        if (save_full_fit) {
            saveRDS(list(idata = idata_full, fit = fit_full),
                    sprintf('%s/%s_full_fit.rds', results_dir, tid))
        }
    } else {
        log_message('Reusing a cached full-tree reference fit.', logfile, verbose = verbose)
        cached <- if (is.character(full_fit)) readRDS(full_fit) else full_fit
        fit_full <- if (is.null(cached$fit)) cached else cached$fit
    }
    fits_full <- index_fits(fit_full)
    want <- as.vector(outer(gene_cols, regimes, paste, sep = '||'))
    if (!all(want %in% names(fits_full))) {
        stop(sprintf('The reference fit is missing %d of %d gene x regime combinations.',
                     sum(!want %in% names(fits_full)), length(want)))
    }

    # X for every gene, in the tip order the full tree's model objects are indexed by.
    X_all <- list()
    for (gene in gene_cols) {
        xv <- setNames(as.numeric(idata_full$meta_data[, gene]), idata_full$meta_data[, 'species'])
        X_all[[gene]] <- xv[full_tip_order]
    }

    best_model_full <- sapply(gene_cols, function(gene) {
        aics <- sapply(regimes, function(reg) fit_aic(fits_full[[paste(gene, reg, sep = '||')]], n_tips)$aic)
        if (all(is.na(aics))) NA_character_ else regimes[which.min(aics)]
    })

    # Each gene's reference regime -- the one `keep_predictions = 'ref_model'` retains. On
    # simulated data it is encoded in the gene name; on real data there is no such thing, so the
    # full-tree AIC winner stands in. (This is only a reporting convenience: the headline result
    # keys off `is_selected`, the regime the MASKED fit chose, which never saw the held-out tips.)
    ref_model <- if (mode == 'sim') {
        setNames(str_extract(gene_cols, paste(regimes, collapse = '|')), gene_cols)
    } else {
        best_model_full
    }
    if (anyNA(ref_model)) {
        log_message(sprintf('%d genes have no model encoded in their name; using best_model_full for those.',
                            sum(is.na(ref_model))), logfile, verbose = verbose)
        ref_model[is.na(ref_model)] <- best_model_full[is.na(ref_model)]
    }
    ref_source <- if (mode == 'sim') 'simulated_generator' else 'best_model_full'

    ############### Masks ###############
    if (!is.null(mask_nodes)) {
        mask_nodes <- as.integer(mask_nodes)
        n_replicates <- length(mask_nodes)
        clade_list <- lapply(mask_nodes, function(nd) list(node = nd, tips = clade_tips(phy, nd)))
        log_message(sprintf('Using %d caller-supplied clades, sizes %s.', n_replicates,
                            paste(range(lengths(lapply(clade_list, `[[`, 'tips'))), collapse = '-')),
                    logfile, verbose = verbose)
    } else if (mask_draw == 'upfront') {
        # Drawn without replacement before any fitting, so no clade repeats and the draw does not
        # depend on the RNG runSCOUT consumes in proportion to the gene count.
        cand <- eligible_clades(phy, clade_min, clade_max, states = tip_states,
                                exclude_root_children = TRUE)
        take <- if (nrow(cand) == 0) integer(0) else sample(nrow(cand), min(n_replicates, nrow(cand)))
        clade_list <- lapply(take, function(i) list(node = cand$node[i], tips = attr(cand, 'tips')[[i]]))
        if (length(clade_list) < n_replicates) {
            log_message(sprintf('Only %d of %d requested clades are eligible; running %d replicates.',
                                length(clade_list), n_replicates, length(clade_list)),
                        logfile, verbose = TRUE)
            n_replicates <- length(clade_list)
        }
    } else {
        # Historical behaviour: one independent draw per replicate, interleaved with fitting.
        clade_list <- NULL
    }

    ############### Replicates ###############
    clade_rows <- list()
    pred_rows <- list()
    metric_rows <- list()
    param_rows <- list()
    select_rows <- list()

    for (r in seq_len(n_replicates)) {
        clade <- if (is.null(clade_list)) {
            sample_clade(phy, clade_min, clade_max, states = tip_states, exclude_root_children = TRUE)
        } else {
            clade_list[[r]]
        }
        if (is.null(clade)) {
            clade_rows[[length(clade_rows) + 1]] <- data.frame(
                replicate = r, arm = NA_character_, node = NA_integer_, n_masked = NA_integer_,
                masked_tips = NA_character_, status = 'no_eligible_clade')
            log_message(sprintf('Replicate %d: no clade satisfies the size and regime constraints.', r),
                        logfile, verbose = verbose)
            next
        }

        masks <- list()
        if ('clade' %in% arms) {
            masks[['clade']] <- list(node = clade$node, tips = clade$tips)
        }
        if ('random' %in% arms) {
            # NA means the historical behaviour: k matched to this replicate's clade, which is the
            # like-for-like control for it. A fraction instead fixes k across replicates, which is
            # what a dose-response over mask size needs -- clade sizes cannot be dialled freely,
            # since the topology only offers the clade sizes it happens to contain.
            fracs <- if (is.null(random_frac)) NA_real_ else random_frac
            for (p in fracs) {
                nm <- if (is.na(p)) 'random' else sprintf('random_%02.0fpct', 100 * p)
                k_r <- if (is.na(p)) length(clade$tips) else max(1L, round(p * n_tips))
                rtips <- sample_random_tips(phy, k_r, states = tip_states)
                if (is.null(rtips)) {
                    log_message(sprintf('Replicate %d | arm %s: no draw of %d tips keeps every regime; skipping.',
                                        r, nm, k_r), logfile, verbose = TRUE)
                    next
                }
                masks[[nm]] <- list(node = NA_integer_, tips = rtips)
            }
        }

        for (arm in names(masks)) {
            held_tips <- masks[[arm]]$tips
            k <- length(held_tips)
            held_idx <- match(held_tips, full_tip_order)
            log_message(sprintf('Replicate %d | arm %s | masking %d of %d tips.', r, arm, k, n_tips),
                        logfile, verbose = verbose)

            phy_m <- drop.tip(phy, held_tips)
            phy_m$states <- NULL # formatSCOUT re-derives this from the metadata

            # The whole method rests on the masked tree being on the same time scale as the full
            # one, so parameters fit here transfer onto the full tree's covariance. drop.tip
            # collapses degree-two nodes, which is fine, but if it ever removed the root the
            # depths would shift silently.
            d_f <- setNames(node.depth.edgelength(phy)[seq_len(n_tips)], phy$tip.label)
            d_m <- setNames(node.depth.edgelength(phy_m)[seq_along(phy_m$tip.label)], phy_m$tip.label)
            if (!isTRUE(all.equal(d_m, d_f[phy_m$tip.label]))) {
                stop(sprintf('Replicate %d arm %s: dropping the mask changed root-to-tip distances, so masked-tree parameters no longer transfer.', r, arm))
            }

            idata_m <- formatSCOUT(tree_path = phy_m, metadata_path = meta, outpath = results_dir,
                                   species_key = species_key, quant_traits = gene_cols,
                                   regimes = regimes, anc_infer = 'ape', normalize = normalize,
                                   blacklist = blacklist, logfile = logfile)
            # A gene can go all-zero on the retained tips; formatSCOUT drops those.
            dropped <- setdiff(gene_cols, idata_m$gene_cols)
            if (length(dropped) > 0) {
                log_message(sprintf('  %d genes are all-zero on the retained tips and were dropped.',
                                    length(dropped)), logfile, verbose = verbose)
            }
            fit_m <- runSCOUT(idata_m, lambda1 = lambda1, lambda2 = lambda2, fixed.root = fixed_root,
                              tau_prior_mean = tau_prior_mean, tau_prior_sd = tau_prior_sd,
                              scaleHeight = FALSE, runid = sprintf('%s_r%03d_%s', tid, r, arm),
                              cores = cores, logfile = logfile, verbose = FALSE)
            fits_m <- index_fits(fit_m)
            n_kept <- length(phy_m$tip.label)

            tally <- c(ok = 0, fit_missing = 0, regime_lost = 0, prediction_failed = 0,
                       gene_dropped = 0)
            for (gene in gene_cols) {
                ref <- unname(ref_model[gene])
                X <- X_all[[gene]]
                if (gene %in% dropped) {
                    tally['gene_dropped'] <- tally['gene_dropped'] + 1
                    next
                }

                aics_m <- rep(NA_real_, length(regimes))
                names(aics_m) <- regimes
                scored <- list()

                for (reg in regimes) {
                    key <- paste(gene, reg, sep = '||')
                    res_m <- fits_m[[key]]
                    res_f <- fits_full[[key]]
                    defaults_full <- tree_info_full[[reg]]$defaults

                    aic_m <- fit_aic(res_m, n_kept)
                    aics_m[reg] <- aic_m$aic
                    converge_m <- if (is.null(res_m$final$converge)) NA_character_ else as.character(res_m$final$converge)

                    param_rows[[length(param_rows) + 1]] <- data.frame(
                        replicate = r, arm = arm, gene = gene, regime = reg, true_model = ref,
                        alpha = res_m$paras$alpha %||% NA_real_,
                        sigma = res_m$paras$sigma %||% NA_real_,
                        tau = res_m$paras$tau %||% NA_real_,
                        ll_total = aic_m$ll, param_count = aic_m$k, aic = aic_m$aic, aicc = aic_m$aicc,
                        converge = converge_m, n_kept = n_kept,
                        alpha_full = res_f$paras$alpha %||% NA_real_,
                        sigma_full = res_f$paras$sigma %||% NA_real_,
                        tau_full = res_f$paras$tau %||% NA_real_,
                        stringsAsFactors = FALSE)

                    paras_m <- fit_paras(res_m)
                    paras_f <- fit_paras(res_f)
                    if (is.null(paras_m) || is.null(paras_f)) {
                        tally['fit_missing'] <- tally['fit_missing'] + 1
                        aics_m[reg] <- NA_real_
                        next
                    }

                    # A regime with no retained tips has no estimable theta.
                    uncovered <- check_regime_coverage(defaults_full, paras_m$theta)
                    if (length(uncovered) > 0) {
                        tally['regime_lost'] <- tally['regime_lost'] + 1
                        aics_m[reg] <- NA_real_
                        next
                    }

                    pred_m <- predict_masked_latent(paras_m, defaults_full, held_idx, X)
                    pred_f <- predict_masked_latent(paras_f, defaults_full, held_idx, X)
                    if (is.null(pred_m) || is.null(pred_f)) {
                        tally['prediction_failed'] <- tally['prediction_failed'] + 1
                        aics_m[reg] <- NA_real_
                        next
                    }

                    z_true <- X[held_idx]
                    if (!identical(names(z_true), names(pred_m$mean)) ||
                        !identical(names(z_true), names(pred_f$mean))) {
                        stop('Predicted and true values are not aligned to the same tips.')
                    }
                    tally['ok'] <- tally['ok'] + 1
                    scored[[reg]] <- list(pred_m = pred_m, pred_f = pred_f, z_true = z_true,
                                          converge = converge_m, tau_m = paras_m$tau,
                                          tau_f = paras_f$tau,
                                          X_kept_mean = mean(X[pred_m$kept_idx]))
                }

                # Model selection on the MASKED fit only, so it never sees the held-out tips.
                # This is what `is_selected` marks and what the headline result filters on.
                best_m <- if (all(is.na(aics_m))) NA_character_ else regimes[which.min(aics_m)]
                select_rows[[length(select_rows) + 1]] <- data.frame(
                    replicate = r, arm = arm, gene = gene, true_model = ref,
                    best_full = unname(best_model_full[gene]), best_masked = best_m,
                    stringsAsFactors = FALSE)

                for (reg in names(scored)) {
                    s <- scored[[reg]]
                    sel <- identical(reg, best_m)
                    sd_lat <- sqrt(s$pred_m$var)
                    preds <- list(masked = s$pred_m$mean,
                                  prior = s$pred_m$prior_mean,
                                  oracle = s$pred_f$mean,
                                  global = rep(s$X_kept_mean, k))
                    # Each predictor is scored with its OWN posterior SD; the observed-scale
                    # variant adds tau^2, which is the right interval when the target carries
                    # tip fog (always, in data mode).
                    sds <- list(masked = sd_lat, prior = NULL,
                                oracle = sqrt(s$pred_f$var), global = NULL)
                    taus <- list(masked = s$tau_m, prior = NA_real_,
                                 oracle = s$tau_f, global = NA_real_)

                    for (pname in names(preds)) {
                        psd <- sds[[pname]]
                        sc <- score_recovery(s$z_true, preds[[pname]], pred_sd = psd,
                                             pred_sd_obs = if (is.null(psd)) NULL else sqrt(psd^2 + taus[[pname]]^2))
                        metric_rows[[length(metric_rows) + 1]] <- cbind(
                            data.frame(replicate = r, arm = arm, gene = gene, regime = reg,
                                       true_model = ref, is_selected = sel, predictor = pname,
                                       n_masked = k, n_kept = n_kept, converge = s$converge,
                                       stringsAsFactors = FALSE),
                            sc)
                    }

                    if (keep_predictions == 'all' || sel ||
                        (keep_predictions == 'ref_model' && identical(reg, ref))) {
                        pred_rows[[length(pred_rows) + 1]] <- data.frame(
                            replicate = r, arm = arm, gene = gene, regime = reg, true_model = ref,
                            is_selected = sel,
                            tip = names(s$z_true), z_true = as.numeric(s$z_true),
                            z_masked = as.numeric(s$pred_m$mean), z_prior = as.numeric(s$pred_m$prior_mean),
                            z_oracle = as.numeric(s$pred_f$mean), pred_sd = sd_lat,
                            pred_sd_obs = sqrt(sd_lat^2 + s$tau_m^2),
                            stringsAsFactors = FALSE)
                    }
                }
            }

            dominant <- names(tally)[which.max(replace(tally, 'ok', -1))]
            clade_rows[[length(clade_rows) + 1]] <- data.frame(
                replicate = r, arm = arm, node = masks[[arm]]$node, n_masked = k,
                masked_tips = paste(held_tips, collapse = ';'),
                status = if (max(tally[names(tally) != 'ok']) == 0) 'ok' else dominant,
                n_ok = unname(tally['ok']), n_fit_missing = unname(tally['fit_missing']),
                n_regime_lost = unname(tally['regime_lost']),
                n_prediction_failed = unname(tally['prediction_failed']),
                n_gene_dropped = unname(tally['gene_dropped']),
                stringsAsFactors = FALSE)
        }

        if (checkpoint) {
            saveRDS(list(clades = clade_rows, predictions = pred_rows, metrics = metric_rows,
                         params = param_rows, model_select = select_rows,
                         replicate_done = r, n_replicates = n_replicates),
                    sprintf('%s/%s_checkpoint.rds', results_dir, tid))
        }
    }

    ############### Collect and write ###############
    # bind_rows, not do.call(rbind, ): the metric list runs to tens of thousands of one-row
    # data.frames, where rbind is quadratic.
    bind <- function(rows) if (length(rows) == 0) NULL else as.data.frame(dplyr::bind_rows(rows))

    out <- list(clades = bind(clade_rows),
                predictions = bind(pred_rows),
                metrics = bind(metric_rows),
                params = bind(param_rows),
                model_select = bind(select_rows),
                settings = list(testid = tid, randseed = randseed, sim_key = sim_key,
                                mode = mode, ref_source = ref_source,
                                normalize = normalize, blacklist = blacklist,
                                species_key = species_key, n_genes = length(gene_cols),
                                mask_draw = mask_draw,
                                mask_nodes = if (is.null(clade_list)) NA_integer_ else
                                    sapply(clade_list, `[[`, 'node'),
                                data = if (mode != 'data') NA_character_ else
                                    if (is.data.frame(data)) '<data.frame>' else data,
                                regimes = regimes, arms = arms,
                                random_frac = if (is.null(random_frac)) NA_real_ else random_frac,
                                n_replicates = n_replicates,
                                clade_min = clade_min, clade_max = clade_max,
                                n_tips = n_tips, fixed_root = fixed_root,
                                lambda1 = lambda1, lambda2 = lambda2))

    stamp <- format(Sys.Date(), "%Y%m%d")
    saveRDS(out, sprintf('%s/%s_%s_dropout.rds', results_dir, tid, stamp))
    for (nm in c('clades', 'predictions', 'metrics', 'params', 'model_select')) {
        if (!is.null(out[[nm]])) {
            write.csv(out[[nm]], sprintf('%s/%s_%s_%s.csv', results_dir, tid, stamp, nm), row.names = FALSE)
        }
    }
    log_message(sprintf('Dropout validation complete. %d metric rows written.',
                        if (is.null(out$metrics)) 0 else nrow(out$metrics)), logfile, verbose = TRUE)

    return(out)
}
