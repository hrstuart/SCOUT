# scout_perturb_lib.R -- tree perturbation (shuffle / nni / collapse) and classification scoring
# for the tree-robustness analyses. Uses only SCOUT's API; does not modify the package.
#
# Usage:
#   source('scout_perturb_lib.R')
#   pt  <- perturb_tree(phy, type = 'nni', intensity = 0.1, randseed = 1)
#   res <- score_classification(pt, dat, dataset_id = 'x', regimes = ...)

suppressPackageStartupMessages({
    library(SCOUT); library(ape); library(dplyr); library(stringr)
})

MODEL_PATTERN <- 'BM1|OU1|OUM'


# ---------------------------------------------------------------------------------------------
# Tree perturbation
# ---------------------------------------------------------------------------------------------

#' Number of internal (non-terminal) edges.
n_internal_edges <- function(phy) sum(phy$edge[, 2] > length(phy$tip.label))


#' Resolve polytomies and floor zero-length edges, so the saved tree is the fitted tree.
resolve_and_floor <- function(phy, tol = 1e-9) {
    if (!ape::is.binary(phy)) phy <- ape::multi2di(phy)
    if (is.null(phy$edge.length)) phy$edge.length <- rep(1, ape::Nedge(phy))
    phy$edge.length[phy$edge.length <= tol] <- 1e-7
    phy
}


#' Degrade a lineage tree in one controlled way.
#'
#' `shuffle` breaks the association between cells and tree positions without touching the tree
#' itself, giving a no-signal floor. `nni` rearranges local topology while keeping every branch
#' length, isolating topological error. `collapse` deletes internal edges, so both the topology and
#' the root-to-tip path lengths change.
#'
#' intensity is a fraction of tips (shuffle) or of internal edges (nni, collapse), rounded to the
#' nearest whole unit with a floor of 1 whenever intensity > 0 -- a count of 0 would silently turn an
#' arm into a no-op.
perturb_tree <- function(phy, type = c('shuffle', 'nni', 'collapse'),
                         intensity = 0.01, randseed = NULL) {
    type <- match.arg(type)
    if (!inherits(phy, 'phylo')) stop('phy must be a phylo object.')
    if (intensity < 0 || intensity > 1) stop('intensity must lie in [0, 1].')
    if (!is.null(randseed)) set.seed(randseed)

    if (is.null(phy$edge.length)) phy$edge.length <- rep(1, ape::Nedge(phy))
    unit_n <- if (type == 'shuffle') length(phy$tip.label) else n_internal_edges(phy)
    k <- if (intensity == 0) 0L else max(1L, as.integer(round(intensity * unit_n)))

    out <- switch(type,
        shuffle  = perturb_shuffle(phy, k),
        nni      = perturb_nni(phy, k),
        collapse = perturb_collapse(phy, k))

    out <- resolve_and_floor(out)
    attr(out, 'perturb_type') <- type
    attr(out, 'perturb_intensity') <- intensity
    attr(out, 'perturb_k') <- k
    out
}


#' Permute tip labels among themselves.
#'
#' Uses a random cyclic shift, which is a derangement for any k >= 2, so no selected tip keeps its
#' original label. A plain sample() leaves roughly one tip in place per permutation, which would
#' quietly weaken the no-signal floor.
perturb_shuffle <- function(phy, k) {
    if (k < 2) return(phy)
    idx <- sample(seq_along(phy$tip.label), k)
    shift <- sample.int(k - 1, 1)
    phy$tip.label[idx] <- phy$tip.label[idx[(seq_len(k) + shift - 1L) %% k + 1L]]
    phy
}


#' Apply k random nearest-neighbour interchanges, preserving branch lengths.
perturb_nni <- function(phy, k) {
    if (k < 1) return(phy)
    if (!requireNamespace('phangorn', quietly = TRUE)) {
        stop('perturb_tree(type = "nni") requires the phangorn package.')
    }
    phangorn::rNNI(phy, moves = k)
}


#' Delete k internal edges, creating polytomies that resolve_and_floor then re-resolves.
#'
#' The arbitrary re-resolution is part of the perturbation, not an artefact -- it is what a
#' reconstruction method does when it cannot resolve a node.
perturb_collapse <- function(phy, k) {
    if (k < 1) return(phy)
    internal <- which(phy$edge[, 2] > length(phy$tip.label))
    phy$edge.length[sample(internal, min(k, length(internal)))] <- 0
    ape::di2multi(phy, tol = 1e-9)
}


#' Generate and save a set of perturbed trees.
#'
#' Each replicate gets its own seed derived from randseed, so any single tree can be regenerated
#' without rerunning the set, and the set is reproducible across machines.
#'
#' Returns a manifest with `rf_dist` (Robinson-Foulds distance to baseline) as a direct check that
#' the perturbation did what it claims. Note RF is defined over label-bearing splits, so the shuffle
#' arm registers a large RF even though it leaves the unlabelled tree shape untouched -- that is the
#' correct reading, since shuffling is precisely a change in which cell sits at which tip.
perturb_tree_set <- function(phy, type, intensity, n, randseed = 1, outdir, prefix = NULL) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    tag <- if (is.null(prefix)) sprintf('%s_i%03d', type, round(intensity * 100)) else prefix
    has_phangorn <- requireNamespace('phangorn', quietly = TRUE)

    rows <- lapply(seq_len(n), function(r) {
        seed <- randseed + r
        pt <- perturb_tree(phy, type = type, intensity = intensity, randseed = seed)
        f <- file.path(outdir, sprintf('%s_rep%03d.nwk', tag, r))
        ape::write.tree(pt, f)
        rf <- if (has_phangorn) suppressWarnings(as.numeric(phangorn::RF.dist(phy, pt))) else NA_real_
        data.frame(arm = type, intensity = intensity, replicate = r, seed = seed,
                   tree_file = f, k = attr(pt, 'perturb_k'), rf_dist = rf,
                   stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}


# ---------------------------------------------------------------------------------------------
# Fit and score
# ---------------------------------------------------------------------------------------------

#' Fit every regime on one tree/matrix pair and return the AIC-selected model per gene.
#'
#' Mirrors what SCOUT() does in main.R:105-121 -- extract the grid-search history, annotate it (which
#' computes AIC = -2*ll + 2*param.count and the truth label from the gene name), and keep the
#' delta_AIC == 0 row per gene. Calling runSCOUT() directly rather than SCOUT() deliberately skips the
#' ~14 MB per-run .rds; across 650 runs that would be ~9 GB of files nothing downstream reads.
#'
#' Exact AIC ties are broken toward the model with fewer parameters, so a tie can never be scored as
#' a correct call for the more complex model by accident.
score_classification <- function(phy, dat, dataset_id,
                                 regimes = c('BM1', 'OU1', 'OUM'),
                                 normalize = TRUE,
                                 lambda1 = 0.2, lambda2 = 0.2, fixed_root = FALSE,
                                 cores = 1, outdir = tempdir(), logfile = NULL) {

    idata <- formatSCOUT(tree_path = phy, metadata_path = dat, outpath = outdir,
                         species_key = 'species', regimes = regimes,
                         anc_infer = 'ape', normalize = normalize, logfile = logfile)

    fit <- runSCOUT(idata, lambda1 = lambda1, lambda2 = lambda2, fixed.root = fixed_root,
                    scaleHeight = FALSE, cores = cores, logfile = logfile, verbose = FALSE)

    history <- SCOUT:::extract_history_grid_search(fit)
    if (!'converge' %in% colnames(history)) history$converge <- NA
    history$dataset <- dataset_id
    history$ntips <- length(phy$tip.label)

    annotated <- annotate_history(history, datasetid = 'dataset')

    # Model selection column, across two SCOUT builds.
    #
    # The SCOUT installed on 2026-08-18 02:06 UTC changed annotate_history: it still computes
    # delta_AIC and still derives best_fit = (delta_AIC == 0) from it, but it now DROPS delta_AIC
    # and delta_AICc before returning (analysis_utils.R:47, the trailing
    # `select(-c(delta_AIC, delta_AICc))`). The previous build (2026-08-16 22:51 UTC, kept at
    # software/R/SCOUT_installed_backup_20260817_prebuild), which produced every array result
    # currently in output/, kept delta_AIC and had no best_fit column at all.
    #
    # A bare `filter(delta_AIC == 0)` therefore fails on the current build with
    # "object 'delta_AIC' not found", taking down the fit AFTER runSCOUT has already run -- about
    # 31 minutes of wasted compute per condition. Accept either column; they are the same
    # predicate by construction, and both are computed under the same group_by(dataset, gene_name).
    if ('best_fit' %in% colnames(annotated)) {
        picked <- annotated %>% filter(best_fit)
    } else if ('delta_AIC' %in% colnames(annotated)) {
        picked <- annotated %>% filter(delta_AIC == 0)
    } else {
        stop('annotate_history() returned neither `best_fit` nor `delta_AIC`; cannot select a ',
             'model per gene. Check the installed SCOUT build.')
    }

    best <- picked %>%
        group_by(gene_name) %>%
        arrange(param.count, .by_group = TRUE) %>%
        slice_head(n = 1) %>%
        ungroup()

    selected <- data.frame(dataset = dataset_id,
                           gene_name = best$gene_name,
                           model = best$model,
                           truth = str_extract(best$gene_name, MODEL_PATTERN),
                           AIC = best$AIC,
                           loglik = best$ll_total,
                           alpha = best$alpha,
                           sigma = best$sigma,
                           tau = best$tau,
                           converge = as.character(best$converge),
                           stringsAsFactors = FALSE)

    list(selected = selected, annotated = annotated,
         n_genes = length(idata$gene_cols),
         min_trait = min(dat[, idata$gene_cols], na.rm = TRUE))
}


#' Accuracy summary for one fitted dataset.
#'
#' Wraps SCOUT's calculate_group_class_accuracy(), which returns overall accuracy with a binomial 95%
#' CI plus per-class sensitivity/specificity. The Overall row is the summary the prompt asks for.
accuracy_summary <- function(selected, tag = NULL) {
    acc <- calculate_group_class_accuracy(selected, dataset_col = 'dataset',
                                          model_col = 'model', truth_col = 'truth')
    if (!is.null(tag)) acc <- cbind(tag[rep(1, nrow(acc)), , drop = FALSE], acc, row.names = NULL)
    acc
}


#' Confusion table (truth x selected) for one fitted dataset.
confusion_table <- function(selected, tag = NULL) {
    tab <- as.data.frame(table(truth = selected$truth, selected = selected$model))
    if (!is.null(tag)) tab <- cbind(tag[rep(1, nrow(tab)), , drop = FALSE], tab, row.names = NULL)
    tab
}


#' Write the three per-task output files for one fitting unit.
write_task_output <- function(res, tag, results_dir, stem) {
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(accuracy_summary(res$selected, tag),
              file.path(results_dir, sprintf('%s_accuracy.csv', stem)), row.names = FALSE)
    write.csv(confusion_table(res$selected, tag),
              file.path(results_dir, sprintf('%s_confusion.csv', stem)), row.names = FALSE)
    write.csv(cbind(tag[rep(1, nrow(res$selected)), , drop = FALSE], res$selected, row.names = NULL),
              file.path(results_dir, sprintf('%s_selected.csv', stem)), row.names = FALSE)
    invisible(NULL)
}


#