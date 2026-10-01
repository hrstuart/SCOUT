
####### MODEL-CLASS IDENTIFIABILITY SWEEP #######
#
# How well SCOUT can tell BM1 / OU1 / OUM apart is not a fixed property of the method: it depends
# on the strength of selection relative to the height of the tree, and on the shape of the tree's
# branch lengths. Two distinct signals separate BM from OU.
#
#   1. Correlation shape. Under BM the tip covariance is sigma^2 * s_ij (shared root-to-MRCA
#      time); under OU it decays as exp(-alpha * d_ij). As alpha*H -> 0 the OU covariance tends to
#      the BM covariance and the two models become *mathematically* indistinguishable -- no method
#      can classify in that limit. As alpha*H grows the OU covariance tends to a diagonal, which
#      BM can never reproduce, and classification becomes easy.
#   2. Tip variance. Under BM the variance of a tip is sigma^2 * (its root-to-tip depth), so on a
#      tree whose tips sit at different depths the marginal variances vary; under OU the
#      stationary variance is the same for every tip. On an ultrametric tree this signal is
#      absent entirely and only signal 1 is available.
#
# These functions sweep alpha*H across tree variants built from one fixed topology, so tree shape
# is never confounded with tree height, and report both the empirical classification accuracy and
# an analytic identifiability reference that requires no model fitting.


######## TREE VARIANTS ########

#' Tip-indexed covariance inputs for a tree
#' Mirrors what preprocessTree assembles for parsed_alt_tree, without the regime machinery, so
#' compute_VCV can be called on an arbitrary tree.
#' @param phy Tree of class phylo.
#' @return List with n_leaves, leaf_dists and shared_lengths.
vcv_inputs <- function(phy) {
    n <- length(phy$tip.label)
    root_dists <- node.depth.edgelength(phy)
    list(n_leaves = n,
         leaf_dists = root_dists[seq_len(n)],
         shared_lengths = matrix(root_dists[mrca(phy)[seq_len(n), seq_len(n)]], n, n))
}

#' Build a branch-length variant of a tree, holding the topology fixed
#'
#' Produces the arms of the tree-structure factorial. Because the topology never changes, any
#' difference in downstream results is attributable to branch-length structure alone. Rescaling to
#' a common `height` is what makes alpha*H comparable across arms.
#'
#' @param phy Tree of class phylo.
#' @param branch `'keep'` retains the original branch lengths; `'unit'` sets every edge to 1.
#' @param depths `'keep'` leaves root-to-tip depths alone; `'equalise'` pads terminal branches so
#'   every tip sits at the same depth (ultrametric); `'vary'` equalises and then extends each
#'   terminal branch by a uniform draw, giving a spread of tip depths.
#' @param height Rescale all edge lengths so the tree height equals this. NULL leaves the scale
#'   alone.
#' @param depth_spread For `depths = 'vary'`, the maximum proportional extension of a tip's depth.
#'   The default 0.5 gives depths spanning a 1.5x range, matching the casNJ trees used earlier
#'   (tip depths 6-9 on a height-9 tree).
#' @param randseed Optional seed, used only by `depths = 'vary'`.
#' @return A phylo object.
#'
#' @details On a non-balanced topology "every edge equals 1" and "ultrametric" cannot both hold.
#'   `branch = 'unit'` with `depths = 'equalise'` therefore yields unit *internal* branches with
#'   padded terminal branches; the structural contrast that matters is uniform versus
#'   heterogeneous internal branch lengths, and that is preserved.
#' @import ape
#' @export
tree_variant <- function(phy, branch = c('keep', 'unit'), depths = c('keep', 'equalise', 'vary'),
                         height = NULL, depth_spread = 0.5, randseed = NULL) {
    branch <- match.arg(branch)
    depths <- match.arg(depths)

    if (is.null(phy$edge.length)) {
        phy$edge.length <- rep(1, Nedge(phy))
    }
    if (branch == 'unit') {
        phy$edge.length <- rep(1, Nedge(phy))
    }

    n <- length(phy$tip.label)
    term_row <- match(seq_len(n), phy$edge[, 2])   # the edge subtending each tip

    if (depths %in% c('equalise', 'vary')) {
        d <- node.depth.edgelength(phy)[seq_len(n)]
        # Extend rather than shrink, so a terminal branch can never be driven negative.
        phy$edge.length[term_row] <- phy$edge.length[term_row] + (max(d) - d)
    }
    if (depths == 'vary') {
        if (!is.null(randseed)) set.seed(randseed)
        H0 <- max(node.depth.edgelength(phy)[seq_len(n)])
        phy$edge.length[term_row] <- phy$edge.length[term_row] + runif(n, 0, depth_spread * H0)
    }

    if (!is.null(height)) {
        H0 <- max(node.depth.edgelength(phy)[seq_len(n)])
        phy$edge.length <- phy$edge.length * (height / H0)
    }

    return(phy)
}

#' The four standard arms of the tree-structure factorial
#' @param phy Tree of class phylo, supplying the topology.
#' @param height Common tree height for every arm.
#' @param depth_spread Passed to tree_variant.
#' @param randseed Optional seed.
#' @return Named list of phylo objects.
#' @export
tree_variant_set <- function(phy, height = 7, depth_spread = 0.5, randseed = NULL) {
    list(
        # as generated: heterogeneous branch lengths, all tips at one depth
        yule_ultra = tree_variant(phy, 'keep', 'keep', height, depth_spread, randseed),
        # TedSim-like: uniform internal branches, ultrametric
        unit_ultra = tree_variant(phy, 'unit', 'equalise', height, depth_spread, randseed),
        # casNJ-like: uniform branches, tips at varying depths
        unit_vary  = tree_variant(phy, 'unit', 'vary', height, depth_spread, randseed),
        # heterogeneous branches with varying tip depths
        yule_vary  = tree_variant(phy, 'keep', 'vary', height, depth_spread, randseed))
}


######## ANALYTIC IDENTIFIABILITY ########

#' How distinguishable are BM and OU on a given tree, before any fitting?
#'
#' Compares the BM and OU tip covariance matrices directly. Because this involves no data and no
#' optimisation it gives a clean floor: where the two covariances coincide, classification is
#' impossible for any method, and an empirical accuracy curve that does not track this is a sign
#' the harness is wrong rather than the model.
#'
#' @param phy Tree of class phylo.
#' @param alphaH Vector of alpha*H values to evaluate.
#' @param sigma Diffusion variance (already squared, as elsewhere in SCOUT).
#' @param add.root Root treatment for the OU covariance. BM always uses the fixed-root form,
#'   matching preprocessTree, which forces `root.fixed = TRUE` for BM1.
#' @return Data frame with one row per alpha*H:
#'   \describe{
#'     \item{cor_shape}{Correlation between the off-diagonal entries of the two *correlation*
#'       matrices. Tends to 1 as alpha*H -> 0, i.e. OU and BM predict the same correlation
#'       structure and are unidentifiable.}
#'     \item{resid_sd}{The part of the OU correlation structure that no affine rescaling of the BM
#'       correlation structure can reproduce, i.e. `sd(cor_ou) * sqrt(1 - cor_shape^2)`. Zero when
#'       the two are interchangeable and growing as they separate. Well behaved across the whole
#'       range, unlike a Gaussian divergence, which is dominated by the near-singularity of the
#'       small-alpha covariance rather than by anything that limits real inference.}
#'     \item{cv_var_bm, cv_var_ou}{Coefficient of variation of the diagonal. On an ultrametric
#'       tree both are 0 and the tip-variance signal is unavailable; where tips sit at different
#'       depths cv_var_bm > 0 while cv_var_ou stays at 0, which is an extra discriminator.}
#'   }
#'
#' @details The `alpha -> 0` limit depends on the root treatment. With `add.root = TRUE` (fixed
#'   root) the OU covariance converges exactly to the BM covariance and the models become strictly
#'   unidentifiable. With `add.root = FALSE` (stationary root, SCOUT's default when
#'   `fixed_root = FALSE`) it converges instead to a uniformly correlated matrix; that is not BM,
#'   but at small alpha both are near-affine in the shared time `s_ij`, which is why `cor_shape`
#'   still approaches 1 and the two remain hard to separate in practice.
#' @export
identifiability_reference <- function(phy, alphaH = c(0.25, 0.5, 1, 2, 5, 10, 20, 30, 50),
                                      sigma = 1, add.root = FALSE) {
    info <- vcv_inputs(phy)
    H <- max(info$leaf_dists)

    V_bm <- compute_VCV(info, 1e-10, sigma, add.root = TRUE)
    cv <- function(x) if (mean(x) == 0) 0 else stats::sd(x) / mean(x)
    ut <- upper.tri(V_bm)
    cor_bm <- stats::cov2cor(V_bm)[ut]

    out <- lapply(alphaH, function(ah) {
        alpha <- ah / H
        V_ou <- compute_VCV(info, alpha, sigma, add.root = add.root)
        cor_ou <- stats::cov2cor(V_ou)[ut]
        r <- stats::cor(cor_bm, cor_ou)
        data.frame(alphaH = ah, alpha = alpha, height = H,
                   cor_shape = r,
                   resid_sd = stats::sd(cor_ou) * sqrt(max(0, 1 - r^2)),
                   cv_var_bm = cv(diag(V_bm)),
                   cv_var_ou = cv(diag(V_ou)))
    })
    return(do.call(rbind, out))
}


######## EMPIRICAL SWEEP ########

#' Simulate on a tree at one alpha, fit all model classes, and score the classification
#' @param phy Tree with `states` and `node.label` set.
#' @param alpha Selection strength.
#' @param matrix_type `'evf'` scores the latent OU factors, `'counts'` the beta-Poisson counts.
#' @param sigma,t0,theta_step,ngenes,nevfs Simulation settings.
#' @param regimes,lambda1,lambda2,fixed_root,cores Fitting settings passed through to runSCOUT.
#' @param outdir Directory for the simulator's intermediate CSVs.
#' @param verbose Whether to log.
#' @return List with `selected` (named vector of chosen model per gene), `truth`, and `params`.
classification_accuracy <- function(phy, alpha, matrix_type = c('evf', 'counts'),
                                    sigma = 1, t0 = 8, theta_step = 2, ngenes = 30, nevfs = 30,
                                    regimes = c('BM1', 'OU1', 'OUM'),
                                    lambda1 = 0.2, lambda2 = 0.2, fixed_root = FALSE,
                                    outdir = tempdir(), cores = 1, verbose = FALSE) {
    matrix_type <- match.arg(matrix_type)

    sim <- simulate_test_data(ngenes = ngenes, ncells = length(phy$tip.label), tree = phy,
                              outdir = outdir, out_prefix = paste0('idsweep_', matrix_type),
                              a = alpha, s = sigma, t0 = t0, theta_step = theta_step,
                              nevfs = nevfs, sim_lin = FALSE)

    key <- names(sim$simulated_OU$evfs)[1]
    # The EVF matrix is the latent OU process itself and contains negative values, so it must not
    # be log-normalised. The counts matrix reaches tens of thousands, far above the theta upper
    # bound of 1000, so it must be.
    if (matrix_type == 'evf') {
        dat <- sim$simulated_OU$evfs[[key]]
        normalize <- FALSE
    } else {
        dat <- sim$simulated_OU$counts[[key]]
        normalize <- TRUE
    }
    rownames(dat) <- dat$species

    gcols <- setdiff(colnames(dat), c('OUM', 'species'))
    min_trait <- min(dat[, gcols], na.rm = TRUE)
    if (matrix_type == 'evf' && min_trait <= 0) {
        warning(sprintf('Simulated traits reach %.3f; theta is bounded below at 1e-10, so raise t0.', min_trait))
    }

    idata <- formatSCOUT(tree_path = phy, metadata_path = dat, outpath = outdir,
                         species_key = 'species', regimes = regimes,
                         anc_infer = 'ape', normalize = normalize)
    fit <- runSCOUT(idata, lambda1 = lambda1, lambda2 = lambda2, fixed.root = fixed_root,
                    scaleHeight = FALSE, cores = cores, verbose = FALSE)

    fits <- index_fits(fit)
    ntips <- length(phy$tip.label)
    genes <- idata$gene_cols

    # Selection is by AIC, matching SCOUT()'s filter(delta_AIC == 0) in main.R.
    selected <- vapply(genes, function(g) {
        aics <- vapply(regimes, function(r) fit_aic(fits[[paste(g, r, sep = '||')]], ntips)$aic, numeric(1))
        if (all(is.na(aics))) NA_character_ else regimes[which.min(aics)]
    }, character(1))

    params <- do.call(rbind, lapply(fit, function(x) data.frame(
        gene = x$settings[[1]], regime = x$settings[[2]],
        alpha = x$paras$alpha %||% NA_real_, sigma = x$paras$sigma %||% NA_real_,
        tau = x$paras$tau %||% NA_real_,
        converge = if (is.null(x$final$converge)) NA_character_ else as.character(x$final$converge),
        stringsAsFactors = FALSE)))

    return(list(selected = selected,
                truth = str_extract(genes, 'BM1|OU1|OUM'),
                min_trait = min_trait,
                params = params))
}

#' Sweep model-class identifiability over alpha*H and tree shape
#'
#' Simulates and fits across a grid of alpha*H values and branch-length variants of one topology,
#' and returns classification accuracy alongside an analytic identifiability reference. The grid
#' is parameterised by **alpha*H rather than alpha** because that is the quantity that governs
#' identifiability; every tree arm is rescaled to a common height so the two are interchangeable
#' within a run.
#'
#' @param tree Path to a newick file, or a phylo object, supplying the topology and tip states.
#' @param results_dir Directory for output files.
#' @param states Named vector of tip regimes. Optional if the tree carries `tree$states`.
#' @param alphaH Grid of alpha*H values.
#' @param arms Which tree variants to run; see tree_variant_set.
#' @param matrices Which data matrices to score: `'evf'`, `'counts'`, or both.
#' @param height Common tree height imposed on every arm.
#' @param sigma Diffusion variance.
#' @param t0 Root optimum. When NULL it is set to 3 * sigma * sqrt(height) rounded up, which keeps
#'   Brownian trajectories positive; theta cannot be fitted below 1e-10.
#' @param theta_step Spacing between regime optima.
#' @param ngenes,nevfs Genes for the counts matrix and latent factors for the EVF matrix.
#' @param n_replicates Independent simulations per grid cell.
#' @param depth_spread Passed to tree_variant.
#' @param regimes,lambda1,lambda2,fixed_root Fitting settings.
#' @param randseed Seed set once at the start.
#' @param testid Run ID; generated when NULL.
#' @param cores Cores for runSCOUT.
#' @param logfile Optional log file.
#' @param verbose Whether to log progress.
#' @return List with `accuracy` (per cell, with binomial CIs), `confusion`, `reference`
#'   (the analytic floor per arm), `params` and `settings`. All are written to `results_dir`.
#' @import ape
#' @import stringr
#' @export
runSCOUT.identifiability <- function(tree,
    results_dir,
    states = NULL,
    alphaH = c(0.5, 1, 2, 5, 10, 20, 30),
    arms = c('yule_ultra', 'unit_ultra', 'unit_vary', 'yule_vary'),
    matrices = c('evf', 'counts'),
    height = 7,
    sigma = 1,
    t0 = NULL,
    theta_step = 2,
    ngenes = 30,
    nevfs = 30,
    n_replicates = 2,
    depth_spread = 0.5,
    regimes = c('BM1', 'OU1', 'OUM'),
    lambda1 = 0.2,
    lambda2 = 0.2,
    fixed_root = FALSE,
    randseed = NULL,
    testid = NULL,
    cores = 1,
    logfile = NULL,
    verbose = TRUE) {

    if (!is.null(randseed)) set.seed(randseed)
    matrices <- match.arg(matrices, several.ok = TRUE)
    if (is.null(t0)) t0 <- ceiling(3 * sqrt(sigma) * sqrt(height))

    create_directory_if_not_exists(results_dir)
    tid <- if (is.null(testid)) {
        paste0('SCOUTID_', paste(sample(c(0:9, letters, LETTERS), 8, replace = TRUE), collapse = ""))
    } else testid
    log_message(sprintf('Identifiability sweep run ID = %s', tid), logfile, verbose = TRUE)
    log_message(sprintf('height = %g | sigma = %g | t0 = %g | %d alpha*H x %d arms x %d matrices x %d reps',
                        height, sigma, t0, length(alphaH), length(arms), length(matrices), n_replicates),
                logfile, verbose = verbose)

    phy0 <- prep_dropout_tree(tree, states, logfile, verbose)
    variants <- tree_variant_set(phy0, height = height, depth_spread = depth_spread,
                                 randseed = randseed)
    bad <- setdiff(arms, names(variants))
    if (length(bad) > 0) stop(sprintf('Unknown arm(s): %s', paste(bad, collapse = ', ')))
    variants <- variants[arms]

    sim_dir <- file.path(tempdir(), paste0('idsweep_', tid))
    create_directory_if_not_exists(sim_dir, message = FALSE)

    ## analytic reference -- no fitting, so evaluate on a denser grid
    dense <- sort(unique(c(alphaH, exp(seq(log(0.1), log(60), length.out = 25)))))
    reference <- do.call(rbind, lapply(names(variants), function(arm) {
        r <- identifiability_reference(variants[[arm]], dense, sigma = sigma, add.root = fixed_root)
        cbind(data.frame(arm = arm, stringsAsFactors = FALSE), r)
    }))

    ## empirical sweep
    acc_rows <- list(); conf_rows <- list(); par_rows <- list()
    for (arm in names(variants)) {
        phy <- variants[[arm]]
        for (ah in alphaH) {
            alpha <- ah / height
            for (mt in matrices) {
                for (rep in seq_len(n_replicates)) {
                    log_message(sprintf('arm %s | alpha*H %.2f | %s | rep %d', arm, ah, mt, rep),
                                logfile, verbose = verbose)
                    cell <- tryCatch(
                        classification_accuracy(phy, alpha, matrix_type = mt, sigma = sigma,
                                                t0 = t0, theta_step = theta_step,
                                                ngenes = ngenes, nevfs = nevfs, regimes = regimes,
                                                lambda1 = lambda1, lambda2 = lambda2,
                                                fixed_root = fixed_root, outdir = sim_dir,
                                                cores = cores, verbose = FALSE),
                        error = function(err) {
                            log_message(sprintf('  FAILED: %s', err$message), logfile, verbose = TRUE)
                            NULL
                        })
                    if (is.null(cell)) next

                    tag <- data.frame(arm = arm, alphaH = ah, alpha = alpha, matrix = mt,
                                      replicate = rep, stringsAsFactors = FALSE)
                    df <- data.frame(dataset = sprintf('%s|%g|%s|%d', arm, ah, mt, rep),
                                     model = unname(cell$selected), truth = cell$truth,
                                     stringsAsFactors = FALSE)
                    acc_rows[[length(acc_rows) + 1]] <- cbind(tag, calculate_group_class_accuracy(df))
                    conf_rows[[length(conf_rows) + 1]] <- cbind(tag, as.data.frame(
                        table(truth = df$truth, selected = df$model)))
                    par_rows[[length(par_rows) + 1]] <- cbind(tag, cell$params)
                }
            }
        }
    }

    bind <- function(x) if (length(x) == 0) NULL else do.call(rbind, x)
    out <- list(accuracy = bind(acc_rows), confusion = bind(conf_rows),
                reference = reference, params = bind(par_rows),
                trees = variants,
                settings = list(testid = tid, randseed = randseed, alphaH = alphaH, arms = arms,
                                matrices = matrices, height = height, sigma = sigma, t0 = t0,
                                theta_step = theta_step, ngenes = ngenes, nevfs = nevfs,
                                n_replicates = n_replicates, fixed_root = fixed_root))

    stamp <- format(Sys.Date(), "%Y%m%d")
    saveRDS(out, sprintf('%s/%s_%s_identifiability.rds', results_dir, tid, stamp))
    for (nm in c('accuracy', 'confusion', 'reference', 'params')) {
        if (!is.null(out[[nm]])) {
            write.csv(out[[nm]], sprintf('%s/%s_%s_%s.csv', results_dir, tid, stamp, nm),
                      row.names = FALSE)
        }
    }
    unlink(sim_dir, recursive = TRUE)
    log_message('Identifiability sweep complete.', logfile, verbose = TRUE)

    return(out)
}
