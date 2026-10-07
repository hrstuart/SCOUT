#!/usr/bin/env Rscript
# ==============================================================================
# simulated_masking.R
#
# SCOUT dropout (held-out masking) validation on simulated beta-Poisson counts, swept over
# tree size (128-1024 tips) x alpha (0.5-3), with both a known latent trait and real
# observation noise.
#
# Masked-cell prediction (SCOUT::runSCOUT.dropout) is based on mvMORPH::estim()
# (Clavel, Escarguel & Merceron 2015, Methods Ecol. Evol. 6:1311-1319).
#
# Usage:
#   SCOUT_SMOKE=1 SCOUT_CELL=1 Rscript simulated_masking.R     # smoke test (~3-6 min)
#   SCOUT_CELL=1 SCOUT_CORES=16 Rscript simulated_masking.R    # one grid cell (or SLURM_ARRAY_TASK_ID)
#   Rscript simulated_masking.R                                 # all 16 cells serially (days)
#
#   Grid: cells 1-4 = 128 tips, 5-8 = 256, 9-12 = 512, 13-16 = 1024;
#         alpha cycles 0.5 / 1 / 2 / 3 within each block.
# ==============================================================================

.libPaths(c(Sys.getenv('SCOUT_LIB',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/software/R/R-4.4.2/library'),
    .libPaths()))

# runSCOUT parallelises with future::multisession; leaving BLAS threads unpinned
# oversubscribes every worker. Same convention as the bandwidth counts sweep.
Sys.setenv(OMP_NUM_THREADS = '1')

suppressPackageStartupMessages({
    library(SCOUT)
    library(ape)
    library(castor)
    library(dplyr)
    library(tidyr)
})

t_start <- Sys.time()

## ---------------------------------------------------------------- config ----

BASE    <- Sys.getenv('SCOUT_BASE',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/2_dropout_output')
OUT_DIR <- file.path(BASE, 'output_counts_loo')       # new dir: v2 EVF outputs untouched
SIM_DIR <- file.path(BASE, 'simulated_data_counts')

SMOKE <- identical(Sys.getenv('SCOUT_SMOKE'), '1')

cfg <- list(
    sizes        = c(128, 256, 512, 1024),
    alphas       = c(0.5, 1, 2, 3),      # spans the requested 0.5-3; superset of the
                                         # 1_tree_robustness 0.5/1/3 grid; a=3 == v2
    sigma        = 1,
    t0           = 8,        # not the default 3: low t0 drives BM1 traits non-positive
    theta_step   = 2,
    nevfs        = 25,       # ngenes == nevfs is REQUIRED for identity effects -> 75 traits
    regimes      = c('BM1', 'OU1', 'OUM'),
    arms         = c('clade', 'random'),
    random_frac  = c(NA, 0.10),          # NA = size-matched to the clade
    clade_frac   = c(0.02, 0.08),        # 2%-8% of tips == v2's 5-20 of 256, held constant
    scale_s      = 10,
    burst_scale  = 1,
    signal_gain  = c(1, 0, 0),           # kon keeps the OU trait; koff and s collapse
    lambda1      = 0.2,
    lambda2      = 0.2,
    fixed_root   = FALSE,
    tau_prior_mean = 0.2,
    tau_prior_sd   = 0.1,
    tree_seed    = 1,
    randseed     = 42,
    cores        = as.integer(Sys.getenv('SCOUT_CORES', '16')),
    gate1_median_min = 0.70,   # abort the cell below this
    gate1_frac_min   = 0.60,   # fraction of genes required above 0.8
    do_calibration   = FALSE)  # isotonic calibration: exploratory, off by default

# Replicate budget tapers with tree size to bound wall time.
reps_for <- function(n) switch(as.character(n), '128' = 20L, '256' = 20L,
                                                '512' = 10L, '1024' = 6L, 10L)

if (SMOKE) {
    cfg <- modifyList(cfg, list(sizes = c(64), alphas = c(3), nevfs = 10,
                                random_frac = c(NA, 0.20), cores = 4))
    reps_for <- function(n) 2L
}

# Regime transition matrix for tip states -- identical to the v2 run.
P <- matrix(c(0.90, 0.08, 0.02,
              0.10, 0.80, 0.10,
              0.02, 0.08, 0.90), 3, 3, byrow = TRUE)

for (d in c(OUT_DIR, SIM_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

## The affine path in Get_params only engages when param_mode == 'affine' AND
## !is.na(signal_gain[i]) (sim_utilities.R:141). c(1, NA, NA) would SILENTLY fall
## back to the rank path for koff and s -- the exact opposite of the intent.
stopifnot(!any(is.na(cfg$signal_gain)), length(cfg$signal_gain) == 3)

## ------------------------------------------------------------------ grid ----

grid <- expand.grid(n_tips = cfg$sizes, alpha_true = cfg$alphas)
grid <- grid[order(grid$n_tips, grid$alpha_true), ]
rownames(grid) <- NULL
grid$cell <- seq_len(nrow(grid))

cell_env <- Sys.getenv('SCOUT_CELL', Sys.getenv('SLURM_ARRAY_TASK_ID', ''))
cells <- if (nzchar(cell_env)) as.integer(cell_env) else grid$cell
stopifnot(all(cells %in% grid$cell))

LOGFILE <- file.path(OUT_DIR, sprintf('counts_loo_cell%s_run.log',
                                      if (nzchar(cell_env)) cell_env else 'ALL'))
say <- function(...) {
    msg <- sprintf('[%s] %s', format(Sys.time(), '%H:%M:%S'), sprintf(...))
    cat(msg, '\n', sep = '')
    cat(msg, '\n', sep = '', file = LOGFILE, append = TRUE)
}

say('=== SCOUT counts LOO sweep%s ===', if (SMOKE) ' | SMOKE' else '')
say('SCOUT %s | R %s | %d cores of %d | cells: %s',
    as.character(packageVersion('SCOUT')), getRversion(), cfg$cores,
    parallel::detectCores(), paste(cells, collapse = ','))

# Preflight: fail now, not hours in. This script depends on two things the
# installed SCOUT must already have -- it installs nothing itself.
if (!'random_frac' %in% names(formals(SCOUT::runSCOUT.dropout))) {
    stop('Installed SCOUT has no `random_frac` argument on runSCOUT.dropout.')
}
if (!'param_mode' %in% names(formals(SCOUT:::Get_params))) {
    stop('Installed SCOUT has no `param_mode` argument on Get_params; the affine ',
         'path this analysis depends on is absent.')
}
if (!'signal_gain' %in% names(formals(SCOUT:::Get_params))) {
    stop('Installed SCOUT has no `signal_gain` argument on Get_params.')
}

## ============================================================================
## BLOCK 2 -- identity-effect counts simulator (spec changes 1 + 2)
## ============================================================================
# Mirrors generate_EvoCounts (simulate.R:117) and simulate_OU_genes
# (simulate.R:145) exactly, with ONE substitution: the random sparse gene-effect
# matrix from GeneEffects() is replaced by the identity. Because
# params[[i]] <- evf[[i]] %*% t(gene_effects[[i]]), identity effects make gene g
# a function of EVF g ALONE -- which is what makes a known latent truth possible.

match_param_densities <- function() {
    e <- new.env()
    utils::data('match_params', package = 'SCOUT', envir = e)
    mp <- get('match_params', envir = e)
    for (j in 1:3) mp[, j] <- log(base = 10, mp[, j])
    lapply(1:3, function(i) stats::density(mp[, i], n = 2000))
}

simulate_counts_identity <- function(tree, metadata, nevfs, a, s, t0, theta_step,
                                     nstates, scale_s, burst_scale, signal_gain) {
    ngenes <- nevfs                     # identity requires squareness
    ncells <- length(tree$tip.label)

    eevfs <- SCOUT:::generate_EvoEVFs(tree, metadata, nevfs, a, s, t0, theta_step, nstates)
    ge    <- rep(list(diag(nevfs)), 3)  # <-- the identity substitution
    dens  <- match_param_densities()

    kinetics <- list()
    ecounts <- lapply(seq_along(eevfs), function(m) {
        imod <- eevfs[[m]]
        params <- SCOUT:::Get_params(ge, imod, dens, bimod = 0, scale_s,
                                     'affine', signal_gain, burst_scale)
        kinetics[[m]] <<- params
        cnt <- lapply(seq_len(ngenes), function(i) {
            sapply(seq_len(ncells), function(j) {
                y <- rbeta(1, params[[1]][i, j], params[[2]][i, j])
                rpois(1, y * params[[3]][i, j])
            })
        })
        do.call(rbind, cnt)
    })

    counts <- t(do.call(rbind, ecounts))
    states <- metadata$cluster; names(states) <- metadata$cellID
    gnames <- c(paste0('BM1_', 1:ngenes), paste0('OU1_', 1:ngenes), paste0('OUM_', 1:ngenes))

    counts_df <- as.data.frame(counts)
    colnames(counts_df) <- gnames
    counts_df$OUM     <- states[tree$tip.label]
    counts_df$species <- tree$tip.label

    # The saved truth is the kon block only -- eevfs[[model]][[1]] -- which is
    # exactly the block signal_gain = c(1,0,0) keeps alive in the counts.
    ou_evfs <- cbind(eevfs[[1]][[1]], eevfs[[2]][[1]], eevfs[[3]][[1]])
    evf_df  <- as.data.frame(ou_evfs)
    colnames(evf_df) <- gnames
    evf_df$OUM     <- states[tree$tip.label]
    evf_df$species <- tree$tip.label

    list(counts = counts_df, evfs = evf_df, kinetics = kinetics, gene_names = gnames)
}

## ============================================================================
## BLOCK 3 -- Gates 1 and 2
## ============================================================================

run_gate1 <- function(sim, prefix) {
    cnt <- as.matrix(sim$counts[, sim$gene_names, drop = FALSE])
    z   <- as.matrix(sim$evfs[,  sim$gene_names, drop = FALSE])

    rho <- vapply(seq_along(sim$gene_names), function(g)
        suppressWarnings(stats::cor(log1p(cnt[, g]), z[, g], method = 'spearman')),
        numeric(1))

    # Scrambled control: shuffling the gene index must destroy the correlation.
    # Without this the gate could pass on a tautology.
    set.seed(99)
    perm <- sample(ncol(z))
    rho_scram <- vapply(seq_along(sim$gene_names), function(g)
        suppressWarnings(stats::cor(log1p(cnt[, g]), z[, perm[g]], method = 'spearman')),
        numeric(1))

    # koff and s must be per-gene constants under signal_gain = c(1,0,0).
    sd_koff <- max(vapply(sim$kinetics, function(p) max(apply(p[[2]], 1, stats::sd)), numeric(1)))
    sd_s    <- max(vapply(sim$kinetics, function(p) max(apply(p[[3]], 1, stats::sd)), numeric(1)))

    tab <- data.frame(gene = sim$gene_names, spearman = rho, spearman_scrambled = rho_scram)
    write.csv(tab, file.path(OUT_DIR, paste0(prefix, '_gate1.csv')), row.names = FALSE)

    med  <- stats::median(rho, na.rm = TRUE)
    frac <- mean(rho > 0.8, na.rm = TRUE)
    say('  Gate 1: median rho = %.3f | frac>0.8 = %.2f | scrambled median = %.3f',
        med, frac, stats::median(rho_scram, na.rm = TRUE))
    say('  Gate 1: max sd(koff) = %.3g | max sd(s) = %.3g (both must be ~0)', sd_koff, sd_s)

    ok <- med >= cfg$gate1_median_min && frac >= cfg$gate1_frac_min &&
          abs(stats::median(rho_scram, na.rm = TRUE)) < 0.3 &&
          sd_koff < 1e-8 && sd_s < 1e-8
    list(ok = ok, median = med, frac_above_0.8 = frac,
         scrambled_median = stats::median(rho_scram, na.rm = TRUE),
         sd_koff = sd_koff, sd_s = sd_s)
}

run_gate2 <- function(sim, prefix) {
    cnt <- as.matrix(sim$counts[, sim$gene_names, drop = FALSE])
    zf  <- mean(cnt == 0)
    gm  <- stats::median(colMeans(cnt))
    out <- data.frame(zero_fraction = zf, median_gene_mean = gm,
                      max_count = max(cnt), scale_s = cfg$scale_s,
                      burst_scale = cfg$burst_scale)
    write.csv(out, file.path(OUT_DIR, paste0(prefix, '_sparsity.csv')), row.names = FALSE)
    # v2 counts were 7.2% zeros vs ~83% for real C. elegans: this is OPTIMISTIC
    # relative to real scRNA-seq and must be reported as such. scale_s sparsifies.
    say('  Gate 2: zeros = %.1f%% | median gene mean = %.1f | max = %d  (v2 ref: 7.2%%)',
        100 * zf, gm, max(cnt))

    # CAVEAT on sparsifying. runSCOUT.dropout will warn "Simulated traits reach
    # 0.0000 ... raise t0". On COUNTS that warning is a false alarm in its stated
    # form -- log1p(0) is exactly 0, so every zero count trips it, and raising t0
    # does not remove count zeros. It does flag a real limit though: m_step_tipfog
    # bounds theta below at 1e-10, so a gene that is mostly zeros has its theta
    # pinned at the bound and is fit badly. At the default scale_s the median gene
    # mean is in the hundreds and this is harmless; it is the binding constraint on
    # how far scale_s can be lowered in pursuit of realistic sparsity.
    lg <- log1p(cnt)
    near_bound <- mean(colMeans(lg) < 0.1)
    if (near_bound > 0) {
        say('  Gate 2: WARNING %.1f%% of genes have mean log1p < 0.1 -- theta will pin at its',
            100 * near_bound)
        say('          lower bound (1e-10) for these. Raise scale_s, or drop them.')
    }
    out$frac_genes_near_theta_bound <- near_bound
    out
}

## ============================================================================
## BLOCK 4 -- tree construction
## ============================================================================

build_tree <- function(n) {
    set.seed(cfg$tree_seed)
    Q <- expm::logm(P)
    obj <- generate_tree_hbd_reverse(n, rho = 1, lambda = 1, mu = 0)
    phy <- obj$trees[[1]]
    phy$tip.label <- paste0('t', phy$tip.label)
    st  <- simulate_mk_model(phy, Q = Q, include_nodes = FALSE)
    phy$states <- as.factor(paste0('state', st$tip_states))
    phy <- infer_anc(phy)
    # Normalise through the same path runSCOUT.dropout uses internally, so the
    # tree the counts are simulated on is byte-for-byte the tree that gets fitted.
    ts  <- setNames(as.character(phy$states), phy$tip.label)
    phy <- prep_dropout_tree(phy, states = ts, verbose = FALSE)
    phy
}

## ============================================================================
## BLOCK 6 -- post-hoc latent scoring (spec changes 4 + 5)
## ============================================================================
# res$predictions carries z_true (= held-out log1p count, the OBSERVED scale) and
# z_masked / z_prior / z_oracle. The LATENT truth lives in the EVF matrix. Join on
# (tip, gene) and take Spearman: rank correlation is invariant to the monotone
# count transform, so it is directly comparable to the v2 EVF-scale Spearman
# (random arm 0.666, clade arm 0.262).

score_latent <- function(res, evf_df, gene_names, prefix) {
    if (is.null(res$predictions) || nrow(res$predictions) == 0) {
        say('  WARNING: no predictions retained; skipping latent scoring.')
        return(NULL)
    }
    z_long <- evf_df %>%
        select(species, all_of(gene_names)) %>%
        pivot_longer(-species, names_to = 'gene', values_to = 'z_latent')

    # INNER join: formatSCOUT drops all-zero-sum genes (SCOUT_EM_utils.R:1009),
    # so the counts run can carry fewer genes than the EVF truth.
    pred <- res$predictions %>%
        inner_join(z_long, by = c('tip' = 'species', 'gene' = 'gene'))

    n_missing <- length(setdiff(unique(res$predictions$gene), unique(pred$gene)))
    if (n_missing > 0) say('  NOTE: %d gene(s) present in predictions but not in truth.', n_missing)

    safe_cor <- function(a, b, m) {
        if (length(a) < 3 || stats::sd(a) == 0 || stats::sd(b) == 0) return(NA_real_)
        suppressWarnings(stats::cor(a, b, method = m))
    }

    # Per-mask. On a near-ultrametric tree V_ho is close to rank one, so within a
    # single masked CLADE the predictions are nearly constant and this is unstable
    # (dropout_validation.R:519-521). Hence the pooled version below.
    per_mask <- pred %>%
        group_by(replicate, arm, gene, regime, true_model, is_selected) %>%
        summarise(n_masked = n(),
                  spearman_latent_masked = safe_cor(z_masked, z_latent, 'spearman'),
                  spearman_latent_prior  = safe_cor(z_prior,  z_latent, 'spearman'),
                  spearman_latent_oracle = safe_cor(z_oracle, z_latent, 'spearman'),
                  pearson_latent_masked  = safe_cor(z_masked, z_latent, 'pearson'),
                  spearman_obs_masked    = safe_cor(z_masked, z_true,   'spearman'),
                  .groups = 'drop')

    # Pooled across replicates within an arm: the number to trust for the clade arm.
    pooled <- pred %>%
        group_by(arm, regime, true_model, is_selected) %>%
        summarise(n_points = n(),
                  spearman_latent_masked = safe_cor(z_masked, z_latent, 'spearman'),
                  spearman_latent_prior  = safe_cor(z_prior,  z_latent, 'spearman'),
                  spearman_latent_oracle = safe_cor(z_oracle, z_latent, 'spearman'),
                  pearson_latent_masked  = safe_cor(z_masked, z_latent, 'pearson'),
                  .groups = 'drop')

    write.csv(per_mask, file.path(OUT_DIR, paste0(prefix, '_latent_spearman.csv')),
              row.names = FALSE)
    write.csv(pooled, file.path(OUT_DIR, paste0(prefix, '_latent_spearman_pooled.csv')),
              row.names = FALSE)
    list(per_mask = per_mask, pooled = pooled)
}

## ============================================================================
## Per-cell driver
## ============================================================================

run_cell <- function(cell) {
    n     <- grid$n_tips[grid$cell == cell]
    alpha <- grid$alpha_true[grid$cell == cell]
    nrep  <- reps_for(n)
    prefix <- sprintf('counts_loo_n%d_a%s', n, format(alpha, trim = TRUE))
    if (SMOKE) prefix <- paste0(prefix, '_smoke')

    say('--- cell %d: n = %d, alpha = %s, reps = %d, prefix = %s ---',
        cell, n, format(alpha), nrep, prefix)

    ## tree
    phy <- build_tree(n)
    tip_states <- setNames(as.character(phy$states), phy$tip.label)
    stopifnot(Ntip(phy) == n, is.binary(phy), !is.null(phy$node.label))
    say('  Tree: %d tips, ultrametric = %s, regimes = %s', Ntip(phy), is.ultrametric(phy),
        paste(sprintf('%s:%d', names(table(tip_states)), table(tip_states)), collapse = ' '))

    # Clade sizes as a constant FRACTION of tips, so the clade arm means the same
    # thing at every tree size (v2's 5-20 of 256 is 2%-8%).
    clade_min <- max(3L, as.integer(round(cfg$clade_frac[1] * n)))
    clade_max <- max(clade_min + 2L, as.integer(round(cfg$clade_frac[2] * n)))
    elig <- eligible_clades(phy, clade_min, clade_max, states = tip_states,
                            exclude_root_children = TRUE)
    say('  Eligible clades of %d-%d tips: %d (need %d).', clade_min, clade_max, nrow(elig), nrep)
    if (nrow(elig) < nrep) {
        stop(sprintf('Only %d eligible clades for %d replicates at n = %d; widen clade_frac.',
                     nrow(elig), nrep, n))
    }

    write.tree(phy, file.path(SIM_DIR, paste0(prefix, '_tree.nwk')))

    ## simulate counts with identity gene effects
    say('  Simulating: nevfs = ngenes = %d, a = %s, s = %s, t0 = %s, signal_gain = c(%s)',
        cfg$nevfs, format(alpha), format(cfg$sigma), format(cfg$t0),
        paste(cfg$signal_gain, collapse = ','))
    set.seed(cfg$randseed)
    metadata <- data.frame(cellID = phy$tip.label, cluster = phy$states)
    nstates  <- length(unique(phy$states))
    sim <- simulate_counts_identity(phy, metadata, cfg$nevfs, alpha, cfg$sigma,
                                    cfg$t0, cfg$theta_step, nstates,
                                    cfg$scale_s, cfg$burst_scale, cfg$signal_gain)

    write.csv(sim$counts, file.path(SIM_DIR, paste0(prefix, '_counts.csv')), row.names = FALSE)
    write.csv(sim$evfs,   file.path(SIM_DIR, paste0(prefix, '_evf_truth.csv')), row.names = FALSE)

    ## gates
    g1 <- run_gate1(sim, prefix)
    g2 <- run_gate2(sim, prefix)
    if (!g1$ok) {
        say('  GATE 1 FAILED for cell %d -- the count/latent mapping is not usable. Skipping.', cell)
        return(invisible(list(cell = cell, gate1 = g1, status = 'gate1_failed')))
    }

    ## ------------------------------------------------------------------------
    ## BLOCK 5 -- the counts swap (spec change 4)
    ## runSCOUT.dropout reads its OBSERVED data from
    ## sim_res$simulated_OU$evfs[[sim_key]] (dropout_validation.R:756). Putting the
    ## COUNTS matrix in that slot is what makes this a counts run; the EVF matrix
    ## is held back as the scoring truth. normalize = TRUE applies log1p in
    ## formatSCOUT -- required, since counts reach ~1e5, far above the theta upper
    ## bound of 1000 (precedent: classification_accuracy, identifiability.R:204).
    ## ------------------------------------------------------------------------
    key <- sprintf('a%.2f_s%.2f', alpha, cfg$sigma)
    sim_res <- list(scLT = NULL,
                    simulated_OU = list(counts = list(), evfs = list()))
    sim_res$simulated_OU$evfs[[key]]   <- sim$counts   # <-- DELIBERATE
    sim_res$simulated_OU$counts[[key]] <- sim$counts

    say('  Launching %d reps x %d arms = %d masked fits, %d traits x %d regimes each.',
        nrep, 1 + length(cfg$random_frac), nrep * (1 + length(cfg$random_frac)),
        length(sim$gene_names), length(cfg$regimes))

    res <- runSCOUT.dropout(
        tree             = phy,
        results_dir      = OUT_DIR,
        sim_res          = sim_res,
        sim_key          = key,
        testid           = prefix,
        normalize        = TRUE,          # log1p the counts
        regimes          = cfg$regimes,
        n_replicates     = nrep,
        clade_min        = clade_min,
        clade_max        = clade_max,
        arms             = cfg$arms,
        random_frac      = cfg$random_frac,
        mask_draw        = 'upfront',     # sim mode defaults to 'inline', which repeats clades
        keep_predictions = 'true_model',
        lambda1          = cfg$lambda1,
        lambda2          = cfg$lambda2,
        fixed_root       = cfg$fixed_root,
        tau_prior_mean   = cfg$tau_prior_mean,
        tau_prior_sd     = cfg$tau_prior_sd,
        randseed         = cfg$randseed,
        cores            = cfg$cores,
        save_full_fit    = TRUE,
        checkpoint       = TRUE,
        logfile          = LOGFILE,
        verbose          = TRUE)

    ## post-hoc latent scoring
    lat <- score_latent(res, sim$evfs, sim$gene_names, prefix)

    ## counts-scale summaries (commensurable with the C. elegans run)
    arm_levels <- c('clade', 'random',
                    sprintf('random_%02.0fpct', 100 * cfg$random_frac[!is.na(cfg$random_frac)]))
    sel <- res$metrics %>% filter(is_selected)
    by_arm <- sel %>%
        mutate(arm = factor(arm, levels = intersect(arm_levels, unique(arm)))) %>%
        group_by(arm, predictor) %>%
        summarise(n_gene_masks = n(),
                  rmse = mean(rmse, na.rm = TRUE), mae = mean(mae, na.rm = TRUE),
                  pearson = mean(pearson, na.rm = TRUE),
                  spearman = mean(spearman, na.rm = TRUE),
                  coverage = mean(coverage, na.rm = TRUE),
                  coverage_obs = mean(coverage_obs, na.rm = TRUE),
                  mean_n_masked = mean(n_masked, na.rm = TRUE), .groups = 'drop') %>%
        arrange(arm, predictor)
    write.csv(by_arm, file.path(OUT_DIR, paste0(prefix, '_summary_by_arm.csv')), row.names = FALSE)

    model_recovery <- res$model_select %>%
        group_by(arm, true_model) %>%
        summarise(n = n(),
                  acc_masked = mean(best_masked == true_model, na.rm = TRUE),
                  acc_full   = mean(best_full   == true_model, na.rm = TRUE),
                  .groups = 'drop')
    write.csv(model_recovery, file.path(OUT_DIR, paste0(prefix, '_model_recovery.csv')),
              row.names = FALSE)

    mask_sizes <- res$clades %>% filter(!is.na(arm)) %>%
        mutate(frac_masked = n_masked / Ntip(phy))
    write.csv(mask_sizes, file.path(OUT_DIR, paste0(prefix, '_mask_sizes.csv')), row.names = FALSE)

    fail_cols <- intersect(c('n_fit_missing', 'n_regime_lost', 'n_prediction_failed',
                             'n_gene_dropped'), names(mask_sizes))
    fails <- colSums(mask_sizes[fail_cols], na.rm = TRUE)
    say('  Failure tallies: %s', paste(sprintf('%s=%d', names(fails), fails), collapse = ' '))
    if (any(fails > 0)) say('  WARNING: non-zero failure tallies; inspect the clades table.')

    ## manifest -- note the design columns are n_tips / alpha_true, NOT `alpha`.
    ## The params table already carries its own fitted `alpha`; a second column of
    ## that name would recreate the known duplicate-alpha trap in the sweep outputs.
    manifest <- list(prefix = prefix, cell = cell, n_tips = n, alpha_true = alpha,
                     sigma = cfg$sigma, t0 = cfg$t0, theta_step = cfg$theta_step,
                     nevfs = cfg$nevfs, ngenes = cfg$nevfs, n_traits = length(sim$gene_names),
                     n_replicates = nrep, arms = paste(cfg$arms, collapse = ','),
                     random_frac = paste(cfg$random_frac, collapse = ','),
                     clade_min = clade_min, clade_max = clade_max,
                     regimes = paste(cfg$regimes, collapse = ','),
                     scale_s = cfg$scale_s, burst_scale = cfg$burst_scale,
                     signal_gain = paste(cfg$signal_gain, collapse = ','),
                     param_mode = 'affine', effect_mode = 'identity',
                     score_target = 'counts', normalize = TRUE,
                     lambda1 = cfg$lambda1, lambda2 = cfg$lambda2,
                     tree_seed = cfg$tree_seed, randseed = cfg$randseed,
                     gate1_median = g1$median, gate1_frac_above_0.8 = g1$frac_above_0.8,
                     gate1_scrambled_median = g1$scrambled_median,
                     zero_fraction = g2$zero_fraction, median_gene_mean = g2$median_gene_mean,
                     sim_key = key, smoke = SMOKE,
                     scout_version = as.character(packageVersion('SCOUT')),
                     r_version = as.character(getRversion()),
                     script = 'scripts/260817_counts_loo_sweep.R')
    write.csv(data.frame(parameter = names(manifest),
                         value = vapply(manifest, function(x) paste(x, collapse = ','), character(1)),
                         row.names = NULL),
              file.path(OUT_DIR, paste0(prefix, '_manifest.csv')), row.names = FALSE)
    saveRDS(list(manifest = manifest, session = sessionInfo()),
            file.path(OUT_DIR, paste0(prefix, '_manifest.rds')))

    say('  Cell %d complete.', cell)
    invisible(list(cell = cell, n_tips = n, alpha_true = alpha, prefix = prefix,
                   by_arm = by_arm, latent = lat, gate1 = g1, gate2 = g2, status = 'ok'))
}

## ============================================================================
## Run
## ============================================================================

results <- list()
for (cl in cells) {
    results[[as.character(cl)]] <- tryCatch(run_cell(cl), error = function(e) {
        say('  ERROR in cell %d: %s', cl, conditionMessage(e))
        list(cell = cl, status = 'error', message = conditionMessage(e))
    })
}

## combined table across whatever cells this process ran
combined <- bind_rows(lapply(results, function(r) {
    if (is.null(r$latent) || is.null(r$latent$pooled)) return(NULL)
    r$latent$pooled %>% mutate(n_tips = r$n_tips, alpha_true = r$alpha_true,
                               prefix = r$prefix, .before = 1)
}))
if (!is.null(combined) && nrow(combined) > 0) {
    fn <- file.path(OUT_DIR, sprintf('260817_counts_loo_sweep_combined%s.csv',
                                     if (nzchar(cell_env)) paste0('_cell', cell_env) else ''))
    write.csv(combined, fn, row.names = FALSE)
    say('Combined latent-Spearman table -> %s', fn)
}

elapsed <- as.numeric(difftime(Sys.time(), t_start, units = 'hours'))
say('=== done in %.2f h | cells: %s ===', elapsed,
    paste(sprintf('%s:%s', names(results), vapply(results, function(r) r$status, character(1))),
          collapse = ' '))
