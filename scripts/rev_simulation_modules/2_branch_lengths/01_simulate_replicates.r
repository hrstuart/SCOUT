# 01_simulate_replicates.r -- Goal 2, step 01: simulate independent replicate datasets (fresh
# tree, tip states and data each) on true-branch-length trees, for paired fitting on true vs unit
# branch lengths. t0 is set per replicate from tree height.
#
# Usage:
#   Rscript 01_simulate_replicates.r        # optional env: SCOUT_ALPHA_G2 (default 3), SCOUT_NREP, SCOUT_DATADIR

# SCOUT_LIB may be a colon-separated path LIST, so a private build (e.g. the support_clip
# SCOUT) can be prepended while its dependencies still resolve from the shared library.
.libPaths(strsplit(Sys.getenv('SCOUT_LIB',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/software/R/R-4.4.2/library'), ':')[[1]])
ROOT <- Sys.getenv('SCOUT_ROOT',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/1_tree_robustness')
source(file.path(ROOT, 'scripts', 'scout_perturb_lib.R'))
suppressPackageStartupMessages({library(castor); library(expm); library(parallel)})

# DATADIR is overridable so a second signal model can be simulated side by side without touching the
# published rank data. NGENES/NEVFS default to the original values, so a bare run reproduces it.
DATADIR <- Sys.getenv('SCOUT_DATADIR', file.path(ROOT, 'data', '2_branch_lengths'))
dir.create(DATADIR, recursive = TRUE, showWarnings = FALSE)

NCELLS     <- 256
NGENES     <- as.integer(Sys.getenv('SCOUT_NGENES', '50'))   # per model -> 3x this many counts traits
NEVFS      <- as.integer(Sys.getenv('SCOUT_NEVFS',  '20'))   # per model -> 3x this many EVF traits
SIGMA      <- 1
THETA_STEP <- 2
ALPHA      <- as.numeric(Sys.getenv('SCOUT_ALPHA_G2', '3'))
N_REP      <- as.integer(Sys.getenv('SCOUT_NREP', '50'))
BASE_SEED  <- 50000
SIMCORES   <- as.integer(Sys.getenv('SCOUT_SIMCORES', '25'))

# Margin, in units of the BM tip sd, between the root optimum and zero. See the header note.
T0_SD_MARGIN <- as.numeric(Sys.getenv('SCOUT_T0_MARGIN', '4.3'))

# Signal model defaults: the ORIGINAL framework, param_mode = 'rank', burst_scale = 1, scale_s = 10.
#
# KNOWN REALISM COST. These are SymSim TRUE transcript counts, not observed UMIs -- generate_EvoCounts
# is rbeta then rpois and stops, and True2ObservedCounts is only wired into the TedSim path
# (simulate.R:50). So the counts run large (max ~95,611 for one gene in one cell, per-cell totals over
# 100,000) and the zero fraction is implausibly low at ~0.076. Accepted here as the price of the
# original framework; scale_s = 2 would be capture efficiency 0.2, SymSim's own UMI default.
#
# scale_s does NOT trade against fit: at scale_s 10 Poisson log-noise is 0.037 against bursting noise
# of 0.566, under 0.5% of the noise variance. It is a realism knob only. burst_scale is the one knob
# that moves the fit here (+0.073 at burst 10), and it is left at SymSim's calibrated 1.
#
# The alternatives need no code change:
#   SCOUT_PARAM_MODE=affine SCOUT_SCALE_S=2 SCOUT_SUPPORT_CLIP=TRUE   realistic affine model
#   SCOUT_BURST_SCALE=10                                             the +0.073 fit knob
PARAM_MODE  <- Sys.getenv('SCOUT_PARAM_MODE', 'rank')
SIGNAL_GAIN <- as.numeric(Sys.getenv('SCOUT_SIGNAL_GAIN', '0.3'))
BURST_SCALE <- as.numeric(Sys.getenv('SCOUT_BURST_SCALE', '1'))
SCALE_S     <- as.numeric(Sys.getenv('SCOUT_SCALE_S', '10'))
SUPPORT_CLIP <- as.logical(Sys.getenv('SCOUT_SUPPORT_CLIP', 'TRUE'))   # inert unless param_mode=affine
gain_vec    <- if (PARAM_MODE == 'affine') c(NA, NA, SIGNAL_GAIN) else c(NA, NA, NA)

if (!'param_mode' %in% names(formals(simulate_test_data)))
    stop('The installed SCOUT has no param_mode argument -- point SCOUT_LIB at the sim-signal-gain build.')
if (PARAM_MODE == 'affine' && !'support_clip' %in% names(formals(simulate_test_data)))
    stop('param_mode=affine needs the support-clip build -- point SCOUT_LIB at it.')
PASS_CLIP <- 'support_clip' %in% names(formals(simulate_test_data))

P <- matrix(c(0.90, 0.08, 0.02,
              0.10, 0.80, 0.10,
              0.02, 0.08, 0.90), 3, 3, byrow = TRUE)
Q <- expm::logm(P)

cat(sprintf('%d replicates | %d genes + %d EVFs per model | alpha %g | sigma %g\n',
            N_REP, NGENES, NEVFS, ALPHA, SIGMA))
cat(sprintf('t0 per replicate = ceiling(%g * sqrt(sigma) * sqrt(H))\n', T0_SD_MARGIN))
cat(sprintf('signal model: param_mode=%s signal_gain=%s burst_scale=%g scale_s=%g support_clip=%s\n',
            PARAM_MODE, if (PARAM_MODE == 'affine') format(SIGNAL_GAIN) else 'n/a',
            BURST_SCALE, SCALE_S,
            if (PARAM_MODE == 'affine') SUPPORT_CLIP else 'n/a (rank)'))

sim_one <- function(r) {
    seed <- BASE_SEED + r
    set.seed(seed)
    prefix <- sprintf('rep%03d', r)

    phy <- castor::generate_tree_hbd_reverse(NCELLS, rho = 1, lambda = 1, mu = 0)$trees[[1]]
    phy$tip.label <- paste0('t', phy$tip.label)
    st <- castor::simulate_mk_model(phy, Q = Q, include_nodes = FALSE)
    phy$states <- as.factor(paste0('state', st$tip_states))
    phy <- infer_anc(phy)

    H <- max(node.depth.edgelength(phy)[seq_len(NCELLS)])

    # Per-replicate root optimum: T0_SD_MARGIN BM tip sds above zero, using THIS tree's height.
    t0_r <- ceiling(T0_SD_MARGIN * sqrt(SIGMA) * sqrt(H))

    sim_args <- list(ngenes = NGENES, ncells = NCELLS, tree = phy,
                     outdir = DATADIR, out_prefix = prefix,
                     a = ALPHA, s = SIGMA, t0 = t0_r, theta_step = THETA_STEP,
                     nevfs = NEVFS, sim_lin = FALSE,
                     param_mode = PARAM_MODE, signal_gain = gain_vec,
                     burst_scale = BURST_SCALE, scale_s = SCALE_S)
    if (PASS_CLIP) sim_args$support_clip <- SUPPORT_CLIP
    sim <- do.call(simulate_test_data, sim_args)

    # True-branch-length tree, and its unit-branch-length counterpart. Both are written so the array
    # tasks only read files, and so the pair actually fitted is auditable.
    tree_bl <- phy; tree_bl$node.label <- NULL
    f_bl <- file.path(DATADIR, sprintf('%s_truebl.nwk', prefix))
    ape::write.tree(tree_bl, f_bl)

    tree_unit <- tree_variant(tree_bl, branch = 'unit', depths = 'keep')
    f_unit <- file.path(DATADIR, sprintf('%s_unitbl.nwk', prefix))
    ape::write.tree(tree_unit, f_unit)

    key <- names(sim$simulated_OU$counts)[1]
    cn  <- sim$simulated_OU$counts[[key]]
    ev  <- sim$simulated_OU$evfs[[key]]
    gc_ <- setdiff(colnames(cn), c('OUM', 'species'))
    ge_ <- setdiff(colnames(ev), c('OUM', 'species'))

    data.frame(
        replicate   = r,
        seed        = seed,
        tree_bl     = f_bl,
        tree_unit   = f_unit,
        counts_file = file.path(DATADIR, sprintf('%s_alpha_%s_sigma_%s_OUEVF_counts.csv', prefix, ALPHA, SIGMA)),
        evf_file    = file.path(DATADIR, sprintf('%s_alpha_%s_sigma_%s_OUEVF_evf_counts.csv', prefix, ALPHA, SIGMA)),
        height_bl   = H,
        height_unit = max(node.depth.edgelength(tree_unit)[seq_len(NCELLS)]),
        t0          = t0_r,
        evf_margin_sd = min(ev[, ge_]) / sqrt(SIGMA * H),   # how close the worst draw came to zero
        n_states    = length(unique(phy$states)),
        min_counts  = min(cn[, gc_]),
        min_evf     = min(ev[, ge_]),
        min_evf_bm  = min(ev[, grep('^BM1_', ge_, value = TRUE)]),
        mean_count  = mean(as.matrix(cn[, gc_])),
        frac_zero   = mean(as.matrix(cn[, gc_]) == 0),
        stringsAsFactors = FALSE)
}

man <- do.call(rbind, mclapply(seq_len(N_REP), sim_one, mc.cores = SIMCORES))
man <- man[order(man$replicate), ]

manfile <- file.path(DATADIR, 'replicate_manifest.csv')
write.csv(man, manfile, row.names = FALSE)

cat('\n-- sanity --\n')
cat(sprintf('tree heights (true BL) : %.3f - %.3f\n', min(man$height_bl), max(man$height_bl)))
cat(sprintf('tree heights (unit BL) : %g - %g   (alpha*H differs by ~%.1fx between arms)\n',
            min(man$height_unit), max(man$height_unit),
            mean(man$height_unit) / mean(man$height_bl)))
cat(sprintf('t0 across replicates   : %d - %d\n', min(man$t0), max(man$t0)))
cat(sprintf('states per replicate   : %s\n', paste(sort(unique(man$n_states)), collapse = ',')))
cat(sprintf('min counts across reps : %.1f\n', min(man$min_counts)))
cat(sprintf('min EVF across reps    : %.3f  (BM class %.3f)\n', min(man$min_evf), min(man$min_evf_bm)))
cat(sprintf('worst EVF margin       : %.2f BM tip sd above zero (want > 0)\n', min(man$evf_margin_sd)))
cat(sprintf('mean count / frac zero : %.0f / %.3f\n', mean(man$mean_count), mean(man$frac_zero)))

stopifnot(all(file.exists(man$tree_bl)), all(file.exists(man$tree_unit)),
          all(file.exists(man$counts_file)), all(file.exists(man$evf_file)),
          all(man$min_counts >= 0), nrow(man) == N_REP)

# This is the check the previous run failed. theta is bounded below at 1e-10 (SCOUT_EM.R:213), so a
# non-positive EVF makes that replicate uninterpretable in the EVF arm. Fail here, not 100 array
# tasks later.
if (min(man$min_evf) <= 0) {
    bad <- man$replicate[man$min_evf <= 0]
    print(man[man$min_evf <= 0, c('replicate', 'height_bl', 't0', 'min_evf', 'min_evf_bm')])
    stop(sprintf('EVFs reach <= 0 in %d replicate(s): %s. Raise SCOUT_T0_MARGIN above %g.',
                 length(bad), paste(bad, collapse = ','), T0_SD_MARGIN))
}
cat('\nEVFs are strictly positive in every replicate -- t0 is adequate for the EVF arm.\n')
if (any(man$n_states < 3)) cat('WARNING: some replicates have fewer than 3 tip states; OUM is degenerate there.\n')

cat(sprintf('\n%d replicates written to %s\nmanifest: %s\n', nrow(man), DATADIR, manfile))
