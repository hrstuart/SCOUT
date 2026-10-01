# Goal 1 on a BRANCH-LENGTH baseline (1_tree_perturb_bl) step 00:
# build the baseline tree and simulate the counts/EVF matrices every perturbation arm is fitted
# against, at each alpha in the grid.
#
# WHAT IS DIFFERENT FROM 1_tree_perturb/00_simulate_baseline.r, AND WHY IT IS THE WHOLE POINT.
#
# That script stamps `simtree$edge.length <- rep(1, Nedge(simtree))` BEFORE simulating. This one
# does not. That single line is the entire change; everything downstream follows from it.
#
# The consequence is not cosmetic. generate_tree_hbd_reverse returns an ULTRAMETRIC tree -- every
# tip at depth 4.7388. Overwriting the edge lengths with 1 turns it into a NON-ultrametric tree of
# height 18 with tip depths spanning 2-18, so BM tip variance varies ~9-fold across tips
# (cv_var_bm = 0.285) while a stationary OU at alphaH = 54 is flat (cv_var_ou = 0). Tip-depth
# heteroscedasticity therefore becomes BM's dominant marginal signature -- a per-tip variance
# gradient that model selection can lean on and that NO amount of topological damage removes. That
# is the most plausible explanation for why the published perturbation sweep came out flat: NNI and
# collapse held 0.822 (counts) / 0.947 (EVF) from 1% intensity all the way to 75%.
#
# Keeping the true branch lengths sets cv_var_bm = 0 at every alpha (checked analytically, see the
# identifiability reference this script writes). BM/OU separation then has to come from the
# covariance structure -- which is exactly what NNI, collapse and shuffle damage. This run re-asks
# the perturbation question in the regime where the topology should actually matter.
#
# Two smaller consequences, both handled below:
#   t0    The tree is now height 4.7388, not 18, so BM tip sd is 2.177 instead of 4.243 and the
#         5-sd margin gives t0 = 11 rather than 22. Derived from H the way
#         2_branch_lengths/01_simulate_replicates.r derives it, never hard-coded.
#   alpha alphaH is now 4.7388 * alpha instead of 18 * alpha, so the SAME alpha lands in a very
#         different identifiability regime. alpha 3 here has cor_shape 0.269; the published run's
#         alpha 3 had 0.190. That is why alpha is swept rather than fixed.
#
# THE ALPHA GRID IS LOCKED. simulate_test_data shares ONE RNG stream across the alpha vector, so
# changing the length or the order of this grid changes the simulated data at every alpha after the
# first -- including the EVFs. Adding an alpha later silently re-rolls the run. The grid is
# therefore simulated in full ONCE and alphas are only ever DROPPED downstream (02_build_jobs.r
# picks which ones are fitted). Do not edit SCOUT_ALPHAS after the data exists.
#
# Signal model is the affine framework the published v2 run used -- param_mode = 'affine',
# scale_s = 2, support_clip = TRUE, 30 genes / 20 EVFs -- so the only moving part between this run
# and that one is the branch lengths.

# SCOUT_LIB may be a colon-separated path LIST, so a private build (e.g. the support_clip
# SCOUT) can be prepended while its dependencies still resolve from the shared library.
.libPaths(strsplit(Sys.getenv('SCOUT_LIB',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/software/R/R-4.4.2/library'), ':')[[1]])
ROOT <- Sys.getenv('SCOUT_ROOT',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/1_tree_robustness')
source(file.path(ROOT, 'scripts', 'scout_perturb_lib.R'))
source(file.path(ROOT, 'scripts', '1_tree_perturb_bl', 'bl_paths.R'))
suppressPackageStartupMessages({library(castor); library(expm)})

P <- bl_paths(ROOT)
DATADIR <- Sys.getenv('SCOUT_DATADIR', P$datadir)

# This run must not be able to disturb the published data. The target directory has to be absent or
# empty unless the overwrite is deliberate.
guard_empty(DATADIR, 'baseline data directory')
dir.create(DATADIR, recursive = TRUE, showWarnings = FALSE)

PREFIX     <- Sys.getenv('SCOUT_BL_PREFIX', '260826_simtree_256cells_truebl')
NCELLS     <- 256
NGENES     <- as.integer(Sys.getenv('SCOUT_NGENES', '30'))   # per model -> 3x this many traits
NEVFS      <- as.integer(Sys.getenv('SCOUT_NEVFS',  '20'))   # per model -> 3x this many EVF traits
SIGMA      <- 1
THETA_STEP <- 2
ALPHAS     <- as.numeric(strsplit(Sys.getenv('SCOUT_ALPHAS', '0.25,0.5,1,3,6'), ',')[[1]])
SEED       <- 1

PARAM_MODE   <- Sys.getenv('SCOUT_PARAM_MODE', 'affine')
SIGNAL_GAIN  <- as.numeric(Sys.getenv('SCOUT_SIGNAL_GAIN', '0.3'))
BURST_SCALE  <- as.numeric(Sys.getenv('SCOUT_BURST_SCALE', '1'))
SCALE_S      <- as.numeric(Sys.getenv('SCOUT_SCALE_S', '2'))
SUPPORT_CLIP <- as.logical(Sys.getenv('SCOUT_SUPPORT_CLIP', 'TRUE'))   # inert unless param_mode=affine
gain_vec     <- if (PARAM_MODE == 'affine') c(NA, NA, SIGNAL_GAIN) else c(NA, NA, NA)

if (!'param_mode' %in% names(formals(simulate_test_data)))
    stop('The installed SCOUT has no param_mode argument -- point SCOUT_LIB at the sim-signal-gain build.')
if (PARAM_MODE == 'affine' && !'support_clip' %in% names(formals(simulate_test_data)))
    stop('param_mode=affine needs the support-clip build -- point SCOUT_LIB at it.')

set.seed(SEED)

# --- tree -------------------------------------------------------------------------------------
# Identical generator, seed and transition matrix to 1_tree_perturb/00_simulate_baseline.r, so the
# topology and the tip states are the SAME draw the published run used. That is deliberate: it makes
# this run a matched control rather than an independent one, and it is what allows the perturbed
# trees built in 01_make_trees.r to be edge-matched to the published set.
P_mk <- matrix(c(0.90, 0.08, 0.02,
                 0.10, 0.80, 0.10,
                 0.02, 0.08, 0.90), 3, 3, byrow = TRUE)
Q <- expm::logm(P_mk)

simtree <- castor::generate_tree_hbd_reverse(NCELLS, rho = 1, lambda = 1, mu = 0)$trees[[1]]
simtree$tip.label <- paste0('t', simtree$tip.label)
st <- castor::simulate_mk_model(simtree, Q = Q, include_nodes = FALSE)
simtree$states <- as.factor(paste0('state', st$tip_states))

# >>> NO `simtree$edge.length <- rep(1, Nedge(simtree))` HERE. That omission IS this run. <<<

depths <- node.depth.edgelength(simtree)[seq_len(NCELLS)]
H <- max(depths)

# Margin, in units of the BM tip sd, between the root optimum and zero. The EVF matrix is fitted
# directly as its own arm, and BM EVFs are the raw latent process, so a non-positive EVF is an
# unfittable gene (theta is bounded below at 1e-10, SCOUT_EM.R:213). 5 sd is what the published v2
# run used. Derived from THIS tree's height, never carried over from the unit-branch-length tree.
T0_SD_MARGIN <- as.numeric(Sys.getenv('SCOUT_T0_MARGIN', '5'))
T0 <- ceiling(T0_SD_MARGIN * sqrt(SIGMA) * sqrt(H))

cat(sprintf('tips              : %d\n', NCELLS))
cat(sprintf('height            : %.4f   ultrametric = %s\n', H, is.ultrametric(simtree)))
cat(sprintf('tip depths        : %.4f - %.4f (median %.4f)\n', min(depths), max(depths), median(depths)))
cat(sprintf('edge lengths      : min %.3e  median %.4f  max %.4f\n',
            min(simtree$edge.length), median(simtree$edge.length), max(simtree$edge.length)))
cat(sprintf('t0                : %d   (%.1f BM tip sd; BM tip sd here is %.3f)\n', T0, T0_SD_MARGIN, sqrt(H)))
cat(sprintf('OUM thetas        : %s\n', paste(seq(T0, by = THETA_STEP, length.out = 3), collapse = ', ')))
cat(sprintf('alpha simulated   : %s\n', paste(ALPHAS, collapse = ', ')))
cat(sprintf('alpha*H simulated : %s\n', paste(sprintf('%.2f', ALPHAS * H), collapse = ', ')))
cat('tip states        : '); print(table(simtree$states))

# --- baseline sanity --------------------------------------------------------------------------
# These are the properties the whole design rests on. If any of them moved, the run is not the run
# this was planned as -- fail now, not 630 tasks later.
stopifnot(is.ultrametric(simtree))
if (abs(H - 4.738808) > 1e-5)
    stop(sprintf('baseline height %.6f != expected 4.738808 -- the tree generator moved.', H))
tip_tab <- table(as.character(simtree$states))
if (!identical(as.integer(tip_tab[c('state1', 'state2', 'state3')]), c(106L, 61L, 89L)))
    stop(sprintf('tip state counts are %s, expected 106/61/89 -- the Mk simulation moved.',
                 paste(tip_tab, collapse = '/')))
# di2multi(tol = 1e-9) is how perturb_collapse turns chosen edges into polytomies. If any edge were
# already at or below that tolerance it would be collapsed too, so the collapse arm would remove
# more than its k edges and stop being edge-matched to the published run.
if (min(simtree$edge.length) <= 1e-9)
    stop('an edge is <= di2multi tolerance (1e-9); the collapse arm would over-collapse.')

simtree <- infer_anc(simtree)   # node labels required by simulate_test_data(sim_lin = FALSE)
stopifnot(!is.null(simtree$node.label), length(simtree$node.label) == simtree$Nnode)

treefile <- file.path(DATADIR, paste0(PREFIX, '.nwk'))
tree_out <- simtree; tree_out$node.label <- NULL   # SCOUT re-infers node states per regime
ape::write.tree(tree_out, treefile)

# The tree is regenerated from SEED = 1 on every run, and the perturbed trees are only valid against
# THIS tree. If a library or RNG change shifts the generator, the perturbed-tree manifest silently
# starts pointing at trees derived from a different baseline -- fail here instead. Recorded
# 2026-08-26 from the seed-1 tree with its true branch lengths kept.
BASELINE_TREE_MD5 <- '6bf1b5f56cc5d96295661f46825ab22d'
tree_md5 <- unname(tools::md5sum(treefile))
cat(sprintf('baseline tree md5 : %s (%s)\n', tree_md5,
            if (identical(tree_md5, BASELINE_TREE_MD5)) 'matches the pinned true-BL baseline'
            else 'MISMATCH'))
if (!identical(tree_md5, BASELINE_TREE_MD5))
    stop(sprintf('baseline tree md5 %s != expected %s -- the generator moved and nothing downstream is comparable.',
                 tree_md5, BASELINE_TREE_MD5))
saveRDS(simtree, file.path(DATADIR, paste0(PREFIX, '_with_states.rds')))

# --- analytic identifiability reference (no fitting, seconds) ---------------------------------
# This is where cv_var_bm = 0 is recorded for the record: the quantity whose removal is the reason
# for the run. It also gives cor_shape at each simulated alpha, which is how 02_build_jobs.r's alpha
# subset is justified.
ref <- identifiability_reference(simtree, alphaH = sort(unique(c(ALPHAS * H, 1, 2, 5, 10, 20, 30, 54))),
                                 sigma = SIGMA, add.root = FALSE)
write.csv(ref, file.path(DATADIR, paste0(PREFIX, '_identifiability_reference.csv')), row.names = FALSE)
cat('\n-- analytic BM/OU separability on this tree (cv_var_bm = 0 is the point) --\n')
print(ref[, c('alphaH', 'alpha', 'cor_shape', 'resid_sd', 'cv_var_bm', 'cv_var_ou')], row.names = FALSE)

# --- simulate ---------------------------------------------------------------------------------
cat('\n-- simulating --\n')
cat(sprintf('signal model: param_mode=%s signal_gain=%s burst_scale=%g scale_s=%g support_clip=%s\n',
            PARAM_MODE, if (PARAM_MODE == 'affine') format(SIGNAL_GAIN) else 'n/a',
            BURST_SCALE, SCALE_S,
            if (PARAM_MODE == 'affine') SUPPORT_CLIP else 'n/a (rank)'))

sim_args <- list(ngenes = NGENES, ncells = NCELLS, tree = simtree,
                 outdir = DATADIR, out_prefix = PREFIX,
                 a = ALPHAS, s = SIGMA, t0 = T0, theta_step = THETA_STEP,
                 nevfs = NEVFS, sim_lin = FALSE,
                 param_mode = PARAM_MODE, signal_gain = gain_vec,
                 burst_scale = BURST_SCALE, scale_s = SCALE_S)
if ('support_clip' %in% names(formals(simulate_test_data)))
    sim_args$support_clip <- SUPPORT_CLIP

sim <- do.call(simulate_test_data, sim_args)

# --- sanity: traits must stay positive --------------------------------------------------------
cat('\n-- trait sanity (theta cannot be fitted below 1e-10) --\n')
EVF_STRICT <- as.logical(Sys.getenv('SCOUT_EVF_STRICT', 'TRUE'))
bad <- 0L
bad_evf <- character(0)
for (k in names(sim$simulated_OU$counts)) {
    cn <- sim$simulated_OU$counts[[k]]; ev <- sim$simulated_OU$evfs[[k]]
    gc_ <- setdiff(colnames(cn), c('OUM', 'species'))
    ge_ <- setdiff(colnames(ev), c('OUM', 'species'))
    min_bm_evf <- min(ev[, grep('^BM1_', ge_, value = TRUE)])
    cat(sprintf('%s : counts min %.1f max %.0f | EVF min %.3f (BM class min %.3f)\n',
                k, min(cn[, gc_]), max(cn[, gc_]), min(ev[, ge_]), min_bm_evf))
    if (min(cn[, gc_]) < 0) { cat('   FAIL negative counts\n'); bad <- bad + 1L }
    if (min_bm_evf <= 0) {
        cat(sprintf('   %s BM EVFs reach <= 0 (min %.3f, %.2f BM tip sd below the root optimum)\n',
                    if (EVF_STRICT) 'FAIL' else 'WARN', min_bm_evf, (T0 - min_bm_evf) / sqrt(SIGMA * H)))
        if (EVF_STRICT) bad_evf <- c(bad_evf, k)
    }
}
stopifnot(bad == 0L)

# The EVF matrix is fitted directly as its own arm, so a non-positive EVF is not a footnote -- it is
# an unfittable gene. Fail here rather than 630 tasks later. Costs one re-simulation; the
# alternative costs the arm.
if (EVF_STRICT && length(bad_evf))
    stop(sprintf('BM EVFs reach <= 0 at alpha %s. Raise SCOUT_T0_MARGIN above %g and re-simulate.',
                 paste(bad_evf, collapse = ', '), T0_SD_MARGIN))
if (EVF_STRICT)
    cat('EVFs are strictly positive at every alpha -- t0 is adequate for the EVF arm.\n')

# --- count realism ------------------------------------------------------------------------------
cat('\n-- count realism (targets: per-gene-cell max ~1e4, per-cell total < 1e5) --\n')
cat(sprintf('%-10s %8s %8s %9s %10s %8s\n', 'alpha', 'median', 'p99', 'max', 'lib med', 'zeros'))
for (k in names(sim$simulated_OU$counts)) {
    cn <- sim$simulated_OU$counts[[k]]
    M  <- as.matrix(cn[, setdiff(colnames(cn), c('OUM', 'species'))])
    cat(sprintf('%-10s %8.0f %8.0f %9.0f %10.0f %8.3f\n', k, median(M), quantile(M, 0.99),
                max(M), median(rowSums(M)), mean(M == 0)))
}
cat(sprintf('support ceiling on a single gene-cell count: %.0f\n', 10^4.043 * SCALE_S))

# --- what downstream reads ----------------------------------------------------------------------
# Written so 02_build_jobs.r resolves data files from a manifest rather than re-deriving filename
# conventions, and so the exact (alpha, t0, alphaH) triple each matrix was simulated under is on
# disk and auditable.
man <- data.frame(
    alpha       = ALPHAS,
    alphaH      = ALPHAS * H,
    height      = H,
    t0          = T0,
    sigma       = SIGMA,
    theta_step  = THETA_STEP,
    tree_file   = treefile,
    counts_file = file.path(DATADIR, sprintf('%s_alpha_%s_sigma_%s_OUEVF_counts.csv',     PREFIX, ALPHAS, SIGMA)),
    evf_file    = file.path(DATADIR, sprintf('%s_alpha_%s_sigma_%s_OUEVF_evf_counts.csv', PREFIX, ALPHAS, SIGMA)),
    stringsAsFactors = FALSE)
stopifnot(all(file.exists(man$counts_file)), all(file.exists(man$evf_file)))
manfile <- file.path(DATADIR, 'baseline_manifest.csv')
write.csv(man, manfile, row.names = FALSE)

cat(sprintf('\nBaseline written to %s\n', DATADIR))
cat(sprintf('tree     : %s\n', treefile))
cat(sprintf('manifest : %s\n', manfile))
print(man[, c('alpha', 'alphaH', 't0')], row.names = FALSE)
