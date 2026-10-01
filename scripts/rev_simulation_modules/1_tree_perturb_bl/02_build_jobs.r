# Goal 1 on a BRANCH-LENGTH baseline (1_tree_perturb_bl) step 02: build the job manifest.
#
# Same contract as build_jobs_v2.r: everything the runner needs is resolved and VALIDATED here, on
# the machine that holds the data, so a task only has to read a CSV and index into it. Every
# referenced path is proven to exist before the manifest is written -- a missing path found now
# costs nothing, found 40 hours into a nohup run it costs a dead task.
#
# TWO PHASES, because alpha is swept and fitting every simulated alpha would be wasteful.
#
#   SCOUT_PHASE=screen   Every simulated alpha, unperturbed baseline only, both layers.
#                        5 alphas x 2 layers = 10 tasks, about an hour. This is the gate: read the
#                        baseline accuracies before committing the long run. They should sit in the
#                        ~0.70-0.90 band and must not be at ceiling -- a saturated baseline has no
#                        headroom for a perturbation effect to show in -- and EVF must come out
#                        ABOVE counts at every alpha (EVF below counts would mean a normalisation
#                        bug, since EVFs are the latent process the counts are generated from).
#
#   SCOUT_PHASE=sweep    The alphas chosen from the screen (SCOUT_ALPHAS_KEEP) x every perturbation
#                        condition x replicates x both layers, PLUS the baseline row for each kept
#                        alpha.
#
# The baseline rows are deliberately repeated in the sweep manifest. They resolve to the same paths
# the screen already wrote, so SCOUT_SKIP_DONE skips them at zero cost -- and in exchange the sweep
# manifest is self-contained and 04_collate.r gets its per-alpha reference point without having to
# know that a separate phase ever existed.
#
# ALPHA IS NEVER RE-SIMULATED HERE. simulate_test_data shares one RNG stream across its alpha
# vector, so alphas can only be dropped, never added, without re-rolling the whole dataset. This
# script therefore reads the alphas that exist from baseline_manifest.csv and subsets them.

.libPaths(strsplit(Sys.getenv('SCOUT_LIB',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/software/R/R-4.4.2/library'), ':')[[1]])
ROOT <- Sys.getenv('SCOUT_ROOT',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/1_tree_robustness')
source(file.path(ROOT, 'scripts', '1_tree_perturb_bl', 'bl_paths.R'))

P       <- bl_paths(ROOT)
DATADIR <- Sys.getenv('SCOUT_DATADIR', P$datadir)
PHASE   <- Sys.getenv('SCOUT_PHASE', 'screen')
N_REP   <- as.integer(Sys.getenv('SCOUT_NREP', '15'))
LAYERS  <- strsplit(Sys.getenv('SCOUT_LAYERS', 'counts,evf'), ',')[[1]]
stopifnot(PHASE %in% c('screen', 'sweep'), all(LAYERS %in% c('counts', 'evf')))

base_man <- file.path(DATADIR, 'baseline_manifest.csv')
if (!file.exists(base_man))
    stop(sprintf('%s not found -- run 00_simulate_baseline.r first.', base_man))
bm <- read.csv(base_man, stringsAsFactors = FALSE)

# Which alphas get fitted. The screen fits them all; the sweep fits the subset the screen justified.
if (PHASE == 'screen') {
    keep <- bm$alpha
} else {
    keep <- as.numeric(strsplit(Sys.getenv('SCOUT_ALPHAS_KEEP', '0.5,1,3'), ',')[[1]])
    missing_a <- setdiff(alpha_str(keep), alpha_str(bm$alpha))
    if (length(missing_a))
        stop(sprintf('alpha %s was never simulated (have: %s). Alphas cannot be added without ',
                     paste(missing_a, collapse = ', '), paste(alpha_str(bm$alpha), collapse = ', ')),
             're-running 00_simulate_baseline.r, which re-rolls the data at EVERY alpha.')
}
bm <- bm[alpha_str(bm$alpha) %in% alpha_str(keep), , drop = FALSE]
stopifnot(nrow(bm) > 0)

# ---- the condition grid ------------------------------------------------------------------------
# The unperturbed baseline is carried as a row like any other (arm='baseline', intensity=0,
# replicate=0) so it flows through the manifest, the tracker, the skip-if-done check and the
# collation with no special-casing anywhere.
baseline_cond <- data.frame(arm = 'baseline', intensity = 0, replicate = 0,
                            tree_file = bm$tree_file[1], stringsAsFactors = FALSE)

if (PHASE == 'screen') {
    conds <- baseline_cond
} else {
    tman <- file.path(DATADIR, 'perturbed_tree_manifest_bl.csv')
    if (!file.exists(tman))
        stop(sprintf('%s not found -- run 01_make_trees.r first.', tman))
    pt <- read.csv(tman, stringsAsFactors = FALSE)
    pt <- pt[pt$replicate <= N_REP, c('arm', 'intensity', 'replicate', 'tree_file')]
    pt <- pt[order(pt$arm, pt$intensity, pt$replicate), ]
    if (!nrow(pt)) stop('the perturbed tree manifest has no replicates <= ', N_REP)
    conds <- rbind(baseline_cond, pt)
}

# ---- cross alpha x condition x layer -----------------------------------------------------------
jobs <- do.call(rbind, lapply(seq_len(nrow(bm)), function(i) {
    a <- bm[i, ]
    do.call(rbind, lapply(LAYERS, function(ly) data.frame(
        goal      = '1_tree_perturb_bl',
        layer     = ly,
        alpha     = a$alpha,
        alphaH    = a$alphaH,
        arm       = conds$arm,
        intensity = conds$intensity,
        replicate = conds$replicate,
        tree_file = conds$tree_file,
        data_file = if (ly == 'counts') a$counts_file else a$evf_file,
        stringsAsFactors = FALSE)))
}))
jobs <- jobs[order(jobs$layer, jobs$alpha, jobs$arm, jobs$intensity, jobs$replicate), ]
jobs <- cbind(task_id = seq_len(nrow(jobs)), jobs)

# ---- validate ----------------------------------------------------------------------------------
check <- function(label, paths) {
    missing <- unique(paths[!file.exists(paths)])
    cat(sprintf('%-30s %5d refs, %d missing %s\n', label, length(paths), length(missing),
                if (length(missing)) paste0('<- ', head(missing, 3), collapse = ' ') else ''))
    length(missing) == 0
}
ok <- c(check('trees',        jobs$tree_file),
        check('counts data',  jobs$data_file[jobs$layer == 'counts']),
        check('evf data',     jobs$data_file[jobs$layer == 'evf']))

# A duplicate condition key means two tasks resolve to the SAME output stem: one silently overwrites
# the other and the run reports more replicates than it actually has.
dup <- duplicated(jobs[, c('layer', 'alpha', 'arm', 'intensity', 'replicate')])
if (any(dup)) {
    print(head(jobs[dup, c('layer', 'alpha', 'arm', 'intensity', 'replicate')], 5), row.names = FALSE)
    stop(sprintf('%d duplicate condition key(s) -- these tasks would overwrite each other.', sum(dup)))
}
stopifnot(all(ok), identical(jobs$task_id, seq_len(nrow(jobs))))

# Prove the output stems are unique too, not just the keys -- that is the property that actually
# protects the files, and it goes through the same resolvers the tasks and the tracker use.
stems <- file.path(task_dir_bl(jobs, P$outroot), task_stems_bl(jobs))
stopifnot(!any(duplicated(stems)))

f <- if (PHASE == 'screen') P$jobs_screen else P$jobs_sweep
write.csv(jobs, f, row.names = FALSE)

cat(sprintf('\nphase %s: %d tasks | alphas %s | layers %s\n', PHASE, nrow(jobs),
            paste(alpha_str(sort(unique(jobs$alpha))), collapse = ','), paste(LAYERS, collapse = ',')))
if (PHASE == 'sweep')
    cat(sprintf('conditions: %d (incl. baseline) x %d replicates\n',
                length(unique(paste(jobs$arm, jobs$intensity))), N_REP))
print(as.data.frame(table(layer = jobs$layer, arm = jobs$arm, intensity = jobs$intensity)) |>
      subset(Freq > 0), row.names = FALSE)
cat(sprintf('\nmanifest: %s\n', f))
cat(sprintf('output  : %s\n', P$outroot))
