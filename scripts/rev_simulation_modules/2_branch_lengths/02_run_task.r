# Goal 2 (2_branch_lengths) step 02: fit ONE replicate's matrix on BOTH trees. This is the array
# task body.
#
# Usage:
#   Rscript 02_run_task.r            # task id from SLURM_ARRAY_TASK_ID
#   Rscript 02_run_task.r 7          # task id given explicitly (for smoke tests)
#
# One task = one (replicate, matrix_type) pair, fitted on the true-branch-length tree and on the
# unit-branch-length tree. Pairing them in a single task guarantees both arms see byte-identical data
# and keeps each task under the wall-clock limit; splitting them would double the I/O for no gain.
#
# EVFs are the latent OU process and contain negative values, so they must NOT be log-normalised --
# log1p would give NaN. Counts reach ~1e5, far above the theta upper bound, so they must be.

# SCOUT_LIB may be a colon-separated path LIST, so a private build (e.g. the support_clip
# SCOUT) can be prepended while its dependencies still resolve from the shared library.
.libPaths(strsplit(Sys.getenv('SCOUT_LIB',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/software/R/R-4.4.2/library'), ':')[[1]])
ROOT <- Sys.getenv('SCOUT_ROOT',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/1_tree_robustness')
source(file.path(ROOT, 'scripts', 'scout_perturb_lib.R'))

args    <- commandArgs(trailingOnly = TRUE)
task_id <- if (length(args) > 0) as.integer(args[1]) else
           as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID', NA))
if (is.na(task_id)) stop('No task id: pass one as an argument or set SLURM_ARRAY_TASK_ID.')

CORES <- as.integer(Sys.getenv('SCOUT_CORES', Sys.getenv('SLURM_CPUS_PER_TASK', '16')))
JOBS  <- Sys.getenv('SCOUT_JOBS', file.path(ROOT, 'scripts', 'jobs_2_branch_lengths.csv'))
OUT_ROOT <- Sys.getenv('SCOUT_OUTROOT', file.path(ROOT, 'output', '2_branch_lengths'))
jobs  <- read.csv(JOBS, stringsAsFactors = FALSE)
if (task_id < 1 || task_id > nrow(jobs)) stop(sprintf('task id %d outside 1..%d', task_id, nrow(jobs)))
job <- jobs[task_id, ]

n_genes_cap <- as.integer(Sys.getenv('SCOUT_SMOKE_GENES', '0'))
regimes     <- strsplit(Sys.getenv('SCOUT_REGIMES', 'BM1,OU1,OUM'), ',')[[1]]

normalize <- job$matrix_type == 'counts'   # EVFs must not be log-normalised
results_dir <- task_dir(job, '2_branch_lengths', OUT_ROOT)

# Idempotence: a nohup rerun after a crash must resume, not redo days of finished work. Uses the
# same on-disk predicate as the status tracker (all output CSVs present AND non-empty), so a task
# killed mid-write is correctly seen as incomplete and gets redone. Off by default so the sbatch
# path is unchanged; run_local.sh turns it on.
if (as.logical(Sys.getenv('SCOUT_SKIP_DONE', 'FALSE')) &&
    task_complete(results_dir, task_stems(job, '2_branch_lengths'))) {
    cat(sprintf('[task %d] already complete -- skipping (delete its CSVs to force a redo)\n', task_id))
    quit(save = 'no', status = 0)
}

cat(sprintf('[task %d] matrix=%s rep=%d cores=%d normalize=%s\n',
            task_id, job$matrix_type, job$replicate, CORES, normalize))

dat <- read.csv(job$data_file, row.names = 1)
rownames(dat) <- dat$species

if (n_genes_cap > 0) {
    gcols <- setdiff(colnames(dat), c('OUM', 'species'))
    keep  <- unlist(lapply(c('BM1', 'OU1', 'OUM'), function(m)
        head(grep(sprintf('^%s_', m), gcols, value = TRUE), n_genes_cap)))
    dat <- dat[, c(keep, 'OUM', 'species')]
    cat(sprintf('  SMOKE: capped to %d traits, regimes %s\n', ncol(dat) - 2, paste(regimes, collapse = ',')))
}

arms <- list(true_bl = job$tree_bl, unit_bl = job$tree_unit)
out  <- list()

for (nm in names(arms)) {
    phy <- read.tree(arms[[nm]])
    stopifnot(setequal(phy$tip.label, dat$species))
    H <- max(node.depth.edgelength(phy)[seq_len(length(phy$tip.label))])

    dataset_id <- sprintf('%s|%s|%03d', job$matrix_type, nm, job$replicate)
    t0 <- Sys.time()
    res <- score_classification(phy, dat, dataset_id = dataset_id, regimes = regimes,
                                normalize = normalize, cores = CORES, outdir = tempdir())
    el <- as.numeric(difftime(Sys.time(), t0, units = 'mins'))

    tag <- data.frame(goal = '2_branch_lengths', matrix_type = job$matrix_type, arm = nm,
                      replicate = job$replicate, alpha = job$alpha, height = H,
                      task_id = task_id, elapsed_min = round(el, 2), stringsAsFactors = FALSE)
    write_task_output(res, tag, results_dir,
                      sprintf('%s_%s_rep%03d', job$matrix_type, nm, job$replicate))

    ov <- accuracy_summary(res$selected)
    ov <- ov[ov$metric_type == 'Overall', ]
    out[[nm]] <- ov$accuracy
    cat(sprintf('  %-7s H=%6.2f  accuracy %.3f [%.3f, %.3f]  %.1f min\n',
                nm, H, ov$accuracy, ov$ci_lower, ov$ci_upper, el))
}

# The paired difference is the estimate this goal exists to produce.
cat(sprintf('[task %d] paired difference (true_bl - unit_bl) = %+.3f\n',
            task_id, out$true_bl - out$unit_bl))
