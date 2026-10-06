# 03_run_task.r -- Goal 1, branch-length baseline, step 03: fit SCOUT to one (alpha, layer, tree)
# job. Array task body; writes three CSVs per task.
#
# Usage:
#   Rscript 03_run_task.r 17        # task id given explicitly
#   Rscript 03_run_task.r           # or from SLURM_ARRAY_TASK_ID

.libPaths(strsplit(Sys.getenv('SCOUT_LIB',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/software/R/R-4.4.2/library'), ':')[[1]])
ROOT <- Sys.getenv('SCOUT_ROOT',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/1_tree_robustness')
source(file.path(ROOT, 'scripts', 'scout_perturb_lib.R'))
source(file.path(ROOT, 'scripts', '1_tree_perturb_bl', 'bl_paths.R'))

args    <- commandArgs(trailingOnly = TRUE)
task_id <- if (length(args) > 0) as.integer(args[1]) else
           as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID', NA))
if (is.na(task_id)) stop('No task id: pass one as an argument or set SLURM_ARRAY_TASK_ID.')

P        <- bl_paths(ROOT)
# Stamped into every output CSV so a merged table can always say which run a row came from. The
# ultrametric-restored variant sets it to 1_tree_perturb_bl_um; unset, this is the original run.
GOAL     <- Sys.getenv('SCOUT_GOAL', '1_tree_perturb_bl')
CORES    <- as.integer(Sys.getenv('SCOUT_CORES', Sys.getenv('SLURM_CPUS_PER_TASK', '16')))
JOBS     <- Sys.getenv('SCOUT_JOBS', P$jobs_sweep)
OUT_ROOT <- Sys.getenv('SCOUT_OUTROOT', P$outroot)

jobs <- read.csv(JOBS, stringsAsFactors = FALSE)
if (task_id < 1 || task_id > nrow(jobs)) stop(sprintf('task id %d outside 1..%d', task_id, nrow(jobs)))
job <- jobs[task_id, ]

# Optional smoke-test overrides: fit fewer genes and fewer regimes so a task finishes in minutes.
n_genes_cap <- as.integer(Sys.getenv('SCOUT_SMOKE_GENES', '0'))
regimes     <- strsplit(Sys.getenv('SCOUT_REGIMES', 'BM1,OU1,OUM'), ',')[[1]]

# The data layer being fitted: 'counts' is the observed matrix, 'evf' the latent OU process that
# produced it. They differ in exactly one thing that matters here -- EVFs are a latent trait that can
# go negative, so log-normalising them would give NaN; counts reach ~1e4 and must be normalised.
stopifnot(job$layer %in% c('counts', 'evf'))
normalize <- job$layer == 'counts'

results_dir <- task_dir_bl(job, OUT_ROOT)
stem <- task_stems_bl(job)

# Idempotence: a nohup rerun after a crash must resume, not redo days of finished work. Uses the
# same on-disk predicate as the status tracker (all output CSVs present AND non-empty), so a task
# killed mid-write is correctly seen as incomplete and gets redone.
if (as.logical(Sys.getenv('SCOUT_SKIP_DONE', 'FALSE')) && task_complete(results_dir, stem)) {
    cat(sprintf('[task %d] already complete -- skipping (delete its CSVs to force a redo)\n', task_id))
    quit(save = 'no', status = 0)
}

cat(sprintf('[task %d] layer=%s alpha=%s arm=%s intensity=%.2f rep=%d cores=%d normalize=%s\n',
            task_id, job$layer, alpha_str(job$alpha), job$arm, job$intensity, job$replicate,
            CORES, normalize))
cat(sprintf('  tree  : %s\n  data  : %s\n', job$tree_file, job$data_file))

phy <- read.tree(job$tree_file)
dat <- read.csv(job$data_file, row.names = 1)
rownames(dat) <- dat$species

if (n_genes_cap > 0) {
    gcols <- setdiff(colnames(dat), c('OUM', 'species'))
    keep  <- unlist(lapply(c('BM1', 'OU1', 'OUM'), function(m)
        head(grep(sprintf('^%s_', m), gcols, value = TRUE), n_genes_cap)))
    dat <- dat[, c(keep, 'OUM', 'species')]
    cat(sprintf('  SMOKE: capped to %d traits, regimes %s\n', ncol(dat) - 2, paste(regimes, collapse = ',')))
}

stopifnot(setequal(phy$tip.label, dat$species))

# H is recorded per task, not assumed. The nni arm preserves branch lengths but the collapse arm
# deletes internal edges, which changes root-to-tip path lengths -- on the published unit-branch-
# length tree that only ever changed an edge COUNT, and here it changes real time. Logging the
# realised height per tree is what makes that visible in the results rather than inferred.
H <- max(node.depth.edgelength(phy)[seq_len(length(phy$tip.label))])

dataset_id <- sprintf('%s|a%s|%s|%g|%03d', job$layer, alpha_str(job$alpha), job$arm,
                      job$intensity, job$replicate)
t0 <- Sys.time()
res <- score_classification(phy, dat, dataset_id = dataset_id, regimes = regimes,
                            normalize = normalize, cores = CORES, outdir = tempdir())
el <- as.numeric(difftime(Sys.time(), t0, units = 'mins'))

tag <- data.frame(goal = GOAL, layer = job$layer, alpha = job$alpha,
                  alphaH = job$alphaH, arm = job$arm, intensity = job$intensity,
                  replicate = job$replicate, tree_height = H, task_id = task_id,
                  elapsed_min = round(el, 2), stringsAsFactors = FALSE)
write_task_output(res, tag, results_dir, stem)

acc <- accuracy_summary(res$selected)
ov  <- acc[acc$metric_type == 'Overall', ]
cat(sprintf('[task %d] accuracy %.3f [%.3f, %.3f] on %d genes, H=%.3f, in %.1f min -> %s\n',
            task_id, ov$accuracy, ov$ci_lower, ov$ci_upper, ov$n, H, el, results_dir))
