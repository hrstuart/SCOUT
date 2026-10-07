# 01b_build_jobs_alpha.r -- Goal 2, step 01b: build the job manifest for ONE alpha without touching
# existing manifests. Checks that the data directory was simulated at the requested alpha.
#
# Usage:
#   SCOUT_ALPHA_G2=0.5 SCOUT_DATA_G2=<data dir for that alpha> Rscript 01b_build_jobs_alpha.r
#   # optional: SCOUT_JOBS_OUT; SCOUT_OUTROOT (must be unique per alpha)

.libPaths(strsplit(Sys.getenv('SCOUT_LIB',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/software/R/R-4.4.2/library'), ':')[[1]])
ROOT <- Sys.getenv('SCOUT_ROOT',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/1_tree_robustness')

SCRIPTS <- file.path(ROOT, 'scripts')
N_REP   <- as.integer(Sys.getenv('SCOUT_NREP', '25'))
# Kept as a STRING: it is part of the simulated filenames, and as.numeric would turn 0.50 into 0.5
# and 1 into 1 in ways that may or may not match what simulate_test_data wrote.
ALPHA   <- Sys.getenv('SCOUT_ALPHA_G2', '3')

DATA_G2  <- Sys.getenv('SCOUT_DATA_G2', file.path(ROOT, 'data', '2_branch_lengths', 'affine_ng30'))
JOBS_OUT <- Sys.getenv('SCOUT_JOBS_OUT', file.path(SCRIPTS, sprintf('jobs_2_branch_lengths_v2_a%s.csv', ALPHA)))
OUTROOT  <- Sys.getenv('SCOUT_OUTROOT', file.path(ROOT, 'output', sprintf('2_branch_lengths_v2_a%s', ALPHA)))

man_f <- file.path(DATA_G2, 'replicate_manifest.csv')
if (!file.exists(man_f))
    stop(sprintf('%s not found -- simulate this alpha first (run_local_g2_alpha.sh sim %s).',
                 man_f, ALPHA))

man <- read.csv(man_f, stringsAsFactors = FALSE)
if (max(man$replicate) < N_REP)
    stop(sprintf('only %d replicates simulated in %s but SCOUT_NREP=%d.',
                 max(man$replicate), DATA_G2, N_REP))
man <- man[man$replicate <= N_REP, ]

# See ALPHA CONSISTENCY in the header. Every data file must name the alpha this manifest claims.
tok <- sprintf('_alpha_%s_sigma_', ALPHA)
bad <- c(man$counts_file, man$evf_file)
bad <- bad[!grepl(tok, basename(bad), fixed = TRUE)]
if (length(bad))
    stop(sprintf('SCOUT_ALPHA_G2=%s but %d data file(s) in %s were simulated at a different alpha, e.g.\n  %s\n',
                 ALPHA, length(bad), DATA_G2, basename(bad[1])),
         'The data directory and the alpha disagree -- fix SCOUT_DATA_G2 or SCOUT_ALPHA_G2.')

jobs <- do.call(rbind, lapply(c('counts', 'evf'), function(mt) data.frame(
    goal        = '2_branch_lengths',
    matrix_type = mt,
    replicate   = man$replicate,
    tree_bl     = man$tree_bl,
    tree_unit   = man$tree_unit,
    data_file   = if (mt == 'counts') man$counts_file else man$evf_file,
    alpha       = as.numeric(ALPHA),
    stringsAsFactors = FALSE)))
jobs <- jobs[order(jobs$matrix_type, jobs$replicate), ]
jobs <- cbind(task_id = seq_len(nrow(jobs)), jobs)

# ---- validate ----------------------------------------------------------------------------------
check <- function(label, paths) {
    missing <- unique(paths[!file.exists(paths)])
    cat(sprintf('%-24s %4d refs, %d missing %s\n', label, length(paths), length(missing),
                if (length(missing)) paste0('<- ', head(missing, 3), collapse = ' ') else ''))
    length(missing) == 0
}
ok <- c(check('true-BL trees', jobs$tree_bl),
        check('unit-BL trees', jobs$tree_unit),
        check('data matrices', jobs$data_file))

stopifnot(all(ok), identical(jobs$task_id, seq_len(nrow(jobs))),
          !any(duplicated(jobs[, c('matrix_type', 'replicate')])))

# The output root must be this alpha's own, or 02_run_task.r would overwrite another alpha's
# results: its filename stems carry matrix_type, arm and replicate but never alpha.
existing <- list.files(OUTROOT, pattern = '_accuracy\\.csv$', recursive = TRUE)
if (length(existing) && !identical(Sys.getenv('SCOUT_G2_FORCE', 'FALSE'), 'TRUE'))
    stop(sprintf('%s already holds %d result file(s). Output stems carry no alpha, so this would ',
                 OUTROOT, length(existing)),
         'overwrite them. Use a per-alpha SCOUT_OUTROOT, or SCOUT_G2_FORCE=TRUE to redo deliberately.')

write.csv(jobs, JOBS_OUT, row.names = FALSE)

cat(sprintf('\nalpha %s | %d replicates | %d tasks (each fits BOTH trees)\n',
            ALPHA, N_REP, nrow(jobs)))
print(as.data.frame(table(matrix_type = jobs$matrix_type)), row.names = FALSE)
cat(sprintf('  counts tasks %d-%d | evf tasks %d-%d\n',
            min(jobs$task_id[jobs$matrix_type == 'counts']), max(jobs$task_id[jobs$matrix_type == 'counts']),
            min(jobs$task_id[jobs$matrix_type == 'evf']),    max(jobs$task_id[jobs$matrix_type == 'evf'])))
cat(sprintf('\nmanifest: %s\noutput  : %s\ndata    : %s\n', JOBS_OUT, OUTROOT, DATA_G2))
