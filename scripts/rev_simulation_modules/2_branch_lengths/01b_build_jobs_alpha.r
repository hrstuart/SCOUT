# Goal 2 (2_branch_lengths) step 01b: build the job manifest for ONE alpha.
#
# WHY THIS EXISTS RATHER THAN build_jobs_v2.r.
#
# build_jobs_v2.r builds goal 1 AND goal 2 together and writes both to fixed filenames
# (jobs_1_tree_perturb_v2.csv, jobs_2_branch_lengths_v2.csv). Running it to add an alpha would
# clobber the manifests the published alpha-3 run was produced from, and would rebuild goal 1 for no
# reason. This script does the goal-2 half only, reads its data directory and alpha from the
# environment, and writes wherever SCOUT_JOBS_OUT points -- so each alpha gets its own manifest and
# nothing existing is touched.
#
# The goal-2 block below is a faithful copy of build_jobs_v2.r's: same columns, same ordering
# (matrix_type then replicate, so counts tasks come first), same task_id contract, same duplicate
# checks. It adds one guard build_jobs_v2.r does not need, because that script only ever ran at one
# alpha -- see ALPHA CONSISTENCY below.
#
# ALPHA CONSISTENCY. `alpha` is written into the manifest as a column, but the data files are found
# through the replicate manifest, whose filenames encode the alpha they were simulated at
# (rep001_alpha_0.5_sigma_1_OUEVF_counts.csv). Pointing SCOUT_DATA_G2 at one alpha's directory while
# setting SCOUT_ALPHA_G2 to another would produce a manifest that fits alpha-3 data and labels every
# row alpha 0.5 -- wrong in a way nothing downstream could detect, because every path would exist.
# The two are therefore cross-checked against each other before anything is written.
#
# OUTPUT PATHS CARRY NO ALPHA. 02_run_task.r writes <out_root>/<matrix_type>/<matrix>_<arm>_repNNN,
# with no alpha anywhere in the stem, so two alphas sharing an output root would silently overwrite
# each other. That is handled by giving each alpha its own SCOUT_OUTROOT (run_local_g2_alpha.sh
# does), and asserted here so a hand-run invocation cannot get it wrong.

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
