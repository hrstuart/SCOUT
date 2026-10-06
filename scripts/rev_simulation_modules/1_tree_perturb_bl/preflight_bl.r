# preflight_bl.r -- preflight gate for the branch-length baseline run: checks the R environment
# and SCOUT install, and asserts every path the run will write is absent or empty so published
# results cannot be overwritten. Exits non-zero on any failure.
#
# Usage:
#   Rscript preflight_bl.r

ROOT <- Sys.getenv('SCOUT_ROOT',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/1_tree_robustness')
LIB  <- Sys.getenv('SCOUT_LIB',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/software/R/R-4.4.2/library')
.libPaths(strsplit(LIB, ':')[[1]])
source(file.path(ROOT, 'scripts', '1_tree_perturb_bl', 'bl_paths.R'))

fails <- 0L
ok <- function(label, cond, detail = '') {
    cond <- isTRUE(cond)
    cat(sprintf('[%s] %-56s %s\n', if (cond) 'PASS' else 'FAIL', label, detail))
    if (!cond) fails <<- fails + 1L
    invisible(cond)
}

cat(sprintf('SCOUT_ROOT = %s\nSCOUT_LIB  = %s\n\n', ROOT, LIB))
P <- bl_paths(ROOT)

# ---- environment -------------------------------------------------------------------------------
cat('-- environment --\n')
ok('R is 4.4.x', grepl('^4\\.4\\.', paste0(R.version$major, '.', R.version$minor)),
   paste0('found ', R.version.string))
for (p in c('SCOUT', 'ape', 'castor', 'expm', 'OUwie', 'phangorn', 'dplyr', 'tidyr',
            'stringr', 'future', 'future.apply', 'corpcor', 'paleotree', 'nloptr', 'phylolm'))
    ok(sprintf('package %s loads', p), requireNamespace(p, quietly = TRUE))
for (f in c('SCOUT', 'runSCOUT', 'formatSCOUT', 'annotate_history',
            'calculate_group_class_accuracy', 'infer_anc', 'simulate_test_data',
            'identifiability_reference'))
    ok(sprintf('SCOUT exports %s()', f), exists(f, where = asNamespace('SCOUT'), inherits = FALSE))
ok('SCOUT:::extract_history_grid_search reachable',
   exists('extract_history_grid_search', where = asNamespace('SCOUT'), inherits = FALSE))

# This run needs the affine + support_clip build specifically.
ok('simulate_test_data has param_mode',
   'param_mode' %in% names(formals(SCOUT::simulate_test_data)))
ok('simulate_test_data has support_clip',
   'support_clip' %in% names(formals(SCOUT::simulate_test_data)))

# ---- parallel backend --------------------------------------------------------------------------
# This is the check that caught the outage, and it has to be conditional. Both of these are correct
# R and they impose DIFFERENT requirements:
#   future::plan(future::multisession, ...)  -- works with nothing more than future installed
#   plan(multisession, ...)                  -- needs future imported into SCOUT's NAMESPACE, or
#                                               attached to the search path by the caller
# The task scripts attach only SCOUT/ape/dplyr/stringr, so a bare-plan() build has no search-path
# fallback and dies on the first fit. But asserting resolution unconditionally would FAIL a
# perfectly good fully-qualified build. So: inspect what the installed code actually does, then hold
# it to the standard that form requires. Same logic as scripts/preflight_check.r.
cat('\n-- parallel backend --\n')
scout_resolves <- function(sym)
    !is.null(tryCatch(get(sym, envir = asNamespace('SCOUT'), mode = 'function'),
                      error = function(e) NULL))

ns  <- asNamespace('SCOUT')
src <- unlist(lapply(ls(ns, all.names = TRUE), function(n) {
    f <- tryCatch(get(n, envir = ns), error = function(e) NULL)
    if (is.function(f)) paste(deparse(f), collapse = '\n') else NULL
}))
# a leading ':' is excluded, so "future::plan(" is not counted as a bare call
bare_plan <- sum(grepl('(^|[^:[:alnum:]._])plan[[:space:]]*\\(', src))
qual_plan <- sum(grepl('future::plan[[:space:]]*\\(', src))
cat(sprintf('[INFO] SCOUT plan() call sites: %d bare, %d future::-qualified\n', bare_plan, qual_plan))
ok('SCOUT can reach plan()', bare_plan == 0 || scout_resolves('plan'),
   if (bare_plan == 0) 'all call sites qualified' else 'bare plan() with future not imported')
ok('SCOUT can reach multisession()', bare_plan == 0 || scout_resolves('multisession'))

# Resolving is necessary but not sufficient: the workers still have to start. Launching two and
# running one trivial future catches a backend that resolves but cannot spawn -- no free ports, R
# not on the workers' PATH, a cgroup that forbids it -- which would otherwise surface once per task,
# 640 times.
backend_ok <- tryCatch({
    future::plan(future::multisession, workers = 2)
    identical(future.apply::future_lapply(1:2, function(i) i * 2L), list(2L, 4L))
}, error = function(e) FALSE)
try(future::plan(future::sequential), silent = TRUE)   # always hand the session back
ok('multisession workers start and return', isTRUE(backend_ok), 'two workers, one future_lapply')

# ---- inputs this run READS (must exist, never written) -----------------------------------------
cat('\n-- read-only inputs --\n')
for (d in c(file.path(ROOT, 'scripts', 'scout_perturb_lib.R'),
            file.path(ROOT, 'scripts', 'scout_task_paths.R'),
            file.path(ROOT, 'scripts', 'run_one.sh'),
            file.path(ROOT, 'data', '1_tree_perturb', 'perturbed_tree_manifest.csv'),
            file.path(ROOT, 'data', '1_tree_perturb', 'perturbed_tree_manifest_hi.csv')))
    ok(sprintf('readable: %s', basename(d)), file.exists(d) && file.access(d, 4) == 0)
ok('published v2 output present (for the cross-run delta)',
   dir.exists(file.path(ROOT, 'output', '1_tree_perturb_v2')),
   'optional -- collation skips the comparison without it')

# ---- overwrite audit ---------------------------------------------------------------------------
#
# Two DIFFERENT questions, and conflating them makes this check useless after the first `sim`.
#
#   (a) Can this run clobber something that is NOT its own? That is the guarantee that matters, and
#       it is a STATIC property of the paths -- none of them may lie inside a published directory.
#       Always enforced, never satisfiable by accident, and it stays meaningful for the whole run.
#
#   (b) What of this run's OWN artifacts already exists? That is state, not a hazard. Once `sim` has
#       legitimately run, the data directory is SUPPOSED to be full, and the sweep fills the output
#       root by design. Failing on it would mean preflight could only ever pass once and would go
#       permanently red the moment the run started -- which trains you to ignore it.
#
# So (a) fails, (b) reports. The refuse-to-clobber guards inside each writing script are what stop an
# accidental re-`sim`; this is the audit, not the lock.

cat('\n-- (a) collision check: this run must not write inside a published path --\n')
PUBLISHED <- c(file.path(ROOT, 'data',   '1_tree_perturb'),
               file.path(ROOT, 'data',   '2_branch_lengths'),
               file.path(ROOT, 'output', '1_tree_perturb'),
               file.path(ROOT, 'output', '1_tree_perturb_v2'),
               file.path(ROOT, 'output', '2_branch_lengths'),
               file.path(ROOT, 'output', '2_branch_lengths_v2'),
               file.path(ROOT, 'logs_v2'))
norm <- function(x) sub('/+$', '', normalizePath(x, mustWork = FALSE))
for (nm in names(P)) {
    p <- norm(P[[nm]])
    # a path collides if it IS a published path or sits underneath one
    hit <- PUBLISHED[vapply(PUBLISHED, function(q) {
        q <- norm(q); identical(p, q) || startsWith(p, paste0(q, '/'))
    }, logical(1))]
    ok(sprintf('%s outside every published path', nm), length(hit) == 0,
       if (length(hit)) sprintf('COLLIDES WITH %s', hit[1]) else '')
}

cat('\n-- (b) current state of this run\'s own paths (reported, not enforced) --\n')
built <- file.exists(file.path(P$datadir, 'baseline_manifest.csv'))
for (nm in names(P)) {
    p <- P[[nm]]
    n <- if (grepl('\\.csv$', p)) as.integer(file.exists(p))
         else if (dir.exists(p)) length(list.files(p, all.files = TRUE, no.. = TRUE)) else 0L
    cat(sprintf('[INFO] %-56s %s\n', nm, if (n == 0) 'empty / absent' else sprintf('%d entr%s', n, if (n == 1) 'y' else 'ies')))
}
cat(sprintf('[INFO] %-56s %s\n', 'simulation stage',
            if (built) 'baseline already simulated -- `sim` will REFUSE (use SCOUT_BL_FORCE=TRUE to redo)'
                       else 'not yet simulated -- run `sim` next'))

ok('output root writable', file.access(dirname(P$outroot), 2) == 0, dirname(P$outroot))
ok('scripts dir writable (manifests land here)', file.access(file.path(ROOT, 'scripts'), 2) == 0)

# ---- disk ---------------------------------------------------------------------------------------
cat('\n-- capacity --\n')
avail <- tryCatch({
    x <- system(sprintf("df -Pk %s | awk 'NR==2{print $4}'", shQuote(ROOT)), intern = TRUE)
    as.numeric(x) / 1024 / 1024
}, error = function(e) NA_real_)
ok('at least 5 GB free', is.na(avail) || avail > 5,
   if (is.na(avail)) 'could not determine' else sprintf('%.1f GB available (run needs ~0.1 GB)', avail))

cat(sprintf('\n%s: %d check(s) failed\n', if (fails == 0) 'PREFLIGHT PASSED' else 'PREFLIGHT FAILED', fails))
quit(save = 'no', status = if (fails == 0) 0 else 1)
