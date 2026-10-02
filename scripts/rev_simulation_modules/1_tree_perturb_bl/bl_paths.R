# Path and completeness helpers for the branch-length baseline run (1_tree_perturb_bl).
#
# WHY THIS FILE EXISTS AT ALL. The published pipeline resolves task paths through task_dir()/
# task_stems() in scripts/scout_task_paths.R. Those two functions key off the goal name and have no
# notion of alpha -- this run sweeps alpha, so two alphas of the same (arm, intensity, replicate)
# would resolve to the SAME output path and silently overwrite each other. Rather than edit a file
# the published goal-1 and goal-2 runs depend on, this run carries its own resolvers.
#
# task_complete() is NOT redefined here. It is sourced from scout_task_paths.R and reused verbatim,
# so "done" means exactly what it means for every other run in this project -- every output CSV
# present AND non-empty, derived from files on disk and nothing else. That is the predicate that
# makes `resume` correct after a crash, and there must be only one of it.
#
# Deliberately BASE R ONLY, with no library() calls: 05_status.r refreshes on a 60s loop and cannot
# afford library(SCOUT), which costs ~2m35s over DARTFS.

ROOT <- Sys.getenv('SCOUT_ROOT',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/1_tree_robustness')

# Reused, never redefined -- see the header.
source(file.path(ROOT, 'scripts', 'scout_task_paths.R'))


#' Alpha rendered exactly as it appears in a simulate_test_data() filename.
#'
#' SCOUT builds its output names with sprintf('%s', a) on the numeric alpha, so 0.25 becomes "0.25"
#' and 3 becomes "3" (not "3.0"). Every path that has to line up with a simulated file -- the data
#' file a task reads, the directory its results land in -- goes through this one function, so the
#' manifest and the data on disk cannot drift apart on formatting.
alpha_str <- function(a) sprintf('%s', as.numeric(a))


#' Directory tag for one alpha, e.g. 0.5 -> "alpha0.5".
alpha_tag <- function(a) paste0('alpha', alpha_str(a))


#' Directory a task writes into.
#'
#' <out_root>/<layer>/alpha<a>/<arm>. Alpha sits in the DIRECTORY rather than the filename so
#' 04_collate.r can glob a whole alpha at once, and so the stems stay identical to the published
#' run's -- which is what lets the cross-run comparison against output/1_tree_perturb_v2/ match
#' files up by name.
task_dir_bl <- function(job, out_root) {
    file.path(out_root, job$layer, alpha_tag(job$alpha), job$arm)
}


#' Filename stem a task writes.
#'
#' One uniform form for every row including the unperturbed baseline, which is recorded as
#' arm='baseline', intensity=0, replicate=0 and so lands at baseline_i000_rep000. Keeping the
#' baseline in the same shape as the arms means it flows through the manifest, the tracker, the
#' skip-if-done check and the collation with no special-casing anywhere.
task_stems_bl <- function(job) {
    sprintf('%s_i%03d_rep%03d', job$arm, round(job$intensity * 100), job$replicate)
}


#' Every directory and manifest this run is allowed to write.
#'
#' Single source of truth for the overwrite guards: the preflight asserts none of these exists yet,
#' and each writing script refuses to start if its own target is already populated. Listing them in
#' one place is what makes "this run cannot touch the published results" checkable rather than
#' merely intended.
#'
#' `tag` suffixes the whole path set at once, so a VARIANT of this run gets its own data directory,
#' output tree, log directory and job manifests from a single switch and 02/03/04/05 need no edits
#' to serve it. It defaults to SCOUT_BL_TAG, which is unset (and therefore '') for the original
#' non-ultrametric run -- with an empty tag every path below is exactly what it always was.
#' SCOUT_BL_TAG='_um' selects the ultrametric-restored set built by 01b_ultrametricize_trees.r.
#'
#' Because the default reads the environment, a script that must pin one specific set regardless of
#' how it was invoked should pass `tag` explicitly (01b_ultrametricize_trees.r passes tag = '' to
#' read the non-ultrametric source trees while writing into the tagged target).
bl_paths <- function(root = ROOT, tag = Sys.getenv('SCOUT_BL_TAG', '')) {
    goal <- paste0('1_tree_perturb_bl', tag)
    list(
        datadir  = file.path(root, 'data',   goal),
        treedir  = file.path(root, 'data',   goal, 'perturbed_trees'),
        outroot  = file.path(root, 'output', goal),
        logdir   = file.path(root, paste0('logs_bl', tag)),
        jobs_screen = file.path(root, 'scripts', sprintf('jobs_%s_screen.csv', goal)),
        jobs_sweep  = file.path(root, 'scripts', sprintf('jobs_%s_sweep.csv',  goal)))
}


#' Refuse to write into a directory that already holds results.
#'
#' The whole point of this run is that it cannot disturb the published one, so a populated target is
#' treated as an error and not as something to merge into. SCOUT_BL_FORCE=TRUE is the deliberate
#' override for an intentional redo.
guard_empty <- function(dir, what) {
    if (!dir.exists(dir)) return(invisible(TRUE))
    n <- length(list.files(dir, all.files = TRUE, no.. = TRUE))
    if (n == 0) return(invisible(TRUE))
    if (isTRUE(as.logical(Sys.getenv('SCOUT_BL_FORCE', 'FALSE')))) {
        cat(sprintf('SCOUT_BL_FORCE=TRUE -- overwriting %s (%d existing entries in %s)\n', what, n, dir))
        return(invisible(TRUE))
    }
    stop(sprintf('%s already exists and is not empty (%d entries):\n  %s\nRefusing to overwrite. Set SCOUT_BL_FORCE=TRUE to redo it deliberately.',
                 what, n, dir), call. = FALSE)
}
