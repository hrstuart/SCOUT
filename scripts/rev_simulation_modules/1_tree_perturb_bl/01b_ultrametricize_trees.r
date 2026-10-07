# 01b_ultrametricize_trees.r -- Goal 1, ultrametric-restored set, step 01b: refit the branch
# lengths of every perturbed tree from 01_make_trees.r so all tips share the baseline depth,
# topology unchanged. Removes the tip-depth variance signal that biases fits against BM.
#
# Usage:
#   Rscript 01b_ultrametricize_trees.r      # optional env: SCOUT_UM_METHOD = nnls (default) | extend | chronos
#                                           #               SCOUT_BL_TAG (default '_um')

.libPaths(strsplit(Sys.getenv('SCOUT_LIB',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/software/R/R-4.4.2/library'), ':')[[1]])
ROOT <- Sys.getenv('SCOUT_ROOT',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/1_tree_robustness')
source(file.path(ROOT, 'scripts', 'scout_perturb_lib.R'))
source(file.path(ROOT, 'scripts', '1_tree_perturb_bl', 'bl_paths.R'))
suppressPackageStartupMessages({library(ape); library(phangorn); library(phytools)})

# The SOURCE set: the existing non-ultrametric bl run. Read-only throughout.
SRC <- bl_paths(ROOT, tag = '')
# The TARGET set: whatever SCOUT_BL_TAG points at (the runner sets it to '_um').
TAG <- Sys.getenv('SCOUT_BL_TAG', '_um')
if (!nzchar(TAG))
    stop('SCOUT_BL_TAG is empty, so this script would overwrite the non-ultrametric run in place. ',
         'Set SCOUT_BL_TAG (the runner sets it to "_um").')
DST <- bl_paths(ROOT, tag = TAG)

METHOD  <- Sys.getenv('SCOUT_UM_METHOD', 'nnls')
ARMS    <- strsplit(Sys.getenv('SCOUT_UM_ARMS', 'shuffle,nni,collapse'), ',')[[1]]
N_REP   <- as.integer(Sys.getenv('SCOUT_NREP', '15'))
stopifnot(METHOD %in% c('nnls', 'extend', 'chronos'))

# Tolerances. FLOOR_LEN matches resolve_and_floor() in scout_perturb_lib.R, so a zero-length edge
# is handled the same way it is everywhere else in this project.
#
# Neither ultrametricity nor the height comes out EXACT, because NNLS returns zero-length edges
# (8 on nni@10%, up to 163 on collapse@90%) and flooring them adds up to FLOOR_LEN to some
# root-to-tip paths and not others. Measured across all six nni/collapse conditions, the residual
# is cv_var_bm <= 8.8e-7 and |H - H_TARGET| <= 7e-7. The tolerances below sit two orders above the
# worst observed value and are still four to five orders BELOW the smallest heteroscedasticity this
# run exists to remove (cv_var_bm = 0.111 at nni@10%), so nothing scientific rides on the slack.
FLOOR_LEN <- 1e-7
CV_TOL    <- 1e-4
H_TOL     <- 1e-5

guard_empty(DST$treedir, 'ultrametric-restored tree directory')
dir.create(DST$treedir, recursive = TRUE, showWarnings = FALSE)

# --- source material ----------------------------------------------------------------------------
src_man_f <- file.path(SRC$datadir, 'perturbed_tree_manifest_bl.csv')
src_base_f <- file.path(SRC$datadir, 'baseline_manifest.csv')
for (f in c(src_man_f, src_base_f))
    if (!file.exists(f))
        stop(sprintf('%s not found -- run the non-ultrametric bl pipeline (sim) first.', f))

man <- read.csv(src_man_f, stringsAsFactors = FALSE)
man <- man[man$arm %in% ARMS & man$replicate <= N_REP, , drop = FALSE]
man <- man[order(man$arm, man$intensity, man$replicate), ]
if (!nrow(man)) stop('no rows left in the source manifest after filtering by SCOUT_UM_ARMS/SCOUT_NREP.')

bm <- read.csv(src_base_f, stringsAsFactors = FALSE)
base_tree <- read.tree(bm$tree_file[1])
NT <- length(base_tree$tip.label)
H_TARGET <- max(node.depth.edgelength(base_tree)[seq_len(NT)])
stopifnot(is.ultrametric(base_tree))
if (abs(H_TARGET - 4.738808) > 1e-5)
    stop(sprintf('baseline height %.6f != 4.738808 -- this is not the tree the bl run used.', H_TARGET))

cat(sprintf('source trees : %s\n', SRC$treedir))
cat(sprintf('target trees : %s\n', DST$treedir))
cat(sprintf('baseline     : %d tips, H = %.6f (ultrametric)\n', NT, H_TARGET))
cat(sprintf('method       : %s | arms: %s | replicates: %d | %d trees\n\n',
            METHOD, paste(ARMS, collapse = ','), N_REP, nrow(man)))

# --- the transform --------------------------------------------------------------------------------
tip_depths <- function(phy) node.depth.edgelength(phy)[seq_len(length(phy$tip.label))]
cv_var_bm  <- function(phy) { d <- tip_depths(phy); stats::sd(d) / mean(d) }

#' Restore a common tip depth without touching the topology.
#'
#' Returns the tree plus the method actually used -- an already-ultrametric tree short-circuits to
#' 'noop' so the shuffle arm comes through byte-identical rather than merely numerically close,
#' and a failed NNLS falls back to `extend` rather than taking the whole set down.
#' force.ultrametric prints a fixed-text Note with cat(), which suppressMessages cannot reach, so
#' it is captured and discarded rather than repeated 90 times in the log. The advice in it -- that
#' this is a coercion and not rate smoothing -- is the point: we are not dating a tree, we are
#' removing a tip-depth gradient that is an artefact of the perturbation, and `chronos` (which IS
#' formal rate smoothing) is available via SCOUT_UM_METHOD for the sensitivity check.
quietly <- function(expr) {
    out <- NULL
    invisible(utils::capture.output(out <- suppressWarnings(suppressMessages(expr))))
    out
}

ultrametricize <- function(phy, method) {
    if (is.ultrametric(phy, tol = 1e-8)) return(list(tree = phy, method = 'noop'))
    um <- tryCatch(quietly(switch(method,
            nnls    = phytools::force.ultrametric(phy, method = 'nnls'),
            extend  = phytools::force.ultrametric(phy, method = 'extend'),
            chronos = ape::chronos(phy, quiet = TRUE))),
        error = function(e) NULL)
    if (is.null(um) || is.null(um$edge.length) || any(!is.finite(um$edge.length))) {
        um <- quietly(phytools::force.ultrametric(phy, method = 'extend'))
        return(list(tree = um, method = paste0(method, '_failed->extend')))
    }
    list(tree = um, method = method)
}

#' Floor, rescale to the baseline height, floor again.
#'
#' The order matters. Flooring AFTER the rescale lengthens whichever root-to-tip paths happen to
#' contain floored edges and leaves the height up to 3.2e-6 over target; flooring first means the
#' rescale is computed on the depths the tree will actually have, and the second pass only catches
#' edges the rescale itself pushed under the floor. Measured, that brings |H - H_TARGET| from
#' 3.2e-6 down to <= 7e-7 with no cost to cv_var_bm.
finish <- function(phy, H_target) {
    phy$edge.length[phy$edge.length < FLOOR_LEN] <- FLOOR_LEN
    phy$edge.length <- phy$edge.length * (H_target / max(tip_depths(phy)))
    phy$edge.length[phy$edge.length < FLOOR_LEN] <- FLOOR_LEN
    phy
}

rows <- lapply(seq_len(nrow(man)), function(i) {
    r   <- man[i, ]
    src <- read.tree(r$tree_file)

    res <- ultrametricize(src, METHOD)
    um  <- finish(res$tree, H_TARGET)

    # Topology is the perturbation; if it moved, this tree is not a controlled version of its twin.
    rf_vs_src <- as.numeric(phangorn::RF.dist(src, um))
    if (rf_vs_src != 0)
        stop(sprintf('%s: RF against its own source tree is %g, not 0 -- the topology moved.',
                     basename(r$tree_file), rf_vs_src))
    if (!setequal(um$tip.label, src$tip.label))
        stop(sprintf('%s: tip labels changed.', basename(r$tree_file)))

    cv <- cv_var_bm(um)
    if (!is.finite(cv) || cv > CV_TOL)
        stop(sprintf('%s: cv_var_bm is %.3e after restoration (tol %.0e) -- not ultrametric.',
                     basename(r$tree_file), cv, CV_TOL))
    if (any(um$edge.length <= 0))
        stop(sprintf('%s: non-positive edge length after restoration.', basename(r$tree_file)))

    # How much of the perturbed covariance structure survived the ultrametric constraint.
    pc <- suppressWarnings(stats::cor(as.vector(as.dist(cophenetic(src))),
                                      as.vector(as.dist(cophenetic(um)))))

    out_f <- file.path(DST$treedir,
                       sub('\\.nwk$', '_um.nwk', basename(r$tree_file)))
    ape::write.tree(um, out_f)

    data.frame(arm = r$arm, intensity = r$intensity, replicate = r$replicate, seed = r$seed,
               tree_file = out_f, k = r$k, rf_dist = r$rf_dist,
               src_tree_file = r$tree_file, um_method = res$method,
               height_src = max(tip_depths(src)), height_um = max(tip_depths(um)),
               cv_var_bm_src = cv_var_bm(src), cv_var_bm_um = cv,
               rf_vs_src = rf_vs_src, patristic_cor = pc,
               stringsAsFactors = FALSE)
})
um_man <- do.call(rbind, rows)

# --- set-level assertions -------------------------------------------------------------------------
stopifnot(nrow(um_man) == nrow(man),
          all(file.exists(um_man$tree_file)),
          !any(duplicated(um_man$tree_file)),
          all(um_man$rf_vs_src == 0),
          all(um_man$cv_var_bm_um <= CV_TOL),
          all(abs(um_man$height_um - H_TARGET) < H_TOL))

# The shuffle arm must have come through untouched -- it is the positive control for this whole run.
sh <- um_man[um_man$arm == 'shuffle', ]
if (nrow(sh)) {
    if (!all(sh$um_method == 'noop'))
        stop('shuffle trees were modified; they are already ultrametric and must pass through as no-ops.')
    same <- mapply(function(a, b) identical(tools::md5sum(a)[[1]], tools::md5sum(b)[[1]]),
                   sh$tree_file, sh$src_tree_file)
    cat(sprintf('shuffle no-op check: %d of %d trees byte-identical to their source\n',
                sum(same), length(same)))
    if (!all(same))
        stop('a shuffle tree differs from its source despite the no-op path -- investigate before running.')
}

# 02_build_jobs.r resolves the data files and the baseline tree from THIS directory's
# baseline_manifest.csv. It is copied rather than regenerated, so the um run reads the exact same
# counts/EVF matrices and the exact same baseline tree as the run it is being compared against.
file.copy(src_base_f, file.path(DST$datadir, 'baseline_manifest.csv'), overwrite = TRUE)

# Named as 02_build_jobs.r expects; it lives in a different directory from the non-ultrametric one
# and every tree it points at carries the _um suffix, so nothing can be confused for the other set.
manfile <- file.path(DST$datadir, 'perturbed_tree_manifest_bl.csv')
write.csv(um_man, manfile, row.names = FALSE)

# --- report ---------------------------------------------------------------------------------------
cat('\n-- tip-depth heteroscedasticity removed --\n')
agg <- aggregate(cbind(height_src, height_um, cv_var_bm_src, cv_var_bm_um, patristic_cor) ~
                 arm + intensity, data = um_man, FUN = mean)
agg <- agg[order(agg$arm, agg$intensity), ]
cat(sprintf('%-9s %5s %9s %9s %12s %12s %10s\n',
            'arm', 'int', 'H_src', 'H_um', 'cv_bm_src', 'cv_bm_um', 'patr_cor'))
for (i in seq_len(nrow(agg)))
    cat(sprintf('%-9s %5.2f %9.3f %9.3f %12.4f %12.2e %10.4f\n',
                agg$arm[i], agg$intensity[i], agg$height_src[i], agg$height_um[i],
                agg$cv_var_bm_src[i], agg$cv_var_bm_um[i], agg$patristic_cor[i]))

cat(sprintf('\nmethods used: %s\n',
            paste(sprintf('%s=%d', names(table(um_man$um_method)), table(um_man$um_method)),
                  collapse = ' ')))
cat(sprintf('%d trees written to %s\n', nrow(um_man), DST$treedir))
cat(sprintf('manifest: %s\n', manfile))
cat(sprintf('baseline manifest copied from %s (same data, same baseline tree)\n', src_base_f))
