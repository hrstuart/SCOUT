# Goal 1 on a BRANCH-LENGTH baseline (1_tree_perturb_bl) step 04: collate per-task CSVs.
#
# Produces three things, in increasing order of what the run exists to answer:
#
#   1. ROLLUP -- mean/sd/se accuracy per (layer, alpha, arm, intensity), plus per-class recall and
#      the pooled confusion. The descriptive layer, mirroring 1_tree_perturb/04_collate.r.
#
#   2. ARM-vs-BASELINE DELTA -- mean(arm) - baseline(SAME alpha, SAME layer). The primary
#      within-run result: now that the branch lengths are real and tip-depth heteroscedasticity is
#      gone, does degrading the topology actually cost anything? The baseline is a single
#      unperturbed fit per (alpha, layer), so the uncertainty quoted is the arm's own SE across
#      replicates -- the baseline contributes no replicate variance and it would be wrong to imply
#      otherwise.
#
#   3. CROSS-RUN DELTA against a reference run -- this run minus the reference, paired by
#      replicate, at every condition the two share. This is the contrast the whole design was built
#      for: the perturbed topologies in the two runs are IDENTICAL (01_make_trees.r asserts it via
#      RF), the tip states are identical, the seeds are identical, so the difference isolates one
#      variable and nothing else. Only alphas present in both runs are compared.
#
#      SCOUT_V2_OUTROOT picks the reference and SCOUT_CROSS_TAG names it in the output. For the
#      original run the reference is output/1_tree_perturb_v2 (tag 'v2') and the isolated variable
#      is the BRANCH LENGTHS. For the ultrametric-restored run the reference is
#      output/1_tree_perturb_bl (tag 'noum') and the isolated variable is TIP-DEPTH
#      HETEROSCEDASTICITY -- same topologies, same data, same seeds, differing only in whether the
#      perturbed tree was restored to a common tip depth before fitting.
#
# Safe to run while the sweep is still going -- it reports how many of the expected tasks are
# present and summarises whatever has landed. output/1_tree_perturb_v2/ is read STRICTLY read-only.

.libPaths(strsplit(Sys.getenv('SCOUT_LIB',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/software/R/R-4.4.2/library'), ':')[[1]])
ROOT <- Sys.getenv('SCOUT_ROOT',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/1_tree_robustness')
source(file.path(ROOT, 'scripts', '1_tree_perturb_bl', 'bl_paths.R'))
suppressPackageStartupMessages(library(dplyr))

P      <- bl_paths(ROOT)
OUTDIR <- Sys.getenv('SCOUT_OUTROOT', P$outroot)
SUMDIR <- file.path(OUTDIR, 'summaries')
dir.create(SUMDIR, recursive = TRUE, showWarnings = FALSE)

# What the summary files are called, and what the cross-run reference is called inside them. Both
# default to the original run's values, so an untagged invocation writes exactly the files it
# always did. The ultrametric-restored run overrides them so its summaries are self-describing
# rather than claiming to be a v2 comparison they are not.
PREFIX    <- Sys.getenv('SCOUT_SUMMARY_PREFIX', '1_tree_perturb_bl')
CROSS_TAG <- Sys.getenv('SCOUT_CROSS_TAG', 'v2')
sumf <- function(fmt, ...) file.path(SUMDIR, sprintf(fmt, ...))

# bind_rows, not rbind. The published v2 output tree holds two DIFFERENT schemas: the 900 array-task
# files, and the 6 alpha-screen files under baseline_alpha/ which carry (arm, layer, intensity,
# replicate, alpha, alphaH) instead of (goal, layer, arm, ..., task_id, elapsed_min). rbind() dies on
# "numbers of columns of arguments do not match" the moment both are present. bind_rows fills the
# gaps with NA, and every consumer below selects the columns it needs by name, so a file missing a
# column simply contributes NA there rather than taking the collation down.
read_all <- function(dir, pattern) {
    f <- list.files(dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
    f <- f[!grepl('/summaries/', f)]
    if (length(f) == 0) return(NULL)
    dplyr::bind_rows(lapply(f, read.csv, stringsAsFactors = FALSE))
}

acc <- read_all(OUTDIR, '_accuracy\\.csv$')
if (is.null(acc)) stop('No per-task accuracy files found under ', OUTDIR)

jobs_file <- Sys.getenv('SCOUT_JOBS', P$jobs_sweep)
n_expected <- if (file.exists(jobs_file)) nrow(read.csv(jobs_file)) else NA
overall <- acc %>% filter(metric_type == 'Overall')
cat(sprintf('%d of %s expected tasks present\n', nrow(overall), ifelse(is.na(n_expected), '?', n_expected)))

subtask <- overall %>%
    select(dataset, goal, layer, alpha, alphaH, arm, intensity, replicate, tree_height,
           n, n_correct, accuracy, ci_lower, ci_upper, elapsed_min) %>%
    arrange(layer, alpha, arm, intensity, replicate)
write.csv(subtask, sumf('%s_all_replicates.csv', PREFIX), row.names = FALSE)

# ---- 1. rollup ---------------------------------------------------------------------------------
rollup <- subtask %>%
    group_by(layer, alpha, alphaH, arm, intensity) %>%
    summarise(n_replicates = n(),
              mean_accuracy = mean(accuracy),
              sd_accuracy   = sd(accuracy),
              se_accuracy   = sd(accuracy) / sqrt(n()),
              min_accuracy  = min(accuracy),
              max_accuracy  = max(accuracy),
              mean_height   = mean(tree_height),
              mean_elapsed_min = mean(elapsed_min),
              .groups = 'drop') %>%
    arrange(layer, alpha, arm, intensity)
write.csv(rollup, sumf('%s_summary.csv', PREFIX), row.names = FALSE)

# ---- 2. arm vs its own baseline ----------------------------------------------------------------
base <- subtask %>% filter(arm == 'baseline') %>%
    select(layer, alpha, baseline_accuracy = accuracy)

delta <- rollup %>%
    filter(arm != 'baseline') %>%
    left_join(base, by = c('layer', 'alpha')) %>%
    mutate(delta_vs_baseline = mean_accuracy - baseline_accuracy,
           # t on the ARM's replicates only: the baseline is one unperturbed fit and contributes no
           # replicate variance, so this is a one-sample test of the arm against a fixed reference.
           t_stat = ifelse(se_accuracy > 0, delta_vs_baseline / se_accuracy, NA_real_),
           p_value = ifelse(is.na(t_stat), NA_real_,
                            2 * pt(-abs(t_stat), df = n_replicates - 1))) %>%
    select(layer, alpha, alphaH, arm, intensity, n_replicates, baseline_accuracy,
           mean_accuracy, se_accuracy, delta_vs_baseline, t_stat, p_value) %>%
    arrange(layer, alpha, arm, intensity)
write.csv(delta, sumf('%s_delta_vs_baseline.csv', PREFIX), row.names = FALSE)

# ---- 3. cross-run: bl vs the published v2 run --------------------------------------------------
V2DIR <- Sys.getenv('SCOUT_V2_OUTROOT', file.path(ROOT, 'output', '1_tree_perturb_v2'))
cross <- NULL
if (!dir.exists(V2DIR)) {
    cat(sprintf('\nno %s output at %s -- skipping the cross-run comparison\n', CROSS_TAG, V2DIR))
} else {
    v2 <- read_all(V2DIR, '_accuracy\\.csv$')
    if (is.null(v2)) {
        cat('\nno v2 accuracy files found -- skipping the cross-run comparison\n')
    } else {
        if (!'layer' %in% names(v2)) v2$layer <- 'counts'   # pre-affine runs were counts-only
        v2o <- v2 %>% filter(metric_type == 'Overall') %>%
            select(layer, alpha, arm, intensity, replicate, accuracy_v2 = accuracy)

        # Paired by replicate: replicate r of a condition is the SAME perturbed topology in both
        # runs (asserted in 01_make_trees.r), so the pairing is real and not just a shared index.
        cross <- subtask %>%
            filter(arm != 'baseline') %>%
            select(layer, alpha, arm, intensity, replicate, accuracy_bl = accuracy) %>%
            inner_join(v2o, by = c('layer', 'alpha', 'arm', 'intensity', 'replicate')) %>%
            mutate(delta_bl_minus_v2 = accuracy_bl - accuracy_v2)

        if (!nrow(cross)) {
            cat('\nno conditions shared with the v2 run (alphas differ?) -- nothing to compare\n')
            cross <- NULL
        } else {
            cross_sum <- cross %>%
                group_by(layer, alpha, arm, intensity) %>%
                summarise(n_paired = n(),
                          mean_bl  = mean(accuracy_bl),
                          mean_v2  = mean(accuracy_v2),
                          paired_mean_delta = mean(delta_bl_minus_v2),
                          paired_sd    = sd(delta_bl_minus_v2),
                          paired_se    = sd(delta_bl_minus_v2) / sqrt(n()),
                          p_value = tryCatch(t.test(delta_bl_minus_v2)$p.value,
                                             error = function(e) NA_real_),
                          .groups = 'drop') %>%
                arrange(layer, alpha, arm, intensity)
            # Name the reference columns after the run they actually came from, so a summary
            # cannot be misread as a comparison against something it never touched. The one
            # substitution covers accuracy_v2, mean_v2 and delta_bl_minus_v2, and with the default
            # CROSS_TAG='v2' it rewrites every name to itself -- the original run's files keep the
            # exact column names they have always had.
            retag <- function(d) { names(d) <- sub('_v2$', paste0('_', CROSS_TAG), names(d)); d }
            write.csv(retag(cross),     sumf('%s_vs_%s_replicates.csv', PREFIX, CROSS_TAG), row.names = FALSE)
            write.csv(retag(cross_sum), sumf('%s_vs_%s_summary.csv',    PREFIX, CROSS_TAG), row.names = FALSE)
        }
    }
}

# ---- per-class recall, and confusion -----------------------------------------------------------
sens <- acc %>% filter(metric_type == 'Sensitivity') %>%
    group_by(layer, alpha, arm, intensity, class) %>%
    summarise(mean_recall = mean(accuracy, na.rm = TRUE),
              mean_specificity = mean(specificity, na.rm = TRUE), .groups = 'drop')
write.csv(sens, sumf('%s_per_class.csv', PREFIX), row.names = FALSE)

conf <- read_all(OUTDIR, '_confusion\\.csv$')
if (!is.null(conf)) {
    conf_sum <- conf %>% group_by(layer, alpha, arm, intensity, truth, selected) %>%
        summarise(n = sum(Freq), .groups = 'drop')
    write.csv(conf_sum, sumf('%s_confusion.csv', PREFIX), row.names = FALSE)
}

# ---- report ------------------------------------------------------------------------------------
cat('\n== baseline accuracy per alpha (the alpha screen result) ==\n')
print(as.data.frame(subtask %>% filter(arm == 'baseline') %>%
                    select(layer, alpha, alphaH, accuracy, ci_lower, ci_upper, n)),
      row.names = FALSE, digits = 3)

cat('\n== arm vs its own baseline ==\n')
print(as.data.frame(delta), row.names = FALSE, digits = 3)

if (!is.null(cross)) {
    cat(sprintf('\n== this run - %s at matched conditions (identical topologies) ==\n', CROSS_TAG))
    print(as.data.frame(read.csv(sumf('%s_vs_%s_summary.csv', PREFIX, CROSS_TAG))),
          row.names = FALSE, digits = 3)
}

cat(sprintf('\nSummaries written to %s\n', SUMDIR))
