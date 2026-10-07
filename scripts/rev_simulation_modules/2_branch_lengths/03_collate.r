# 03_collate.r -- Goal 2, step 03: collate per-task CSVs into summary files: per-arm accuracy and
# the paired true-vs-unit difference with a paired t-test.
#
# Usage:
#   Rscript 03_collate.r                    # optional env: SCOUT_OUTROOT, SCOUT_JOBS

# SCOUT_LIB may be a colon-separated path LIST, so a private build (e.g. the support_clip
# SCOUT) can be prepended while its dependencies still resolve from the shared library.
.libPaths(strsplit(Sys.getenv('SCOUT_LIB',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/software/R/R-4.4.2/library'), ':')[[1]])
ROOT <- Sys.getenv('SCOUT_ROOT',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/1_tree_robustness')
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

OUTDIR <- Sys.getenv('SCOUT_OUTROOT', file.path(ROOT, 'output', '2_branch_lengths'))
SUMDIR <- file.path(OUTDIR, 'summaries')
dir.create(SUMDIR, recursive = TRUE, showWarnings = FALSE)

read_all <- function(pattern) {
    f <- list.files(OUTDIR, pattern = pattern, recursive = TRUE, full.names = TRUE)
    f <- f[!grepl('/summaries/', f)]
    if (length(f) == 0) return(NULL)
    do.call(rbind, lapply(f, read.csv, stringsAsFactors = FALSE))
}

acc <- read_all('_accuracy\\.csv$')
if (is.null(acc)) stop('No per-task accuracy files found under ', OUTDIR)

overall <- acc %>% filter(metric_type == 'Overall')
jobs_file <- Sys.getenv('SCOUT_JOBS', file.path(ROOT, 'scripts', 'jobs_2_branch_lengths.csv'))
n_expected <- if (file.exists(jobs_file)) nrow(read.csv(jobs_file)) * 2 else NA
cat(sprintf('%d of %s expected fits present (2 arms per task)\n',
            nrow(overall), ifelse(is.na(n_expected), '?', n_expected)))

# ---- per sub-task summary files (one per matrix type) -----------------------------------------
subtask <- overall %>%
    select(dataset, goal, matrix_type, arm, replicate, alpha, height, n, n_correct,
           accuracy, ci_lower, ci_upper, elapsed_min) %>%
    arrange(matrix_type, arm, replicate)

for (mt in unique(subtask$matrix_type)) {
    sub <- subtask %>% filter(matrix_type == mt)
    f <- file.path(SUMDIR, sprintf('summary_%s.csv', mt))
    write.csv(sub, f, row.names = FALSE)
    cat(sprintf('  %-6s  n=%3d rows -> %s\n', mt, nrow(sub), basename(f)))
}
write.csv(subtask, file.path(SUMDIR, '2_branch_lengths_all_replicates.csv'), row.names = FALSE)

# ---- paired analysis --------------------------------------------------------------------------
paired <- subtask %>%
    select(matrix_type, replicate, arm, accuracy) %>%
    pivot_wider(names_from = arm, values_from = accuracy) %>%
    filter(!is.na(true_bl), !is.na(unit_bl)) %>%
    mutate(diff = true_bl - unit_bl)

write.csv(paired, file.path(SUMDIR, '2_branch_lengths_paired.csv'), row.names = FALSE)

rollup <- paired %>%
    group_by(matrix_type) %>%
    summarise(n_replicates   = n(),
              mean_true_bl   = mean(true_bl),
              mean_unit_bl   = mean(unit_bl),
              mean_diff      = mean(diff),
              se_diff        = sd(diff) / sqrt(n()),
              # Paired t-test: is knowing the true branch lengths worth anything?
              t_stat         = if (n() > 1 && sd(diff) > 0) t.test(true_bl, unit_bl, paired = TRUE)$statistic else NA_real_,
              p_value        = if (n() > 1 && sd(diff) > 0) t.test(true_bl, unit_bl, paired = TRUE)$p.value   else NA_real_,
              ci_lo          = if (n() > 1 && sd(diff) > 0) t.test(true_bl, unit_bl, paired = TRUE)$conf.int[1] else NA_real_,
              ci_hi          = if (n() > 1 && sd(diff) > 0) t.test(true_bl, unit_bl, paired = TRUE)$conf.int[2] else NA_real_,
              n_better       = sum(diff > 0), n_worse = sum(diff < 0), n_tied = sum(diff == 0),
              .groups = 'drop')

write.csv(rollup, file.path(SUMDIR, '2_branch_lengths_summary.csv'), row.names = FALSE)

sens <- acc %>% filter(metric_type == 'Sensitivity') %>%
    group_by(matrix_type, arm, class) %>%
    summarise(mean_recall = mean(accuracy, na.rm = TRUE),
              mean_specificity = mean(specificity, na.rm = TRUE), .groups = 'drop')
write.csv(sens, file.path(SUMDIR, '2_branch_lengths_per_class.csv'), row.names = FALSE)

conf <- read_all('_confusion\\.csv$')
if (!is.null(conf)) {
    conf_sum <- conf %>% group_by(matrix_type, arm, truth, selected) %>%
        summarise(n = sum(Freq), .groups = 'drop')
    write.csv(conf_sum, file.path(SUMDIR, '2_branch_lengths_confusion.csv'), row.names = FALSE)
}

cat('\n== paired rollup (true_bl - unit_bl) ==\n')
print(as.data.frame(rollup), row.names = FALSE, digits = 3)
cat(sprintf('\nSummaries written to %s\n', SUMDIR))
