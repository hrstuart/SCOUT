# 05_status.r -- status tracker for the branch-length baseline run. Rebuilds the status CSV from
# files on disk (manifest, output CSVs, claim files), so it is safe to delete and re-run.
#
# Usage:
#   Rscript 05_status.r                                  # one refresh, print the rollup
#   SCOUT_TRACK_INTERVAL=60 Rscript 05_status.r loop      # refresh until told to stop

ROOT <- Sys.getenv('SCOUT_ROOT',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/1_tree_robustness')
source(file.path(ROOT, 'scripts', '1_tree_perturb_bl', 'bl_paths.R'))

P        <- bl_paths(ROOT)
args     <- commandArgs(trailingOnly = TRUE)
LOOP     <- length(args) > 0 && args[1] == 'loop'
# Names the status CSV, the stop file and the claim-file prefix. Must agree with the GOAL the
# runner passes to run_one.sh, or the tracker and the queue would read different files. Unset, this
# is the original run; the ultrametric-restored variant sets 1_tree_perturb_bl_um.
GOAL     <- Sys.getenv('SCOUT_GOAL', '1_tree_perturb_bl')
JOBS     <- Sys.getenv('SCOUT_JOBS', P$jobs_sweep)
LOGDIR   <- Sys.getenv('SCOUT_LOGDIR', P$logdir)
OUT_ROOT <- Sys.getenv('SCOUT_OUTROOT', P$outroot)
POOL     <- as.integer(Sys.getenv('SCOUT_POOL', '6'))
INTERVAL <- as.integer(Sys.getenv('SCOUT_TRACK_INTERVAL', '60'))
STOPFILE <- file.path(LOGDIR, sprintf('.track_stop_%s', GOAL))

STATUSDIR <- file.path(LOGDIR, 'status')
dir.create(STATUSDIR, recursive = TRUE, showWarnings = FALSE)
CSV    <- file.path(LOGDIR, sprintf('run_status_%s.csv', GOAL))
SUMCSV <- file.path(LOGDIR, sprintf('run_status_summary_%s.csv', GOAL))

if (!file.exists(JOBS)) stop(sprintf('jobs manifest %s not found -- run 02_build_jobs.r first.', JOBS))

read_claim <- function(f) {
    kv <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
    if (!length(kv)) return(NULL)
    parts <- strsplit(kv, '=', fixed = TRUE)
    parts <- parts[vapply(parts, length, 1L) >= 2]
    if (!length(parts)) return(NULL)
    setNames(vapply(parts, function(p) paste(p[-1], collapse = '='), ''),
             vapply(parts, `[`, '', 1))
}

write_atomic <- function(df, path) {
    tmp <- paste0(path, '.', Sys.getpid(), '.tmp')
    write.csv(df, tmp, row.names = FALSE)
    file.rename(tmp, path)
}

refresh <- function() {
    jobs <- read.csv(JOBS, stringsAsFactors = FALSE)
    dirs  <- task_dir_bl(jobs, OUT_ROOT)
    stems <- task_stems_bl(jobs)

    done <- vapply(seq_len(nrow(jobs)), function(i) task_complete(dirs[i], stems[i]), logical(1))

    status   <- ifelse(done, 'done', 'pending')
    started  <- rep(NA_character_, nrow(jobs))
    finished <- rep(NA_character_, nrow(jobs))
    exitcode <- rep(NA_character_, nrow(jobs))
    logfile  <- rep(NA_character_, nrow(jobs))

    claims <- list.files(STATUSDIR, pattern = sprintf('^%s_[0-9]+\\.claim$', GOAL), full.names = TRUE)
    for (f in claims) {
        cl <- read_claim(f)
        if (is.null(cl) || is.na(cl['task_id'])) next
        i <- as.integer(cl[['task_id']])
        if (is.na(i) || i < 1 || i > nrow(jobs)) next
        started[i]  <- if ('started'  %in% names(cl)) cl[['started']]  else NA_character_
        finished[i] <- if ('finished' %in% names(cl)) cl[['finished']] else NA_character_
        exitcode[i] <- if ('exit_code' %in% names(cl)) cl[['exit_code']] else NA_character_
        logfile[i]  <- if ('log_file' %in% names(cl)) cl[['log_file']] else NA_character_
        # A claim never overrides the on-disk verdict: `done` is decided by the output files alone.
        # The claim only distinguishes "not started" from "running" from "failed".
        if (!done[i]) {
            status[i] <- if (is.na(finished[i])) 'running'
                         else if (identical(exitcode[i], '0')) 'incomplete'   # exited 0, wrote nothing usable
                         else 'failed'
        }
    }

    out <- data.frame(task_id = jobs$task_id, layer = jobs$layer, alpha = jobs$alpha,
                      arm = jobs$arm, intensity = jobs$intensity, replicate = jobs$replicate,
                      status = status, started = started, finished = finished,
                      exit_code = exitcode, log_file = logfile, stringsAsFactors = FALSE)
    write_atomic(out, CSV)

    summ <- as.data.frame(table(layer = out$layer, alpha = out$alpha, arm = out$arm,
                                intensity = out$intensity, status = out$status))
    summ <- summ[summ$Freq > 0, ]
    write_atomic(summ, SUMCSV)

    n <- nrow(out); nd <- sum(out$status == 'done'); nr <- sum(out$status == 'running')
    nf <- sum(out$status %in% c('failed', 'incomplete'))

    # ETA from measured throughput, not from a nominal per-task time: mean elapsed over the tasks
    # that actually finished, divided by the pool width.
    eta <- ''
    accf <- list.files(OUT_ROOT, pattern = '_accuracy\\.csv$', recursive = TRUE, full.names = TRUE)
    accf <- accf[!grepl('/summaries/', accf)]
    if (length(accf) >= 3 && nd < n) {
        el <- suppressWarnings(unlist(lapply(head(accf, 400), function(f) {
            d <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
            if (is.null(d) || !'elapsed_min' %in% names(d)) NULL else d$elapsed_min[1]
        })))
        el <- el[is.finite(el)]
        if (length(el)) eta <- sprintf('  | mean %.1f min/task, ETA ~%.1f h',
                                       mean(el), (n - nd) * mean(el) / POOL / 60)
    }

    cat(sprintf('\n[%s] %s : %d/%d done (%.1f%%), %d running, %d failed%s\n',
                format(Sys.time(), '%Y-%m-%d %H:%M:%S'), GOAL, nd, n, 100 * nd / n, nr, nf, eta))

    tab <- table(paste(out$layer, out$alpha, out$arm, sprintf('%.2f', out$intensity)), out$status)
    print(tab)

    if (nf > 0) {
        cat('\nfailed / incomplete tasks:\n')
        bad <- out[out$status %in% c('failed', 'incomplete'), ]
        print(head(bad[, c('task_id', 'layer', 'alpha', 'arm', 'intensity', 'replicate',
                           'exit_code', 'log_file')], 20), row.names = FALSE)
    }
    invisible(nd == n)
}

if (!LOOP) {
    refresh()
} else {
    repeat {
        all_done <- refresh()
        if (file.exists(STOPFILE)) { unlink(STOPFILE); refresh(); break }
        if (all_done) break
        Sys.sleep(INTERVAL)
    }
}
cat(sprintf('\nstatus : %s\nsummary: %s\n', CSV, SUMCSV))
