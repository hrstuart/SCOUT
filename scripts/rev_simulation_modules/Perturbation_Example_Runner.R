.libPaths(c('/dartfs/rc/lab/M/McKennaLab/projects/hannah/software/R/R-4.4.2/library'))

suppressPackageStartupMessages({
  library(SCOUT)
  library(progressr)
  library(future.apply)
  library(corpcor)
  library(paleotree)
  library(nloptr)
  library(dplyr)
  library(stringr)
})


ROOT     <- '/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/1_tree_robustness'
REPS     <- 20
CORES    <- 48
REGIMES  <- c('BM1', 'OU1', 'OU2', 'OUM', 'OU4')   # candidate partitions = counts-CSV columns # OUM = OU3 here. 
TOPRES_DIR  <- file.path(ROOT, 'output', '4_bl_reg_perturb_seeds_v2')  # SCOUT output
DAT_DIR  <- file.path(ROOT, 'data', '4_bl_reg_perturb_seeds_v2')
REPFILE  <- file.path(DAT_DIR, 'replicates_manifest.csv')
ARMS     <- c('evf', 'counts')

# The two expression arms, each fit independently on the same regime draw.
say <- function(...) {
  cat(sprintf('[%s] ', format(Sys.time(), '%Y-%m-%d %H:%M:%S')), ..., '\n', sep = '')
  flush.console()
}

annotate_history_fixed <- function(dataset1, datasetid) {
  grp1 <- c(datasetid, "gene_name", "regime")
  grp2 <- c(datasetid, "gene_name")
  dataset1 %>%
    group_by(!!!rlang::syms(grp1)) %>% arrange(desc(iter)) %>% slice_head(n = 1) %>% ungroup() %>%
    group_by(!!!rlang::syms(grp2)) %>%
    mutate(AIC  = -2 * ll_total + 2 * param.count,
           AICc = -2 * ll_total + (2 * param.count * (ntips / (ntips - param.count - 1))),
           delta_AIC  = AIC - min(AIC),
           delta_AICc = AICc - min(AICc),
           AIC_weight  = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC)),
           AICc_weight = exp(-0.5 * delta_AICc) / sum(exp(-0.5 * delta_AICc)),
           AICc_next_worse = lead(AICc) - AICc,
           AIC_next_worse  = lead(AIC) - AIC,
           best_fit = ifelse(delta_AIC == 0, TRUE, FALSE)) %>%
    ungroup() %>%                      # delta_AIC / delta_AICc deliberately retained
    mutate(truth = str_extract(gene_name, 'BM1|OU1|OUM|OU\\d+'))
}
assignInNamespace("annotate_history", annotate_history_fixed, ns = "SCOUT")


# ----------------------------------------------------------------------------
# Main loop
# ----------------------------------------------------------------------------

say('regime-perturbation seed replicates')
say('  regimes            : ', paste(REGIMES, collapse = ', '), '   (BM1 / OU1 not fitted)')
say('  replicate manifest : ', REPFILE)
say('  results            : ', TOPRES_DIR)
say('  cores              : ', CORES)

replicates <- read.csv(REPFILE)

dir.create(TOPRES_DIR)

manifest <- list()
failed   <- character(0)

for (i in 1:nrow(replicates)) {
  runinfo  <- replicates[i, ]
  TREEFILE <- runinfo$tree 
  RES_DIR  <- file.path(TOPRES_DIR, runinfo$subset)
  COUNTS   <- runinfo$counts 
  EVF      <- runinfo$evf
  REP      <- runinfo$rep
  REPID    <- gsub('_OUEVF.*', '', basename(COUNTS))

  for (arm in ARMS){
      if (arm == 'evf'){
        repfile <- EVF
      } else {
        repfile <- COUNTS
      }

      id   <- sprintf('%s_regperturb_%s', REPID, arm)
      done <- file.path(RES_DIR, paste0(id, '_all_genes_full_history.csv'))

      if (file.exists(done)) {
        say('  [', arm, '] already done, skipping (', basename(done), ')')
        next
      }

      say('  [', arm, '] SCOUT -> ', id)
      t0 <- Sys.time()
      ok <- tryCatch({
        SCOUT(counts.file = repfile,
              testid      = id,
              tree.file   = TREEFILE,
              results_dir = RES_DIR,
              regimes     = REGIMES,
              scale_tree  = FALSE,
              cores       = CORES,
              normalize   = if (arm == 'evf') FALSE else TRUE,  
              logfile     = file.path(RES_DIR, paste0(id, '.log')),
              verbose     = TRUE)
        TRUE
      }, error = function(e) {
        say('  [', arm, '] FAILED: ', conditionMessage(e))
        FALSE
      })

      if (!ok) {
        failed <- c(failed, id)
      } else {
        say('  [', arm, '] done in ',
            round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 1), ' min')
      }
  }
    
}


