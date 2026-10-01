#' Wrapper function to run SCOUT 
#' @param counts.file path to a counts file with a columns for species and relevent columns for evolutionary models being tested if an OUX. Other columns will be genes. 
#' @param tree.file path to newick file. tips match species names in counts.
#' @param results_dir output directory
#' @param regimes a list of regimes to test. This include specifying BM1, OU1, and OUx where OUx correspond to columns in the counts file.
#' @param species_key if the species column is not called 'species' (rather cellBC for example), indicate the correct name here. 
#' @param blacklist genes for testing are inferred automatically. if any columns in the counts file are not for testing and not either the species key or model name, list them here.
#' @param method one of EM (expectation-maximization), SM (smoothing), MTF (tip fog).
#' @param genes optional character vector of genes to fit (e.g. from \code{select_lambda_genes()}).
#'   Default NULL fits every gene column, as before.
#' @param lambda_filter if TRUE, run \code{lambda_screen()} on the same inputs before fitting and
#'   fit only genes passing \code{lambda_min} / \code{lambda_p_max}; the screen is written to
#'   \code{<testid>_lambda_screen.csv}. Default FALSE bypasses the screen (all genes are fit).
#' @param lambda_min,lambda_p_max cutoffs for the lambda filter (lambda >= lambda_min,
#'   LRT P <= lambda_p_max). At least one is required when \code{lambda_filter = TRUE}.
#' @param lambda_engine engine for \code{lambda_screen()}: 'auto' uses the fast engine on
#'   ultrametric trees and reports + falls back to the slow phytools engine otherwise.
#' @param hybrid if TRUE and both BM1 and OU1 are among \code{regimes}, also run
#'   \code{run_hybrid_pipeline()} and write \code{<testid>_hybrid_calls.csv}. The AIC-based
#'   outputs are unchanged.
#' @param hybrid_ic,hybrid_alpha information criterion and significance level for the hybrid
#'   selection.
#' @import paleotree
#' @import corpcor
#' @import nloptr
#' @import future.apply 
#' @import ape
#' @import phylolm
#' @import dplyr stringr
#' @export
SCOUT <- function(counts.file, tree.file, results_dir, 
	regimes, 
	species_key = 'species',
	blacklist = NULL, 
	method = 'EM',
	testid = NULL, 
	normalize = TRUE, 
	scale_tree = FALSE, 
	smoothing_k = NULL, 
	infer_anc = 'ape', 
	tau_prior_sd = 0.1,
	tau_prior_mean=0.2,
	fixed_root = FALSE, 
	lambda1 = 0.2, 
	lambda2 = 0.2, 
	cores = 1, 
	logfile = NULL, 
	verbose = TRUE,
	genes = NULL,
	lambda_filter = FALSE,
	lambda_min = NULL,
	lambda_p_max = NULL,
	lambda_engine = 'auto',
	hybrid = TRUE,
	hybrid_ic = 'AICc',
	hybrid_alpha = 0.05
	){

  create_directory_if_not_exists(results_dir)

	log_message(sprintf('Started logging @ %s', logfile), verbose = verbose)
	log_message(sprintf('Data will be saved to --> %s', results_dir), logfile, verbose=verbose)
	log_message('=========================================================\n', logfile, verbose=verbose)

    counts <- read.csv(counts.file,  row.names=1)
    tree <- ape::read.tree(tree.file)

    # Overcomplicating this but basically, if its 0 or not in the samplesheet then we skip it (make NULL). Otherwise use value. 
    if (!is.null(smoothing_k)) {
      if (smoothing_k == 0){
        ska <- NULL
      } else {
        ska <- smoothing_k
      } 
    } else {
      ska <- smoothing_k
    }
    # defaults just to make things easy for now. 
    if (is.null(testid)){
    	tid <- paste0('SCOUT_', paste(sample(c(0:9, letters, LETTERS), 8, replace = TRUE), collapse = ""))
    } else {
    	tid <- testid
    }
    
    if (method == 'EM'){
    	skipE <- FALSE; ska <- NULL; skipT <- FALSE
    } else if (method == 'MTF') {
    	skipE <- TRUE; ska <- NULL; skipT <- FALSE 
    } else if (method == 'SM') {
    	if (is.null(ska)){
    		log_message('Method set to SM but no value set for smoothing parameter. Defaults to 8.', logfile, verbose = verbose)
    		ska <- 8
    	}
    	skipE <- TRUE; skipT <- TRUE
    }

    # Optional Pagel's lambda screen: drop genes with no phylogenetic signal before fitting.
    lambda_res <- NULL
    if (lambda_filter) {
      if (is.null(lambda_min) && is.null(lambda_p_max)) {
        stop('lambda_filter = TRUE needs lambda_min and/or lambda_p_max (see plot_lambda_screen()).')
      }
      lambda_res <- lambda_screen(counts, tree, species_key = species_key, regimes = regimes,
        blacklist = blacklist, genes = genes, normalize = normalize, engine = lambda_engine,
        cores = cores, lambda_min = lambda_min, p_max = lambda_p_max, logfile = logfile,
        verbose = verbose)
      write.csv(as.data.frame(lambda_res), sprintf('%s/%s_lambda_screen.csv', results_dir, tid),
        row.names = FALSE)
      genes <- lambda_res$gene_name[lambda_res$pass]
      log_message(sprintf('Lambda filter: keeping %d of %d genes (lambda_min = %s, lambda_p_max = %s).',
        length(genes), nrow(lambda_res), format(lambda_min), format(lambda_p_max)), logfile, verbose = TRUE)
      if (!length(genes)) stop('No genes pass the lambda filter; nothing to fit.')
    }

    idata <- formatSCOUT(tree_path = tree, 
      metadata_path = counts, 
      quant_traits = genes, 
      species_key = species_key, 
      anc_infer = infer_anc, 
      outpath = results_dir, 
      regimes = regimes,
      normalize = normalize,
      smoothing_k = ska, 
      blacklist = blacklist, 
      logfile=logfile)

    full.res <- runSCOUT(idata, 
      lambda1=lambda1, 
      lambda2=lambda2, 
      fixed.root=fixed_root, 
      M_only=skipE, 
      runid=tid, 
      skipTau = skipT,
      tau_prior_mean = tau_prior_mean,
      tau_prior_sd = tau_prior_sd,
      scaleHeight = scale_tree, 
      cores = cores, 
      logfile = logfile 
      )

    filename <- sprintf('%s/%s_%s.rds', results_dir, tid, format(Sys.Date(), "%Y%m%d"))
    saveRDS(full.res, filename)
    log_message(sprintf('Done with initial analysis of %s.', tid), logfile, verbose=TRUE)
    log_message('=========================================================\n', logfile, verbose=TRUE)
    log_message(sprintf('Processsing Results... Saving files to %s', results_dir), logfile, verbose=TRUE)

    history <- extract_history_grid_search(full.res)
    if (!'converge' %in% colnames(history)) history$converge <- NA
    history$dataset <- tid

  history$ntips <- length(tree$tip.label)
	annotated <- annotate_history(history, datasetid = 'dataset') %>% 
    group_by(dataset, gene_name) %>% 
	    arrange(AICc) %>% 
      mutate(AICc_next_worse = lead(AICc) - AICc) %>% 
	    arrange(AIC) %>%  
      mutate(AIC_next_worse = lead(AIC)-AIC) %>% ungroup() %>% arrange(gene_name)
	write.csv(annotated, sprintf('%s/%s_all_genes_full_history.csv', results_dir, tid))

	annotate_history_select <- annotated %>% filter(delta_AIC == 0) %>%
    dplyr::select(iter, ll_total, sigma, alpha, tau, converge, gene_name, regime, model, param.count, ntips, AIC, AIC_weight)
  names(annotate_history_select)[which(names(annotate_history_select) == 'll_total')] <- 'loglik'
	write.csv(annotate_history_select, sprintf('%s/%s_annotated_best_fit.csv', results_dir, tid))

  # Hybrid selection: IC for BM1 vs OU1, Bonferroni LRT vs OU1 for multi-regime winners.
  hybrid_res <- NULL
  if (hybrid && all(c('BM1', 'OU1') %in% regimes)) {
    hybrid_res <- tryCatch(
      run_hybrid_pipeline(annotated, ic = hybrid_ic, alpha = hybrid_alpha),
      error = function(e) {
        log_message(sprintf('Hybrid selection failed and was skipped: %s', conditionMessage(e)), logfile, verbose = TRUE)
        NULL
      })
    if (!is.null(hybrid_res)) {
      write.csv(hybrid_res, sprintf('%s/%s_hybrid_calls.csv', results_dir, tid), row.names = FALSE)
    }
  } else if (hybrid) {
    log_message('Hybrid selection skipped: it needs both BM1 and OU1 in regimes.', logfile, verbose = verbose)
  }


  thetas <- list()
  params.full <- extract_parameters(full.res)
  for (r in idata$regimes){
    tmp <- parse_nested_to_df(params.full, r)
    thetas[[r]] <- tmp 
    write.table(tmp, sprintf('%s/%s_all_genes_%s_parameters.csv', results_dir, tid, r))
  }

	return(list(SCOUT_class = annotate_history_select, SCOUT_params = thetas, SCOUT_input = idata,
	  SCOUT_hybrid = hybrid_res, SCOUT_lambda = lambda_res))

}


