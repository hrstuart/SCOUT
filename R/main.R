#' Wrapper function to run SCOUT 
#' @param counts.file path to a counts file with a columns for species and relevent columns for evolutionary models being tested if an OUX. Other columns will be genes. 
#' @param tree.file path to newick file. tips match species names in counts.
#' @param results_dir output directory
#' @param regimes a list of regimes to test. This include specifying BM1, OU1, and OUx where OUx correspond to columns in the counts file.
#' @param species_key if the species column is not called 'species' (rather cellBC for example), indicate the correct name here. 
#' @param blacklist genes for testing are inferred automatically. if any columns in the counts file are not for testing and not either the species key or model name, list them here.
#' @param method one of EM (expectation-maximization), SM (smoothing), MTF (tip fog).
#' @import paleotree
#' @import corpcor
#' @import nloptr
#' @import future.apply 
#' @import ape
#' @import phylolm
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
	tau_prior_sd = 0.1,
	tau_prior_mean=0.2,
	fixed_root = FALSE, 
	lambda1 = 0.2, 
	lambda2 = 0.2, 
	cores = 1, 
	logfile = NULL, 
	verbose = TRUE
	){

	log_message(sprintf('Started logging @ %s', logfile), verbose = verbose)
	log_message(sprintf('Data will be saved to --> %s', results_dir), logfile, verbose=verbose)
	log_message('=========================================================\n', logfile, verbose=verbose)

    counts <- read.csv(counts.file,  row.names=1)
    tree <- ape::read.tree(tree.file)

    # Overcomplicating this but basically, if its 0 or not in the samplesheet then we skip it (make NULL). Otherwise use value. 
    if (!is.null(smoothing_k)) {
      if (ss[i, 'smoothing_k'] == 0){
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

    idata <- formatSCOUT(tree_path = tree, 
      metadata_path = counts, 
      species_key = species_key, 
      anc_infer = 'ape', 
      outpath = results_dir, 
      regimes = regimes,
      normalize = normalize,
      smoothing_k = ska, ## this should be null 
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


	annotated <- annotate_history(history, datasetid = 'dataset') %>% 
	    arrange(AICc) %>% mutate(AICc_next_worse = lead(AICc) - AICc) %>% 
	    arrange(AIC) %>%  mutate(AIC_next_worse = lead(AIC)-AIC) %>% ungroup() %>% 
	                      arrange(gene_name)
	write.csv(annotated, sprintf('%s/%s_all_genes_full_history.csv', results_dir, tid))

	annotate_history_select <- annotated %>% filter(delta_AIC == 0)
	write.csv(annotate_history_select, sprintf('%s/%s_annotated_best_fit.csv', results_dir, tid))


  thetas <- list()
  params.full <- extract_parameters(full.res)
  for (r in idata$regimes){
    tmp <- parse_nested_to_df(params.full, r)
    thetas[[r]] <- tmp 
    write.table(tmp, sprintf('%s/%s_all_genes_%s_parameters.csv', results_dir, tid, r))
  }

	return(list(SCOUT_class = annotate_history_select, SCOUT_params = thetas))

}


