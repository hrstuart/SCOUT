#' Prepare data and inputs for OUwie modeling
#' @param tree_path Path to the phylogenetic tree file in Newick format
#' @param metadata_path Path to the metadata CSV file
#' @param outpath Output directory path
#' @param species_key Column name for species ID in metadata (optional)
#' @param quant_traits Vector or file path of quantitative trait names or 'infer'
#' @param name Prefix for output files
#' @param regimes Vector of regime column names or single regime name
#' @param algorithm OUwie algorithm
#' @param anc_infer Ancestral state inference method: 'ape', 'castor', or 'skip'
#' @param resolve Tree multifurcation resolution method
#' @param normalize Logical, whether to log-normalize data
#' @param smoothing Optional integer smoothing parameter
#' @param floor Numeric floor for zero counts
#' @param scale T/F whether to scale tree height to 1. 
#' @param blacklist Vector of column names to exclude when inferring quantitative traits
#' @param troubleshoot Logical, whether to limit genes for troubleshooting
#' @param write_smooth Logical, whether to save smoothed counts intermediate file
#' @return List with processed inputs, metadata, gene columns, regimes, and attributes
#' @export
prepare_data <- function(tree_path, metadata_path, outpath, species_key = NULL, quant_traits = NULL,
                        name = NULL, regimes = "BM1", algorithm = "three.point", anc_infer = "ape",
                        resolve = "dichotomous", normalize = FALSE, smoothing = NULL, floor = 0, scale = FALSE,
                        blacklist = NULL, troubleshoot = FALSE, write_smooth = FALSE, logfile=NULL) {

    create_directory_if_not_exists(outpath)
    
    if (anc_infer == 'ape' & resolve != 'dichotomous'){
        stop('Mode must be dichotomous if using ape to infer ancestral states.')
    } 

    if (is.null(name)){
        prefix <- gsub('.*$', '', tree_path)
    } else {
        prefix <- name
    }

    if (class(metadata_path) == 'data.frame'){
        meta <- metadata_path
    } else {
        meta = read.csv(metadata_path, row.names=1)
    }

    if (!('species' %in% colnames(meta))){
        if (!is.null(species_key)){
            names(meta)[names(meta) == species_key] <- 'species'
        } else {
            stop('Please provide a column named `species` or use the --species_key option to identify the `species` column')
        }
    }

    reg_oums <- regimes[! regimes %in% c('BM1', 'OU1')]
    if (length(intersect(regimes, colnames(meta))) != length(reg_oums)){
        stop('Specified regimes not found in metadata columns: %s. Please correct. Exiting...', paste(reg_oums, collapse=', '))
    }

    if (is.list(quant_traits) || is.vector(quant_traits)){
        gene_cols <- unlist(quant_traits)
    } else if (is.character(quant_traits)) {
        if (file.exists(quant_traits)){
            gene_cols <- read.table(quant_traits)[[1]]
        }
    } else {
        log_message('Inferring gene list from input data columns.', log_file=logfile,verbose=TRUE)
        all_cols <- colnames(meta)
        gene_cols <- all_cols[! all_cols %in% c('species', regimes, blacklist)]
    }

    # Clean up gene_cols to match R formatted column names 
    gene_cols <- gsub('-', '.' , gene_cols) 

    gene_cols <- intersect(gene_cols, colnames(meta))
    ngenes <- length(gene_cols)
    log_message(paste0(ngenes, ' genes were detected. '), log_file=logfile, verbose=TRUE)
    meta[, gene_cols] <- apply(meta[, gene_cols], 2, as.numeric)
    rownames(meta) <- meta$species

    # Check for all zeros columns 
    meta <- as.data.frame(meta)
    sum_zero_genes <- gene_cols[colSums(meta[, gene_cols]) == 0]
    if (length(sum_zero_genes) > 0) {
        log_message(paste0('Removing ', length(sum_zero_genes), ' from the dataset which have column sumns of zero.'), log_file=logfile, verbose=TRUE)
        gene_cols <- gene_cols[colSums(meta[, gene_cols]) != 0]
    }

    ################### Tree Preprocessing ###################
    # intersect leaves with cells in dataframe and recommend prunning the tree accordingly. 
    if (class(tree_path) == 'phylo'){
        rawtree <- tree_path 
    } else {
        rawtree = ape::read.tree(tree_path)
    }
    
    log_message(paste0(length(rawtree$tip.label), ' leaves were found in the raw tree. ', length(intersect(rawtree$tip.label, meta$species)), ' leaves intersect with species column.'),log_file=logfile,  verbose=TRUE)

    # Check that all cells in a tree are included in the metadata. Force a stop if there is missing data.
    ## This could likely be improved upon later. 

    if (length(rawtree$tip.label) != length(intersect(rawtree$tip.label, meta$species))) {
        missing_cells = rawtree$tip.label[!rawtree$tip.label %in% meta$species]
        print(missing_cells)
        stop(sprintf('Missing metadata for %d leaves in the tree. Please prune the tree or impute relevant metadata', length(missing_cells)))
    } else {
      phy <- rawtree
    }

    log_message('Preprocessing tree',  verbose=TRUE)

    # If the ancestral state annotation is done using APE then the tree must be fully dichotomous. 
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode

    if (anc_infer == 'ape'){
        if (nb.node != nb.tip - 1){
            phy <- multi2di(phy)
        }
    }

    # Check whether there are edge lengths or not. For now, adding everyone in as a 1. 
    ## Improve strategy later. 
    if (is.null(phy$edge.length) ) {
        phy$edge.length <- rep(1, Nedge(phy))
        log_message('No edge lengths, replacing with 1s.', verbose=TRUE)
    } else if (sum(phy$edge.length == 0) > 0){
        log_message('Found 0 length edges. Replacing 0 edge lengths with with 1e-7.', verbose=TRUE)
        phy$edge.length[phy$edge.length == 0] <- 1e-7
    } else {
        log_message('Found edge lengths.', verbose=TRUE)
        #print(fivenum(phy$edge.length))
    }

    # Options for resolving multifurcations (either dichotomous -- above, or bifurcating here.)
    if (resolve == 'bifurcating'){
        phy <- castor::multifurcations_to_bifurcations(phy, dummy_edge_length=0.01) # hard coding this for now. 
    } else if (resolve == 'multifurcating' & nb.node != nb.tip - 1){
        log_message('Leaving multifurcations unresolved.', verbose=TRUE)
    }

    ################### Data Preprocessing ###################

    if (normalize){
        log_message(paste(Sys.time(), 'Log normalizing gene expression counts.'), log_file=logfile,  verbose=TRUE)
        meta[, gene_cols] <- log1p(meta[, gene_cols])
    }

    # Could probably improve upon the distance metric used in the smoothing but alas. 
    if (!is.null(smoothing)) {
        log_message(paste('Smoothing counts at k =', smoothing), log_file=logfile,  verbose=TRUE)

        smooth_meta <- LT.smooth.all(phy, meta[, gene_cols], smoothing)
        g_meta <- smooth_meta[phy$tip.label, ]

    } else {
        g_meta <- meta[phy$tip.label, gene_cols]
    }

    if (floor > 0){
        g_meta[g_meta == 0] <- floor
    }

    cln_meta <- meta[phy$tip.label, ]
    cln_meta[gene_cols] <- g_meta[phy$tip.label, ]

    cln_meta$OU1 <- 'global'#, length(phy$tip.label))

    if (write_smooth & !is.null(smoothing)){
        tmp_fle = paste0(outpath, '/', prefix, '_', resolve, '_smoothed_k', smoothing, '_counts.csv')
        log_message(paste('Writing intermediate smoothed counts to', tmp_fle), log_file=logfile,  verbose=TRUE)
        write.csv(cln_meta, tmp_fle)
    }

    ###################### Prep Inputs #######################

    model_input <- list()
    regime_edits <- c()

    # Kind of convoluted but if there are multiple OUM regimes, need to store a tree for each potential regime. 
    for (r in 1:length(regimes)){
        rname <- regimes[r]
        if (length(unique(cln_meta[[rname]])) > 1) {
            model <- 'OUM'
            if (anc_infer == 'skip'){
                phy$node.label <- rep('ANC', phy$Nnode)
                phy1 <- phy
            } else {
                phy$states <- cln_meta[phy$tip.label, rname]

                ### default is using APE but if you want to get around the fully dichotomous issue, then use castor. 
                phy1 <- infer_anc(phy, anc_infer) 
            }
            rcol <- rname
        } else {
            phy1 <- phy
            ### DONT CHANGE THING FROM ANC OUwie will fail without any regime shifts for some reason.
            phy1$node.label <- rep('ANC', phy$Nnode)
            rcol <- 'OU1'
            model <- rname 
        }
        log_message(paste0('Regime: ', rname, ' Rcol: ', rcol, ' Model: ', model), verbose = TRUE)
        model_input[[rname]] <- list(model = model, tree = phy1, rcol = rcol)
        phy$node.label <- NULL # reset labels 
    }

    rownames(cln_meta) <- NULL

    return(list('inputs'=model_input, 'meta_data' = cln_meta, 'gene_cols' = gene_cols, 'regimes' = regimes, 'attributes' = list('anc_infer' = anc_infer, 'resolve' = resolve, 'scale'=scale)))
}

#' Fit OUwie models in parallel
#' @param inputs List generated from prepare_data()
#' @param cores Number of cores for parallel processing
#' @param write Logical, whether to write per-regime CSV results
#' @param outpath Directory for output files
#' @param prefix File prefix
#' @param testgenes If genes of interest is a subset of the full gene columns list, provide a list of genes to test.
#' @param testregimes If regimes of interest is a subset of the full regime list, provide a list of regimes to test. 
#' @return Data frame of combined model summaries
#' @export
fitModel <- function(inputs, cores = 1, write = TRUE, outpath = NULL, prefix = "OUWIE", logfile=NULL,
    testgenes = NULL, testregimes = NULL) {

    totalCores = parallel::detectCores(logical = FALSE)
    if (cores >= totalCores) {
      stop(paste0('ERROR: too many cores requested. Please specify a number less than ', totalCores))
    }

    if (!is.null(testgenes)){
        genes <- intersect(inputs[['gene_cols']], testgenes)
    } else {
        genes <- inputs[['gene_cols']]
    }

    if (!is.null(testregimes)){
        regimes <- intersect(inputs[['regimes']], testregimes)
    } else {
        regimes <- inputs[['regimes']]
    }

    nReg <- length(regimes)
    nGenes <- length(genes)

    # Create all combinations of regime and gene indices
    combo_df <- expand.grid(i = seq_len(nReg), j = seq_len(nGenes))

    log_message(sprintf('Fitting %d model-gene combinations using %d cores.', nrow(combo_df), cores), verbose=TRUE)
    
    plan(multisession, workers = cores)

    result_list <- future_lapply(seq_len(nrow(combo_df)), function(idx){

        i <- combo_df$i[idx]
        j <- combo_df$j[idx]
        r <- regimes[i]
        g <- genes[j]

        res <- tryCatch({
            return_model(inputs[['inputs']], inputs[['meta_data']], r, g, inputs[['attributes']])
        }, error = function(e){
            log_message(paste0('Error in task regime:', r, ' gene:', g, ' - ', conditionMessage(e), '\n'), log_file=logfile, verbose=TRUE)
            return(NULL)
        })

        return(res)
    }, future.seed = TRUE)
    
    log_message('Done with model fitting.', verbose=TRUE)

    # Gather result 
    models_summary <- do.call(rbind, lapply(result_list, function(x) {x[, c('quant_trait', 'loglik', 'AIC', 'AICc', 'BIC', 'model', 'regime', 'param.count', 'alpha', 'sigma.sq')]}))

    regime_specific_res = list()

    if (write){
      if (is.null(outpath)){
          outpath <- getwd()
      }
      # Write per regime results
      for (r in regimes) {
        # Filter results for this regime
        res_subset <- do.call(rbind, lapply(result_list, function(x) {if (x$regime == r) {return(x)} else {return(NULL)} }))
        regime_specific_res[[r]] = res_subset
        write.csv(res_subset, file = paste0(outpath, '/', prefix, '_', r, '_model_results.csv'), row.names = FALSE)
      }
    log_message(sprintf('Writing files to %s.', outpath), verbose=TRUE)
    write.csv(models_summary, paste0(outpath, '/', prefix, '_model_results.csv'))
    }


    return((list('results' = models_summary, 'per_model' = regime_specific_res)))
}


