#' Log messages to a file and print optionally
#' @param message Text to write
#' @param log_file Path to log file
#' @param verbose Logical, whether to print message
#' @export
log_message <- function(message, log_file = NULL, verbose = FALSE) {
  timestamped_msg <- paste0(Sys.time(), " | ", message)
  
  if (!is.null(log_file)) {
    cat(timestamped_msg, file = log_file, append = TRUE, sep = "\n")
  }
  
  if (verbose) {
    message(timestamped_msg)
  }
}

#' Run OUwie model on a single gene and regime
#' @export
run_one_ouwies_one_gene <- function(phy, meta, regime, gene, model, algorithm, mode, sh) {
    log_message(sprintf('Running... %s ~ %s, Scale? %s', model, gene, sh), verbose=FALSE)
    if (mode != 'dichotomous') {
        if (model == 'BM1' ){
            starting_vals = c(0.1)
        } else {
            starting_vals = rep(0.1, length(unique(indat1[[regime]]))+1)
        }
    } else {
        starting_vals = NULL
    }
    
    indat1 <- meta[,c('species', regime, gene)]
    res1 <- tryCatch(
          {
            OUwie(phy, indat1, 
                 model = model, 
                 algorithm=algorithm,
                 starting.vals = starting_vals,
                 scaleHeight = sh, 
                 ub = c(1e6, 1e6, 1e6), # hard coded, not sure if I want to leave it like that but okay for now. 
                 quiet=TRUE,
                 warn=FALSE)
          }, 
          error = function(e){
            log_message(paste('ERROR - ', model, 'failed for', gene, ' with error: ', e$message), verbose=TRUE)
            return(NULL)
          } ) 

    return(res1)}

#' Process OUwie model results into data frame summarizing parameters and model fit
#' @export
process_ouwie_single <- function(res, gene, regime) {
    selected<-res[c('loglik','AIC','AICc','BIC', 'model','param.count')]
    selected_df <- as.data.frame(selected)
    selected_df$quant_trait <- gene
    selected_df$regime <- regime

    states = names(res$solution['theta', ])
    sln = c(res$solution['alpha', 1], res$solution['sigma.sq', 1], res$solution['theta', ])
    if (is.null(states)){ states = 'theta'}
    names(sln) <- c('alpha', 'sigma.sq', states)

    selected_df <- cbind(selected_df, t(as.data.frame(sln)))
    row.names(selected_df) <- NULL

    return(selected_df)}

#' Infer ancestral state using ape or castor
#' @export
infer_anc <- function(phy, mode='ape') {
    if (mode == 'ape'){
        anc_res <- ape::ace(phy$states, 
        phy, 
        type = "discrete", method = "ML", CI = TRUE, model =  "ER", scaled = TRUE, kappa = 1, corStruct = NULL, ip = 0.1,
        use.expm = FALSE, use.eigen = TRUE, marginal = FALSE)
        most_likely_states <- apply(anc_res$lik.anc, 1, function(x) colnames(anc_res$lik.anc)[which.max(x)])
        names(most_likely_states) <- NULL
    } else if (mode == 'castor') {
        states <- as.factor(phy$states)

        anc_res <- castor::asr_max_parsimony(phy, as.numeric(states))
        most_likely_states <- levels(states)[anc_res$ancestral_states]

        if (anc_res$success == FALSE){
            stop('Unable to estimate ancestral states. Please either check inputs or use `ape` mode.')
        } 
    } 
    
    phy$node.label <- most_likely_states

    return(phy)
}

#' Normalize a numeric vector to 0-1
#' @export
normalize <- function(D) {
  return((D - min(D)) / (max(D) - min(D)))
}

#' Smoothing function for gene expression
#' @export
LT.smooth.all <- function(tree_, lognorm__, ka, s=NULL, v = 'v3', norm = 'B') {
    adj_mat <- as.matrix(cophenetic(tree_), sparse = TRUE)
    
    tip_order <- row.names(adj_mat)
    ln_ordered <- as.matrix(lognorm__[tip_order, ])
    
    m <- nrow(adj_mat)
    sigmas <- rep(0, m)
    windows <- matrix(nrow = nrow(adj_mat), ncol = ncol(adj_mat))
    
    if (v == 'v1') {
        for (i in 1:m){
            dists <- adj_mat[i, ]
            max_dist <- sort(dists)[ka]
            windows[i,] <- dists >= max_dist # store the boolean of whether the distance is greater or less than the max dist. TRUE = make 0. 
            sigmas[i] <- max_dist
        }
    } else if (v == 'v2'){
        for (i in 1:m){
            dists <- adj_mat[i, ]
            max_dist <- sort(dists)[ka]
            windows[i,] <- dists > max_dist # store the boolean of whether the distance is greater or less than the max dist. TRUE = make 0. 
            sigmas[i] <- max_dist
        }
    } else if (v == 'v3') {
        for (i in 1:m){
            dists <- adj_mat[i, ]
            max_dist <- sort(dists)[ka+1]
            windows[i,] <- dists > max_dist # store the boolean of whether the distance is greater or less than the max dist. TRUE = make 0. 
            sigmas[i] <- max_dist
        }
    }
    
    # Overwrite sigmas if there is an established s, otherwise do the dynamic option. 
    if (!is.null(s)){
        sigmas <- rep(s, m)}
    
    sigmas[sigmas == 0] <- 1 # replace zeros with 1 to remove division by zero errors, not sure if this is the way to go though. 

    affinity_matrix <- exp(- (adj_mat^2 / ( sigmas) ))
    affinity_matrix <- affinity_matrix + t(affinity_matrix) # make symmetric
    affinity_matrix[windows] <- 0
    
    if (norm == 'A'){
        normalized_affinity <- as.matrix(normalize(affinity_matrix))
    } else if (norm == 'B') {
        normalized_affinity <- as.matrix(affinity_matrix / rowSums(affinity_matrix))
        diag(normalized_affinity) <- 1 
    }
    
    smoothed_data <- normalized_affinity %*% ln_ordered
    
    if (sum(rownames(smoothed_data) == tree_$tip.label) != length(tip_order)) {
        stop('Data not in the same order at tree. Check results before moving on.')
    }

    return(smoothed_data)
}

#' Add a new column based on condition
#' @export
add_new_column <- function(df, col1_name, col2_name, new_col_name, other_list) {
  df[[new_col_name]] <- ifelse(df[[col1_name]] %in% other_list, df[[col2_name]], df[[col1_name]])
  return(df)
}

#' Wrapper function to return the model results
#' @export
return_model <- function(model_input, meta, regime, candidate, options) {
  res_model <- run_one_ouwies_one_gene(model_input[[regime]]$tree, meta, model_input[[regime]]$rcol, candidate, model_input[[regime]]$model, 'three.point', options$resolve, options$scale)
    if (!is.null(res_model)) {
      res <- process_ouwie_single(res_model, candidate, regime)
    } else {
      res <- NULL
    }
  return(res)
}
