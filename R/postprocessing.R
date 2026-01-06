#' Calculate fit metrics and weights from model summaries
#' @param models_summary Data frame returned from fitModel
#' @param outpath Directory to write weighted results CSV
#' @param prefix File prefix
#' @param write Logical, whether to write weighted results CSV
#' @return Data frame with AICc ranks, delta and weights
#' @export
get_fit_metrics <- function(results, outpath = NULL, prefix = "OUWIE", write = FALSE) {
    
    if (class(results) == 'list') {
        rdf <- results[['results']]
    } else {
        rdf <- results
    }
    results_df <- rdf %>% 
        group_by(quant_trait) %>% mutate(aic_rank = rank(AICc)) %>%  
        mutate(delta_aicc = AICc - min(AICc), 
           aicc_weight = round(exp(-0.5 * delta_aicc) / sum(exp(-0.5 * delta_aicc)), 6)) %>% 
        mutate(aiccw_rank = rank(-aicc_weight)) %>% ungroup() %>% data.frame()
    
    if (write){
        if (is.null(outpath)){
            outpath <- getwd()
        }
        write.csv(results_df, paste0(outpath, '/', prefix, '_model_results_weighted.csv'))
    }

    if (is.list(results)){
        results[['results']] <- results_df 
        return(results)
    } else {
        return(results_df)
    }
}

#' Process results to pick winners based on delta AICc thresholds
#' @export
processes_results <- function(i, d=2, ac=0.1) {
    if (is.character(i)) {
        if (file.exists(i)){
            res_x <- read.csv(i, row.names = 1) %>% distinct()
        }
    print(class(i))
    } else if (class(i) == 'list'){
        res_x <- i[['results']]
    } else if (class(i) == "data.frame"){
        res_x <- i
    } else {
        stop('Input does not match expected. Please check. Should be a list, a file path, or a data.frame.')
    }

    res_i <- tryCatch({
        res_x %>%
            group_by(quant_trait) %>%
            dplyr::arrange(desc(aicc_weight)) %>%
            slice_head(n = 2) %>%
            mutate(full_delta = sum(delta_aicc)) %>%
            filter(aic_rank == 1) %>%
            mutate(selection = ifelse(full_delta < d, "ambiguous", "confident")) %>%
            select(-c(loglik, AIC, AICc, BIC, model, aic_rank, delta_aicc, aiccw_rank)) %>%
            mutate(reliable_alpha = ifelse(regime == "BM1", "BM1", ifelse(alpha < ac, "unstable", "stable"))) %>%ungroup() 
        }, error = function(e) {
          stop(sprintf("Failed with message: %s", conditionMessage(e)))
      })
      

    tot_ambig <- nrow(res_i[which(res_i$selection == 'ambiguous'), ])

    if (tot_ambig > 0) {
        print(sprintf('%d genes have ambiguous evolutionary models.', tot_ambig))
    }
    
    if (nrow(res_i) != nrow(distinct(res_i))){
        stop('Error in post-processing. Models not selected correctly. Aborting.')
    }

    if (class(i) == 'list'){
        i[['model_selection']] <- res_i 
        i[['deltaAIC']] <- d 
        i[['alpha_cutoff']] <- ac

        return(i)
    }

    return(res_i)
}

#' Collapse results files into a single dataframe with parameter estimates. 
#' @export
combined_estimates <- function(x){
    df<- read_files_safe(x)
    if (is.null(df) ){
        return(NULL)
    }
    suppressMessages({
        collect <- melt(df[, c(7, 8, 11:ncol(df))], id_vars = c('quant_trait', 'regime'))
    })

    colnames(collect) <- c('quant_trait', 'regime', 'theta_param', 'theta')
    return(collect)
}

#' Read in files, remove ones that are empty
#' @export
read_files_safe <- function(file){
    tryCatch({
    # Check if file exists and is not empty
    if (file.info(file)$size < 10) {
      message("File is empty")
      return(NULL)
    }
        
    content <- read.csv(file)
    return(content)
    }, error = function(e) {
        message('Error reading file: ', e$message)
        return(NULL)
    
    })
}

#' Calculate accuracy metrics using confusion matrix
#' @export
calc_accuracy <- function(df) {
   a <- as.factor(df$regime)
   b <- as.factor(df$true_regime)

    res <- confusionMatrix(a, b)
    print(res$overall['Accuracy'])
    res$byClass

    return(res)
}

#' Create directory if it does not exist
#' @export
create_directory_if_not_exists <- function(directory_path, message = TRUE) {
  if (!dir.exists(directory_path)) {
    dir.create(directory_path, recursive = TRUE)
    if (message) cat("Directory created:", directory_path, "\n")
  } else {
    if (message) cat("Directory already exists:", directory_path, "\n")
  }
}
