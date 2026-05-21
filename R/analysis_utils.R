extract_history_grid_search <- function(results_list, full.history=FALSE){
    dataset1 <- do.call(rbind, lapply(results_list, function(x){
        if (!full.history){
          history <- x$final 
        } else {
          history <- x$history 
        }
        
        setting <- x$settings
        history$settings <- paste(setting, collapse = '_')
        history$gene_name <- x$settings[[1]] #$gene_name
        history$regime <- x$settings[[2]]
        history$model <- x$defaults$model
        history$param.count <- length(unlist(x$paras))
        if (unique(history$model) == 'BM1'){
            history$param.count <- history$param.count - 1
        }
        return(history)
    }))
    return(dataset1)
}

extract_parameters <- function(results_list){
    per_model_list <- list()
    for (res in results_list){
        gene_name <- res$settings[[1]]
        reg <- res$settings[[2]]
        unq_regs <- res$defaults$parsed_alt_tree$unique_regimes
        
        if (reg == 'BM1'){
            unq_regs <- unq_regs[-1] # strip off root if BM1. 
        }
        tt <- unlist(res$paras$theta)
        names(tt) <- paste0('theta_', unq_regs) 
        res$paras <- c(res$paras, tt)
        res$paras[["theta"]] <- NULL
        if (reg %in% names(per_model_list)){
            per_model_list[[reg]][[gene_name]] <- res$paras
        } else {
            per_model_list[[reg]][[gene_name]] <- res$paras 
        }
    }
    return(per_model_list)
}


annotate_history <- function(dataset1, datasetid){
    grp1 <- c(datasetid, 'gene_name', 'regime')
    grp2 <- c(datasetid, 'gene_name')
    dataset2 <- dataset1 %>% group_by(!!!rlang::syms(grp1)) %>% 
        arrange(desc(iter)) %>% slice_head(n = 1) %>% ungroup() %>% 
        group_by(!!!rlang::syms(grp2))  %>% 
        mutate(ntips = 256) %>%
        mutate(AIC = -2*ll_total+2*param.count, 
               AICc = -2*ll_total+(2*param.count*(ntips/(ntips-param.count-1))), 
               delta_AIC = AIC - min(AIC),
               delta_AICc = AICc - min(AICc), 
               AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC)),
               AICc_weight = exp(-0.5 * delta_AICc) / sum(exp(-0.5 * delta_AICc))) %>% ungroup() %>% 
        mutate(truth = str_extract(gene_name, 'BM1|OU1|OUM'))

    return(dataset2)
}


library(caret)
library(dplyr)

calculate_group_class_accuracy <- function(data, dataset_col = "dataset", 
                                           model_col = "model", truth_col = "truth") {
  
  # Convert to factors to ensure consistent levels
  levs <- unique(c(data[[model_col]], data[[truth_col]]))
  data <- data %>%
    mutate(
      pred = factor(!!sym(model_col), levels = levs),
      ref = factor(!!sym(truth_col), levels = levs)
    )
  
  # Get unique classes
  classes <- unique(c(levels(data$pred), levels(data$ref)))
  
  # Calculate per-group overall accuracy with 95% CI
  overall_accuracy <- data %>%
    group_by(!!sym(dataset_col)) %>%
    summarise(
      n = n(),
      n_correct = sum(pred == ref),
      accuracy = n_correct / n,
      ci_test = list(binom.test(n_correct, n, conf.level = 0.95)),
      ci_lower = ci_test[[1]]$conf.int[1],
      ci_upper = ci_test[[1]]$conf.int[2],
      ci_test = NULL,
      metric_type = "Overall",
      class = NA_character_,
      .groups = 'drop'
    ) %>%
    select(!!sym(dataset_col), metric_type, class, accuracy, n, n_correct, ci_lower, ci_upper)
  
  # Calculate per-class accuracy (Sensitivity/Recall for each class)
  per_class_accuracy <- data %>%
    group_by(!!sym(dataset_col)) %>%
    reframe(
      class = classes,
      sensitivity = sapply(classes, function(c) {
        mask <- ref == c
        if (sum(mask) == 0) return(NA)
        sum(pred[mask] == c) / sum(mask)
      }),
      specificity = sapply(classes, function(c) {
        mask <- ref != c
        if (sum(mask) == 0) return(NA)
        sum(pred[mask] != c) / sum(mask)
      }),
      .groups = 'drop'
    ) %>%
    rename(class = class, accuracy = sensitivity) %>%
    mutate(metric_type = "Sensitivity", n = NA, n_correct = NA, 
           ci_lower = NA, ci_upper = NA) %>%
    select(!!sym(dataset_col), metric_type, class, accuracy, n, n_correct, ci_lower, ci_upper, specificity)
  
  # Combine results
  combined <- bind_rows(
    overall_accuracy %>% mutate(specificity = NA),
    per_class_accuracy
  )
  
  return(combined)
}
