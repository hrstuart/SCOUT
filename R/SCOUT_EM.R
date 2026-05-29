
####### EM #######

e_step <- function(X, phy, edges, paras, defaults, diagnose=TRUE){
    #V <- quickVCV(phy, paras$alpha, paras$sigma, defaults$scaleHeight)
    if (defaults$weight_type == 'alt'){
      V <- compute_VCV(defaults$parsed_alt_tree, paras$alpha, paras$sigma, add.root=defaults$add.root) # stop gap will fix. 
    } else {
      V <- varcov.ou(phy, edges, paras$alpha, paras$sigma, defaults$root.state, defaults$scaleHeight, defaults$assume.station)
    }
    
    Psi <- paras$tau^2
    if (paras$tau < 0){
      stop('Error, tau is < 0.')
    } 
    ridge <- 1e-8 * diag(nrow(V))
    V_S <- V + diag(Psi, nrow(V)) + ridge 

    L <- tryCatch(chol(V_S), error = function(err){
      if (defaults$verbose) print('Warning: tree likely bad. Cannot calculate inverse with cholesky decomp.')
      return(NULL)
    })
    if (is.null(L)){
      V_S_inv <- tryCatch(pseudoinverse(V_S), error = function(err){
        if (defaults$verbose) print(err$message)
        return(NULL)})
      #V_S_inv <- tryCatch(solve(V_S), error = function(err){
      #  if (defaults$verbose) print(err$message)
      #  return(NULL)})
    } else {
      V_S_inv <- chol2inv(L)
    }
    if (is.null(V_S_inv)){
      return(NULL)
    }
    K <- V %*% V_S_inv
    if (defaults$verbose){
      eig_K <- eigen(K, symmetric = TRUE)$values
      # Small numerical error is OK
      if (max(eig_K) > 1 + 1e-6) {
        if (defaults$verbose) cat("ERROR: Genuine eigenvalue > 1\n")
      } else if (max(eig_K) > 1) {
        if (defaults$verbose) cat("WARNING: Eigenvalue slightly > 1 (numerical error)\n")
      }
    }
    
    K <- (K + t(K)) / 2 # symmetric to avoid numerical errors. # REMOVING ONLY FOR PSEUDO INVERSE # ADDING AGAIN 
    #eig_K <- eigen(K, symmetric = TRUE)
    #d_clamped <- pmax(pmin(eig_K$values, 0.95), 0.05)  # Keep in [0.1, 0.9]
    K_final <- K #eig_K$vectors %*% diag(d_clamped) %*% t(eig_K$vectors)
    #if (defaults$verbose) cat('Checking eigenvalues: ', max(eig_K$values), "\n")
#    varZ <- V - V %*% V_S_inv %*% V
    varZ <- V - K_final %*% V

    if (defaults$add.root == TRUE & defaults$model == 'BM1'){
        theta <- c(paras$theta, paras$theta)
    } else {
        theta <- paras$theta
    }

    if (defaults$weight_type == 'alt'){
        W_ <- compute_W_matrix(defaults$parsed_alt_tree, paras$alpha, add.root=defaults$add.root )
    } else {
      W_ <- weight.mat(phy, edges, paras$alpha, defaults$root.state, assume.station=defaults$assume.station)
    }
    if (dim(W_)[2] != length(theta)){
        print(dim(W_))
        print(length(theta))
        print(head(W_))
        print(theta)
        stop('Warning: Dimensions of weights matrix does not match inputted theta.\n')
    }
    E_Z <- W_ %*% theta
    E_Z <- E_Z + K_final %*%(X - E_Z)
    row.names(E_Z) <- names(X)

    if (diagnose){
      check_kalman_gain(K)
    }
    return(list(Var_Z_given_X = varZ, E_Z = E_Z))
}

loglik_OU_threepoint <- function(X, tp_paras){
  transformed.tree <- tp_paras[[1]]
  expected.vals <- tp_paras[[2]]
  N <- tp_paras[[3]]
  #print(transformed.tree$diag)
  comp <- NA
  #comp <- phylolm::three.point.compute(transformed.tree$tree, X, expected.vals, transformed.tree$diag)
  try(comp <- phylolm::three.point.compute(transformed.tree$tree, X, expected.vals, transformed.tree$diagWeight), silent=TRUE)
  if(is.na(comp[1])){
    return(-10000000) # making it negative so it switches to positive when negative returned in obejctive.
  }else{
    logl <- -as.numeric(N * log(2 * pi) + comp$logd + (comp$PP - 2 * comp$QP + comp$QQ))/2
  }
  return(logl) # returns log-likelihood NOT negative log likelihood 
}

loglik_OU_inverse <- function(X, W, V){
  N <- length(X)
  theta <- Inf
  #try(theta <- pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)%*%t(W)%*%pseudoinverse(V)%*%X, silent=TRUE)
  theta <- pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)%*%t(W)%*%pseudoinverse(V)%*%X

  if(any(theta==Inf)){
    return(-10000000)
  }
  DET <- sum(log(abs(Re(diag(qr(V)$qr)))))
  if(!is.finite(DET)){
    DET <- determinant(V, logarithm=TRUE)
    logl <- -.5*(t(X-W%*%theta)%*%pseudoinverse(V)%*%(X-W%*%theta))-.5*as.numeric(DET$modulus)-.5*(N*log(2*pi))
  }else{
    logl <- -.5*(t(X-W%*%theta)%*%pseudoinverse(V)%*%(X-W%*%theta))-.5*as.numeric(DET)-.5*(N*log(2*pi))
  }
  return(logl)
}

#### DRAFT of Joint M-step! ####
# M-step for (phi, alpha, sigma, theta)
m_step_tipfog <- function(Z, X, params_init, phy, edges, defaults, lambda1, lambda2) {

  objective <- function(params_latent) {
      par <- exp(params_latent)    
      if (defaults$model == 'BM1'){
        alpha <- 1e-10
        sigma <- par[1]
        if (defaults$assume.station == FALSE){
            theta <- c(par[2], par[2])
        } else {
            theta <- par[2]
        }
        if (defaults$skipTau){
            tau <- 0
        } else {
            tau <- par[3]
        }
      } else {
        alpha <- par[1]
        sigma <- par[2]
        if (defaults$skipTau){
            theta <- par[3:(length(par))]
            tau <- 0
        } else {
          theta <- par[3:(length(par)-1)]
          tau <- par[length(par)]
        }
        
      }

      # OU model 
      if (defaults$weight_type == 'alt'){
          W_ <- compute_W_matrix( defaults$parsed_alt_tree, alpha, add.root=defaults$add.root )
      } else {
          W_ <- weight.mat(phy, edges, alpha, defaults$root.state, assume.station=defaults$assume.station)
      }

      N <- length(phy$tip.label)
      if (defaults$algorithm == 'three.point'){
        if (defaults$transf == 'ha2014'){
            transformed.tree <- HoAne_transform_tree(phy, defaults$model, defaults$add.root, alpha, sigma)
        } else {
            transform <- defaults[['tp_const']]
            transformed.tree <- transformPhy.new(phy, transform$map, alpha, sigma, transform$tippath)
        }

        if (!defaults$skipTau){
          TIPS <- transformed.tree$tree$edge[,2] <= length(transformed.tree$tree$tip.label)
          transformed.tree$tree$edge.length[TIPS] <- transformed.tree$tree$edge.length[TIPS] + (tau^2/transformed.tree$diag/transformed.tree$diag)
        }
        

        expected.vals <- W_ %*% theta
        row.names(expected.vals) <- phy$tip.label

        tp_paras <- list(transformed.tree, expected.vals, N)
        
        ll_ou <- loglik_OU_threepoint(Z, tp_paras)
      } else {
        if (defaults$weight_type == 'alt'){
          V_ <- compute_VCV(defaults$parsed_alt_tree, alpha, sigma, add.root=defaults$add.root) # stop gap will fix. 
        } else {
          V_ <- varcov.ou(phy, edges, alpha, sigma, defaults$root.state, defaults$scaleHeight, defaults$assume.station)
        }
        
        if (!defaults$skipTau) diag(V_) <- diag(V_) + diag(tau^2, nrow(V_)) # adding tip fog in here. 
        ll_ou <- loglik_OU_inverse(Z, W_, V_) 
        
      }
    

      # Makes it comparable to likelihoods
      if (defaults$model == 'BM1'){
          lambda2_weight <- 0 
      } else {
         lambda2_weight <- N*lambda2
      }

      lambda1_weight <- N*lambda1 

      alpha_reg <- lambda2_weight * dgamma(alpha, shape = defaults$alpha_shape , rate = defaults$alpha_rate, log = TRUE)
      
      ll_prior_tau <- dnorm(tau, mean = sd(X)*defaults$tau_prior_mean, sd = defaults$tau_prior_sd, log = TRUE) # small prior on tau. -- prior parameters are hard coded here. 
      tau_reg <- lambda1_weight * ll_prior_tau

      ll_combined <- ll_ou + alpha_reg + tau_reg 

      return(-ll_combined) #minimize NLL 
    }

    p0 <- unlist(params_init)
    np <- length(p0)
    
    theta_min <- 1e-10   # theta bounds
    theta_max <- 1000
    tau_min <- 0.05
    tau_max <- 2
    if (defaults$model == 'BM1') {
        alpha_min <- 1e-10    # Only for printing info below. 
        alpha_max <- 1e-10
        sigma_min <- 0.01    # Minimum signal variance
        sigma_max <- 100.0   # Maximum signal variance
        
        if (defaults$skipTau){
            lb <- c(sigma_min, theta_min)
            ub <- c(sigma_max, theta_max)
        } else {
            lb <- c(sigma_min, theta_min, tau_min)
            ub <- c(sigma_max, theta_max, tau_max)
        }
    } else {
        alpha_min <- 0.001   # Minimum OU rate (fast selection)
        alpha_max <- 100.0   # Maximum OU rate (slow selection)
        sigma_min <- 0.01    # Minimum signal variance 
        sigma_max <- 100.0   # Maximum signal variance
        if (defaults$skipTau){
            lb <- c(alpha_min, sigma_min, rep(theta_min, np - 2)) # minus 3 for alpha, sigma, tau 
            ub <- c(alpha_max, sigma_max, rep(theta_max, np - 2))
        } else {
            lb <- c(alpha_min, sigma_min, rep(theta_min, np - 3), tau_min) # minus 3 for alpha, sigma, tau 
            ub <- c(alpha_max, sigma_max, rep(theta_max, np - 3), tau_max)          
        }
    }

    log_p0 <- check_bounds(log(p0), log(lb), log(ub), verbose=defaults$verbose) # this should return params on the log scale. 
      
    if (defaults$verbose){
        cat(sprintf("M-step 1: Optimizing with bounds...\n"))
        cat(sprintf("  Sigma bounds: [%.6f, %.6f]\n", sigma_min, sigma_max))
        cat(sprintf("  Alpha bounds: [%.6f, %.6f]\n", alpha_min, alpha_max))
        cat('Parameters:', p0, '\n')
    }
    fit <- tryCatch({
        nloptr(x0 = log_p0, # log(p0)
             eval_f = objective,
             opts = list('algorithm' = "NLOPT_LN_SBPLX", 
                         'maxeval' = '200',      # Increase from 100--brought back down to 100--back up to 200. 
                         'ftol_rel' = 1e-8,      # Tighter tolerance
                         'ftol_abs' = 1e-10),
             lb = log(lb),
             ub = log(ub))
      }, error = function(err){
        print(paste("Optimization error:", err$message))
        return(NULL)
      })
    
    if (is.null(fit)) {
        if (defaults$verbose) cat("Optimization failed, returning previous parameters\n")
        return(list(alpha = params_init$alpha, 
                  sigma = params_init$sigma, 
                  theta = params_init$theta, 
                  loglik = 10000000))
    } else {
        if (defaults$verbose) cat('Optimzation completed.\n')
    }
    
    solution <- exp(fit$solution)
    loglik <- -fit$objective
    
    # Check if hitting bounds
    if (defaults$skipTau){
        if (defaults$model == 'BM1') {
            sigma_new <- solution[1]
            theta_new <- solution[2]
            if (sigma_new <= sigma_min * 1.01) {
                if (defaults$verbose) cat("WARNING: Sigma hit lower bound\n")
            }
            if (sigma_new >= sigma_max * 0.99) {
                if (defaults$verbose) cat("WARNING: Sigma hit upper bound\n")
            }
              alpha_new <- 1e-10
          } else {
              alpha_new <- solution[1]
              sigma_new <- solution[2]
              theta_new <- solution[3:np]
              # Report what happened    
              if (defaults$verbose & sigma_new <= sigma_min * 1.01) {
                  cat("WARNING: Sigma hit lower bound\n")
              }
          }
          tau_new <- 0
      } else {
        if (defaults$model == 'BM1') {
            sigma_new <- solution[1]
            theta_new <- solution[2]
            tau_new <- solution[3]
            if (sigma_new <= sigma_min * 1.01) {
                if (defaults$verbose) cat("WARNING: Sigma hit lower bound\n")
            }
            if (sigma_new >= sigma_max * 0.99) {
                if (defaults$verbose) cat("WARNING: Sigma hit upper bound\n")
            }
            alpha_new <- 1e-10
      } else {
            alpha_new <- solution[1]
            sigma_new <- solution[2]
            theta_new <- solution[3:(np-1)]
            tau_new <- solution[np]
            # Report what happened    
            if (defaults$verbose & sigma_new <= sigma_min * 1.01) {
              cat("WARNING: Sigma hit lower bound\n")
            }
        }
    }
  

    if (defaults$verbose) cat(sprintf("M-step 1 result: alpha=%.6f, sigma=%.6f, loglik=%.2f\n",alpha_new, sigma_new, loglik))

    if (defaults$weight_type == 'alt'){
        W_ <- compute_W_matrix( defaults$parsed_alt_tree, alpha_new, add.root=defaults$add.root )
    } else {
        W_ <- weight.mat(phy, edges, alpha_new, defaults$root.state, assume.station=defaults$assume.station)
    }
    if (defaults$add.root == TRUE & defaults$model == 'BM1'){
        theta.tmp <- c(theta_new, theta_new)
    } else {
        theta.tmp <-  theta_new
    }

    N <- length(phy$tip.label)
    if (defaults$algorithm == 'three.point'){
        if (defaults$transf == 'ha2014'){
            transformed.tree <- HoAne_transform_tree(phy, defaults$model, defaults$add.root, alpha_new, sigma_new)
        } else {
            transform <- defaults[['tp_const']]
            transformed.tree <- transformPhy.new(phy, transform$map, alpha_new, sigma_new, transform$tippath)
        }

        TIPS <- transformed.tree$tree$edge[,2] <= length(transformed.tree$tree$tip.label)
        transformed.tree$tree$edge.length[TIPS] <- transformed.tree$tree$edge.length[TIPS] + (tau_new^2/transformed.tree$diag/transformed.tree$diag)

        expected.vals <- W_ %*% theta.tmp
        row.names(expected.vals) <- phy$tip.label
        
        tp_paras_final <- list(transformed.tree, expected.vals, N)
        
        ll_ou_final <- loglik_OU_threepoint(Z, tp_paras_final)
    } else {
        if (defaults$weight_type == 'alt'){
          V_ <- compute_VCV(defaults$parsed_alt_tree, alpha_new, sigma_new, add.root=defaults$add.root) # stop gap will fix. 
        } else {
          V_ <- varcov.ou(phy, edges, alpha_new, sigma_new, defaults$root.state, defaults$scaleHeight, defaults$assume.station)
        }
       # V_ <- compute_VCV(defaults$parsed_alt_tree, alpha_new, sigma_new, add.root=defaults$add.root) # stop gap will fix. 
        diag(V_) <- diag(V_) + diag(tau_new^2, nrow(V_)) # adding tip fog in here. 
        ll_ou_final <- loglik_OU_inverse(Z, W_, V_) 
    }

    return(list(alpha = alpha_new, 
      sigma = sigma_new, 
      theta = theta_new, 
      tau=tau_new, 
      ll_ou = ll_ou_final, 
      ll_total = ll_ou_final, 
      objective = loglik))
}

# ================================================================
# FULL EM
# ================================================================
run_em <- function(phy, X, params_init, edges, defaults,  
  max_iter = 150, 
  tol_abs = 0.1, 
  tol_rel = 1e-4, 
  lambda1=0.2, # weight on log_prior for noise likelihood in optimization  
  lambda2=0.1, 
  tau_prior_mean = 0.2, 
  tau_prior_sd = 0.1, 
  alpha_shape = 2, 
  alpha_rate = 1, 
  diagnose=TRUE, 
  save=TRUE, 
  save_pre=NULL,
  skipE = FALSE, 
  save_history = FALSE, 
  verbose=TRUE) {

# Defaults should contain: root state, assume station, scale height 
  params <- params_init

  defaults$verbose <- verbose
  defaults$tau_prior_sd <- tau_prior_sd
  defaults$tau_prior_mean <- tau_prior_mean
  defaults$inits <- params_init
  defaults$alpha_shape <- alpha_shape 
  defaults$alpha_rate <- alpha_rate 

  if (verbose) cat('Starting parameters:\n', paste(params, collapse = " "))
  if (any(names(X) != phy$tip.label)){
    stop('WARNING: tip labels dont match data labels. Redo.')
  }

  state_names <- names(params_init$theta)
  history <- data.frame(
    iter = integer(),
    ll_total = numeric(),
    ll_ou = numeric(),
    ll_noise = numeric(),
    sigma = numeric(),
    alpha = numeric(),
    tau = numeric()
  )
  monitor_history <- NULL

  if (defaults$model == 'BM1'){
      params$alpha <- 1e-10
  }

  if (skipE){
    max_iter <- 1 
  }

  for (iter in seq_len(max_iter)) {
    #cat('iter:', iter, '\n')
    tmp <- params

    # E-step
    if (skipE){
      e_result <- list(E_Z = X)
    } else {
      if (verbose) cat('Starting E-step...\n')
      e_result <- e_step(X, phy, edges, params, defaults, diagnose)
      if (is.null(e_result)){
        if (verbose) cat('E_step failed for iter', iter, '. ')
        next
      }
    }
          
    if (defaults$model == 'BM1'){
        pset <- c("sigma", "theta", "tau")
    } else {
        pset <- c("alpha", "sigma", "theta", "tau")
    }

    if (defaults$skipTau){
      pset <- pset[1:(length(pset)-1)] # remove tau 
    }

    m_step1_result <- m_step_tipfog(e_result$E_Z, X, params[pset], phy, edges, defaults, lambda1, lambda2)
    params$alpha <- m_step1_result$alpha
    params$sigma <- m_step1_result$sigma
    params$theta <- m_step1_result$theta 
    params$tau <- m_step1_result$tau  

    ll_total <- m_step1_result$ll_total
    ll1 <- m_step1_result$ll_ou
    
    # Recover tp likelihood given params 
    if (any(is.null(params)) || any(is.na(params))){
        if (verbose){
          cat('Iteration failed to produce estimates..\n')
          print(params)          
        }
        break
    } 

    history <- rbind(history, data.frame(
      iter = iter,
      ll_total = ll_total,
      objective = m_step1_result$objective, 
      sigma = params$sigma,
      alpha = params$alpha,
      tau = params$tau #,
    ))
    
        # Print progress
    if (verbose) {
      cat(sprintf("\nIter %3d | LL_total = %10.2f | sigma = %.6f | alpha = %.6f || tau = %.6f\n", iter, ll_total, params$sigma, params$alpha,params$tau))
      cat(paste0(state_names, ': ', round(params$theta, 3) ), '\n')
      flush.console()
    }
      
    if (diagnose){
        diag_result <- diagnose_residuals(X, e_result)
    }
    # Convergence: relative change must be small AND absolute change must be small
    # Let run for at least three iterations. 
    convergence <- NA
    if (iter > 3) {
      abs_change <- ll_total - ll_prev
      rel_change <- abs(abs_change) / (abs(ll_prev) + 1)
      if (verbose) cat('\nIter', iter, '|| Total Change:', abs_change,'| Relative Change:', rel_change,'\n')
      if (abs_change < tol_abs && rel_change < tol_rel) {
        if (verbose) cat(sprintf("Converged at iteration %d\nFinal sigma: %.8f\n", iter, params$sigma))
        convergence <- 'converged'
        break
      }

      if (iter-10 > 0){
        monitor_oscillating <- is_oscillating(history[(iter-10):iter, ]$ll_total, min_cycles=5) # generous amount of cycles. 
        if (monitor_oscillating) {
            if (verbose) cat(sprintf("Early stop due to oscillating likelihoods at iteration %d\nFinal sigma: %.8f\n", iter, params$sigma))
            convergence <- 'early_stop_oscillating'
            break
        }
      } 
    } 

    ll_prev <- ll_total
  }

  defaults$ntips <- length(phy$tip.label)
  history$converge <- convergence
  if (save){
      if (is.null(save_pre)){
        save_pre <- paste(sample(c(0:9, letters, LETTERS), 8, replace = TRUE), collapse = "")
      } 
      result <- list(paras = params, e_result = e_result, history = history,  defaults = defaults)
      filename <- paste0(save_pre, '.rds')
      saveRDS(result, filename)
      if (verbose) cat('Results saved to RDS:', filename, '\n')
  }

  if (save_history){
      return(list(paras = params, final = history[iter, ], e_Z = e_result$E_Z, history = history, defaults = defaults))
  } else {
    return(list(paras = params, final = history[iter, ], e_Z = e_result$E_Z, defaults = defaults))
  }
}
