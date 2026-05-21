
####### EM #######
# Proceeding this: 
# Need: estimates for alpha, sigma, theta , tau 
# Constants/defaults: edges, phy, root.state, assume.station 
e_step <- function(X, phy, edges, paras, defaults, diagnose=TRUE){
    #V <- quickVCV(phy, paras$alpha, paras$sigma, defaults$scaleHeight)
    V <- varcov.ou(phy, edges, paras$alpha, paras$sigma, defaults$root.state, defaults$scaleHeight, defaults$assume.station)
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
    
    #K <- (K + t(K)) / 2 # symmetric to avoid numerical errors. # REMOVING ONLY FOR PSEUDO INVERSE 
    varZ <- V - V %*% V_S_inv %*% V

    if (defaults$assume.station == FALSE & defaults$model == 'BM1'){
        theta <- c(paras$theta, paras$theta)
    } else {
        theta <- paras$theta
    }

    if (defaults$weight_type == 'alt'){
        W_ <- compute_W_matrix( defaults$parsed_alt_tree, paras$alpha)
    } else {
      W_ <- weight.mat(phy, edges, paras$alpha, defaults$root.state, assume.station=defaults$assume.station)
    }
    if (dim(W_)[2] != length(theta)){
        print(dim(W_))
        print(length(theta))
        stop('Warning: Dimensions of weights matrix does not match inputted theta.\n')
    }
    E_Z <- W_ %*% theta
    E_Z <- E_Z + K%*%(X - E_Z)
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
  try(comp <- phylolm::three.point.compute(transformed.tree$tree, X, expected.vals, transformed.tree$diag), silent=TRUE)
  if(is.na(comp[1])){
    return(-10000000) # making it negative so it switches to positive when negative returned in obejctive.
  }else{
    logl <- -as.numeric(N * log(2 * pi) + comp$logd + (comp$PP - 2 * comp$QP + comp$QQ))/2
  }
  return(logl) # returns log-likelihood NOT negative log likelihood 
}

loglik_Gaus_noise <- function(X, Z, tau) {
  # X: n x p observed data
  # Z: n x p latent true values
  # tau: noise std deviation (scalar or p-vector)
  residuals <- X - Z
  ll <- sum(dnorm(residuals, mean = 0, sd = tau, log = TRUE))
  return(ll)
}

loglik_NegBin_noise <- function(X, Z, r) {
  # X: n x p observed data
  # Z: n x p latent true values
  # r: noise dispersion parameter for NB. 
  ll <- lgamma(X + r) - lgamma(r) - lgamma(X + 1) + X * log(Z / (r + Z)) + r * log(r / (r + Z))
  return(ll)
}


# M-step for (phi, alpha, sigma, theta)
m_step_latent <- function(E_Z, params_init, phy, edges, defaults) {
  
  objective <- function(params_latent, Z) {
    par <- exp(params_latent)    
    if (defaults$model == 'BM1'){
      alpha <- 1e-10
      sigma <- par[1]
      if (defaults$assume.station == FALSE){
          theta <- c(par[2], par[2])
      } else {
          theta <- par[2]
      }
    } else {
      alpha <- par[1]
      sigma <- par[2]
      theta <- par[3:length(par)]
    }

    if (defaults$weight_type == 'alt'){
        W_ <- compute_W_matrix( defaults$parsed_alt_tree, alpha)
    } else {
      W_ <- weight.mat(phy, edges, 
        alpha, defaults$root.state, assume.station=defaults$assume.station)
    }
    transform <- defaults[['tp_const']]
    transformed.tree <- transformPhy.new(phy, transform$map, alpha, sigma, transform$tippath)
    #expected.vals <- W_ %*% theta
    #row.names(expected.vals) <- phy$tip.label
    expected.vals <- colSums(t(W_) * theta)
    names(expected.vals) <- phy$tip.label 
    N <- length(phy$tip.label)
    tp_paras <- list(transformed.tree, expected.vals, N)
    
    ll <- loglik_OU_threepoint(Z, tp_paras)
    #cat('Monitoring - sigma:', sigma, 'theta:', theta, '|| likelihood:', ll, '\n')
    return(-ll)
  }

    p0 <- unlist(params_init)
    np <- length(p0)
  
    # FIXED: More reasonable bounds
    # These should match biological understanding
    theta_min <- 1e-10   # theta bounds
    theta_max <- 100
    if (defaults$model == 'BM1') {
      alpha_min <- 1e-10    # Only for printing info below. 
      alpha_max <- 1e-10
      sigma_min <- 0.01    # Minimum signal variance
      sigma_max <- 100.0   # Maximum signal variance
      lb <- c(sigma_min, theta_min)
      ub <- c(sigma_max, theta_max)
    } else {
      alpha_min <- 0.001   # Minimum OU rate (fast selection)
      alpha_max <- 100.0   # Maximum OU rate (slow selection)
      sigma_min <- 0.01    # Minimum signal variance 
      sigma_max <- 100.0   # Maximum signal variance
      lb <- c(alpha_min, sigma_min, rep(theta_min, np - 2))
      ub <- c(alpha_max, sigma_max, rep(theta_max, np - 2))
    }
    
    if (defaults$verbose){
      cat(sprintf("M-step 1: Optimizing with bounds...\n"))
      cat(sprintf("  Sigma bounds: [%.6f, %.6f]\n", sigma_min, sigma_max))
      cat(sprintf("  Alpha bounds: [%.6f, %.6f]\n", alpha_min, alpha_max))
      cat('Parameters:', p0, '\n')
    }
    fit <- tryCatch({
      nloptr(x0 = log(p0),
           eval_f = objective,
           Z = E_Z,
           opts = list('algorithm' = "NLOPT_LN_SBPLX", 
                       'maxeval' = '1000',      # Increase from 100
                       #'ftol_rel' = 1e-8,      # Tighter tolerance
                       "ftol_rel"=.Machine$double.eps^0.5, # bringing back the OUwie settings for this just to see. 
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
    }
  
    solution <- exp(fit$solution)
    loglik <- -fit$objective
  
    # Check if hitting bounds
    if (defaults$model == 'BM1') {
        sigma_new <- solution[1]
        theta_new <- solution[2]
        if (sigma_new <= sigma_min * 1.01) {
          if (defaults$verbose) cat("WARNING: Sigma hit lower bound\n")
        }
        if (sigma_new >= sigma_max * 0.99) {
          if (defaults$verbose) cat("WARNING: Sigma hit upper bound\n")
        }
        return(list(alpha = 1e-10, sigma = sigma_new, theta = theta_new, loglik = loglik))
    } else {
      alpha_new <- solution[1]
      sigma_new <- solution[2]
      theta_new <- solution[3:np]
    
      # Report what happened
      if (defaults$verbose) cat(sprintf("M-step 1 result: alpha=%.6f, sigma=%.6f, loglik=%.2f\n", alpha_new, sigma_new, loglik))
    
      if (defaults$verbose & sigma_new <= sigma_min * 1.01) {
        cat("WARNING: Sigma hit lower bound\n")
      }
      return(list(alpha = alpha_new, sigma = sigma_new, theta = theta_new, loglik = loglik))
    }
}

#### DRAFT of Joint M-step! ####
# M-step for (phi, alpha, sigma, theta)
m_step_joint <- function(Z, X, params_init, phy, edges, defaults, beta, lambda1) {

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
      tau <- par[3]
    } else {
      alpha <- par[1]
      sigma <- par[2]
      theta <- par[3:(length(par)-1)]
      tau <- par[length(par)]
    }
    # Quick check 
    if (sigma <= 0 || tau <= 0 || alpha < 0) return(1e10)

    # Noise model
    residuals <- X - Z
    ll_noise <- sum(dnorm(residuals, mean = 0, sd = tau, log = TRUE))
    
    # OU model 
    if (defaults$weight_type == 'alt'){
        W_ <- compute_W_matrix( defaults$parsed_alt_tree, alpha)
    } else {
      W_ <- weight.mat(phy, edges, alpha, defaults$root.state, assume.station=defaults$assume.station)
    }
    transform <- defaults[['tp_const']]
    transformed.tree <- transformPhy.new(phy, transform$map, alpha, sigma, transform$tippath)
    expected.vals <- W_ %*% theta
    row.names(expected.vals) <- phy$tip.label
    N <- length(phy$tip.label)
    tp_paras <- list(transformed.tree, expected.vals, N)
    
    ll_ou <- loglik_OU_threepoint(Z, tp_paras)
    #cat('Monitoring - sigma:', sigma, 'theta:', theta, '|| likelihood:', ll, '\n')

    # Combine   
    #ll_combined <- ll_noise + ll_ou
    ll_combined <- ll_ou + beta * ll_noise 
    ll_prior_tau <- dnorm(tau, mean = sd(X)*defaults$tau_prior_mean, sd = defaults$tau_prior_sd, log = TRUE) # small prior on tau. -- prior parameters are hard coded here. 

    return(-(ll_combined + lambda1 * ll_prior_tau ))
  }

    p0 <- unlist(params_init)
    np <- length(p0)
  
    # FIXED: More reasonable bounds
    # These should match biological understanding
    theta_min <- 1e-10   # theta bounds
    theta_max <- 100
    tau_min <- 0.05
    tau_max <- 2
    if (defaults$model == 'BM1') {
      alpha_min <- 1e-10    # Only for printing info below. 
      alpha_max <- 1e-10
      sigma_min <- 0.01    # Minimum signal variance
      sigma_max <- 100.0   # Maximum signal variance
      lb <- c(sigma_min, theta_min, tau_min)
      ub <- c(sigma_max, theta_max, tau_max)
    } else {
      alpha_min <- 0.001   # Minimum OU rate (fast selection)
      alpha_max <- 100.0   # Maximum OU rate (slow selection)
      sigma_min <- 0.01    # Minimum signal variance 
      sigma_max <- 100.0   # Maximum signal variance
      lb <- c(alpha_min, sigma_min, rep(theta_min, np - 3), tau_min) # minus 3 for alpha, sigma, tau 
      ub <- c(alpha_max, sigma_max, rep(theta_max, np - 3), tau_max)
    }

    # Fix bounds for tau 
    if (p0[np] < theta_min){
      p0[np] <- theta_min 
    } else if (p0[np] > tau_max){
      p0[np] <- tau_max
    }
    
    if (defaults$verbose){
      cat(sprintf("M-step 1: Optimizing with bounds...\n"))
      cat(sprintf("  Sigma bounds: [%.6f, %.6f]\n", sigma_min, sigma_max))
      cat(sprintf("  Alpha bounds: [%.6f, %.6f]\n", alpha_min, alpha_max))
      cat('Parameters:', p0, '\n')
    }
    fit <- tryCatch({
      nloptr(x0 = log(p0),
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

    if (defaults$verbose) cat(sprintf("M-step 1 result: alpha=%.6f, sigma=%.6f, loglik=%.2f\n",alpha_new, sigma_new, loglik))

    ### Final likelihood for monitoring 
    residuals <- X - Z
    ll_noise_final <- sum(dnorm(residuals, mean = 0, sd = tau_new, log = TRUE))
    
    if (defaults$weight_type == 'alt'){
        W_ <- compute_W_matrix( defaults$parsed_alt_tree, alpha_new)
    } else {
      W_ <- weight.mat(phy, edges, alpha_new, defaults$root.state, assume.station=defaults$assume.station)
    }
    if (defaults$assume.station == FALSE & defaults$model == 'BM1'){
      theta.tmp <- c(theta_new, theta_new)
    } else {
      theta.tmp <-  theta_new
    }
    expected_vals <- W_ %*% theta.tmp
    row.names(expected_vals) <- phy$tip.label
    transformed.tree <- transformPhy.new(phy, defaults$tp_const$map, alpha_new, sigma_new, defaults$tp_const$tippath)
    expected.vals <- W_ %*% theta.tmp
    row.names(expected.vals) <- phy$tip.label
    N <- length(phy$tip.label)
    tp_paras_final <- list(transformed.tree, expected.vals, N)

    ll_ou_final <- loglik_OU_threepoint(Z, tp_paras_final)


    return(list(alpha = alpha_new, 
      sigma = sigma_new, 
      theta = theta_new, 
      tau=tau_new, 
      ll_ou = ll_ou_final, 
      ll_noise = ll_noise_final, 
      ll_total = ll_ou_final + ll_noise_final))
}


#### A BUNCH OF TAU M-STEPS #### 
# M-step for tau
m_step_noise <- function(X, E_Z, Var_Z_given_X) {
  # Maximize E[log P(X | Z, tau) | X]
  # Under Gaussian model: optimal tau is residual std dev
  
  residuals <- X - E_Z
  # Add correction for uncertainty in Z
  tau_new <- sqrt(mean(residuals^2) + mean(diag(Var_Z_given_X)))
  ll <- sum(dnorm(residuals, mean = 0, sd = tau_new, log = TRUE))

  return(list('tau' = tau_new, 'loglik' = ll))
}

m_step_noise.Reg <- function(X, E_Z, Var_Z_given_X, lambda = 0.1) {
  
  # Original noise estimate
  residuals <- X - E_Z
  tau_unpenalized <- sqrt(mean(residuals^2) + mean(diag(Var_Z_given_X)))
  
  # Add L2 penalty: encourages tau to be small
  # Penalized estimate: lambda * tau^2 term shrinks tau toward 0
  tau_squared <- mean(residuals^2) + mean(diag(Var_Z_given_X))
  
  # Penalized objective: minimize ||residuals||^2 + lambda * tau^2
  # Solution: tau_new = sqrt(tau_squared / (1 + lambda))
  tau_new <- sqrt(tau_squared / (1 + lambda))
  ll <- sum(dnorm(residuals, mean = 0, sd = tau_new, log = TRUE))

  return(list('tau' = tau_new, 'loglik'=ll))
}

m_step_noise_with_prior <- function(X, E_Z, Var_Z_given_X,
                                     tau_prior_mean = 0.2,
                                     tau_prior_sd = 0.05) {
  
  # Data-driven estimate
  residuals <- X - E_Z
  tau_squared <- mean(residuals^2) + mean(diag(Var_Z_given_X))
  tau_data <- sqrt(tau_squared)
  
  # Combine data likelihood + prior belief
  data_weight <- 1 / (tau_prior_sd^2)
  prior_weight <- 1 / (tau_prior_sd^2)
  
  tau_new <- (data_weight * tau_data + prior_weight * tau_prior_mean) / 
             (data_weight + prior_weight)
  
  tau_new <- pmax(tau_new, 0.01)
 # print(residuals)
  
  ll <- sum(dnorm(residuals, mean = 0, sd = tau_new, log = TRUE))
  
  return(list('tau' = tau_new, 'loglik' = ll))
}

m_step_noise_informative <- function(X, E_Z, Var_Z_given_X,
                                      tau_prior_mean,
                                      tau_prior_strength = 100) {
  
  residuals <- X - E_Z
  tau_squared_data <- mean(residuals^2) + mean(diag(Var_Z_given_X))
  
  # Weighted average: strong prior prevents collapse
  tau_squared_new <- (tau_squared_data + tau_prior_strength * tau_prior_mean^2) /
                     (1 + tau_prior_strength)
  
  tau_new <- sqrt(tau_squared_new)
  
  # Enforce floor
  tau_new <- pmax(tau_new, tau_prior_mean * 0.5)
  ll <- sum(dnorm(residuals, mean = 0, sd = tau_new, log = TRUE))

  return(list('tau' = tau_new, 'loglik'=ll))
}


#### DRAFT of optimized noise model ####
# One idea is to keep this design so that we can sub in an NB based M step. 
m_step_noise_optimized <- function(X, E_Z, Var_Z_given_X,
                                    tau_prior_mean,
                                    tau_prior_strength = 100) {
  
  # Objective: maximize LL_noise + log-prior
  objective <- function(tau) {
    
    if (tau <= 0) return(1e10)  # Invalid
    
    residuals <- X - E_Z
    
    # Data likelihood: how well tau explains residuals
    ll_data <- sum(dnorm(residuals, mean = 0, sd = tau, log = TRUE))
    
    # Prior likelihood: how close tau is to prior_mean
    # Using a Gaussian prior centered on tau_prior_mean
    ll_prior <- dnorm(tau, 
                     mean = tau_prior_mean,
                     sd = tau_prior_mean / sqrt(tau_prior_strength),
                     log = TRUE)
    
    # Combined objective
    return(-(ll_data + ll_prior))  # Negative for minimization
  }
  
  # Optimize tau
  fit <- optimize(objective, 
                 interval = c(0.01, 2.0),
                 maximum = FALSE)  # minimize negative LL = maximize LL
  
  tau_new <- fit$minimum
  
  # Compute final likelihood
  residuals <- X - E_Z
  ll_data <- sum(dnorm(residuals, mean = 0, sd = tau_new, log = TRUE))
  ll_prior <- dnorm(tau_new, 
                   mean = tau_prior_mean,
                   sd = tau_prior_mean / sqrt(tau_prior_strength),
                   log = TRUE)
  ll_total <- ll_data + ll_prior
  
  return(list('tau' = tau_new, 'loglik' = ll_data))  # Return data LL
}


# ================================================================
# FULL EM
# ================================================================
log_progress <- function(iter, ll_hist, ll1, ll2) {
  cat('Current likelihood for iter', iter, '=',ll_hist[iter], '|| Noise LL:', ll2 , 'OU LL:', ll1, '\n') 
  flush.console()
}

log_params<- function(iter, params) {
  cat('Current estimates for iter', iter, '|| Alpha:', params$alpha , 'Sigma:', params$sigma, 'Tau:', params$tau, '\n') 
  flush.console()
}

run_em <- function(phy, X, params_init, edges, defaults,  
  max_iter = 100, 
  mmode = 'two_step', 
  tau_prior_strength= 100, # for two step 
  tol_abs = 0.1, 
  tol_rel = 1e-4, 
  burn_in=10, # for joint: how quickly to bring ll_noise into the optimization 
  lambda1=0.2, # for joint: weight on log_prior for noise likelihood in optimization  
  diagnose=TRUE, 
  save=TRUE, 
  save_pre=NULL,
  verbose=TRUE) {

# Defaults should contain: root state, assume station, scale height 
  params <- params_init

  # Adding to defaults 
  defaults$verbose <- verbose
  if (is.null(defaults$tau_prior_mean)) defaults$tau_prior_mean <- 0.2
  if (is.null(defaults$tau_prior_sd)) defaults$tau_prior_mean <- 0.1

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
    tau = numeric(),
    beta = numeric(),
    mse = numeric()
  )
  if (defaults$model == 'BM1'){
      params$alpha <- 1e-10
  }
  for (iter in seq_len(max_iter)) {
    #cat('iter:', iter, '\n')
    tmp <- params

    # E-step
    if (verbose) cat('Starting E-step...\n')
    e_result <- e_step(X, phy, edges, params, defaults, diagnose)
    if (is.null(e_result)){
      if (verbose) cat('E_step failed for iter', iter, '. ')
      next
    }
    #print(head(e_result$E_Z))
    # M-step (A): latent parameters
    #cat('Starting M-step 1...\n')
    if (mmode == 'two_step'){
      if (defaults$model == 'BM1'){
        pset <- c("sigma", "theta")
      } else {
        pset <- c("alpha", "sigma", "theta")
      }
      #print(params)
      m_step1_result <- m_step_latent(e_result$E_Z, params[pset], phy, edges, defaults)
      params[pset] <- m_step1_result[pset]
      params$alpha <- m_step1_result$alpha
      params$sigma <- m_step1_result$sigma
      params$theta <- m_step1_result$theta 

      # M-step (B): noise parameter
      #cat('Starting M-step 2...\n')
      #m_step2_result <- m_step_noise.Reg(X, e_result$E_Z, e_result$Var_Z_given_X)#, lambda = lambda )
      #m_step2_result <- m_step_noise.Reg(X, e_result$E_Z, e_result$Var_Z_given_X, lambda) # Rolling with defaults for tau prior for now. 
      #m_step2_result <- m_step_noise_with_prior(X, e_result$E_Z, e_result$Var_Z_given_X,
      #                                  tau_prior_mean = sd(X) * 0.3, ## Trying with smaller priors for BM1. 0.3
      #                                  tau_prior_sd = sd(X) * 0.1) # 0.1 

      if (verbose) cat('Tau prior mean:', sd(X) *  defaults$tau_prior_mean, '\n')
      m_step2_result <-  m_step_noise_informative(X,  
        e_result$E_Z, 
        e_result$Var_Z_given_X, 
        tau_prior_mean = sd(X) * defaults$tau_prior_mean, 
        tau_prior_strength = tau_prior_strength) # I wonder if this is going to be too strong? 
      params$tau <- m_step2_result$tau

      ll2 <- m_step2_result$loglik; ll1 <- m_step1_result$loglik
      ll_total <- ll2 + ll1
      beta<-1 # just to satisfy stop condition later on. 

    } else if (mmode == 'joint'){
      if (defaults$model == 'BM1'){
        pset <- c("sigma", "theta", "tau")
      } else {
        pset <- c("alpha", "sigma", "theta", "tau")
      }

      #progress <- (iter - max_iter/3) / (max_iter/3)
      #weight_ou <- 0.5 * exp(-progress)
      #weight_noise <- 1 - weight_ou
      # Annealing schedule: gradually shift from equal to noise-favoring
      # annealing schedule
      if (is.null(burn_in)) {
        beta <- 1 # if no burn in then beta stays 1. 
      } else {
        #midpoint <- max_iter / 4
        #beta <- 1 / (1 + exp(-(iter - midpoint) / burn_in)) # sigmoidal burning 
        #if (beta >= 0.999) beta <- 1
        beta <- min(1.0, iter / burn_in) # when iterations get above burn in then weighted to 1. 

      }

      if (verbose) cat(sprintf("Annealing: beta=%.2f\n", beta))


      m_step1_result <- m_step_joint(e_result$E_Z, X, params[pset], phy, edges, defaults, beta, lambda1)
      params$alpha <- m_step1_result$alpha
      params$sigma <- m_step1_result$sigma
      params$theta <- m_step1_result$theta 
      params$tau <- m_step1_result$tau  

      ll_total <- m_step1_result$ll_total
      ll1 <- m_step1_result$ll_ou
      ll2 <- m_step1_result$ll_noise 

    }
    
    # Recover tp likelihood given params 
    if (any(is.null(params)) || any(is.na(params))){
        if (verbose){
          cat('Iteration failed to produce estimates..\n')
          print(params)          
        }
        break
    } 

    MSE <- mean((X - e_result$E_Z)^2)
    history <- rbind(history, data.frame(
      iter = iter,
      ll_total = ll_total,
      ll_ou = ll1,
      sigma = params$sigma,
      alpha = params$alpha,
      ll_noise = ll2, 
      tau = params$tau ,
      beta = ifelse(is.null(beta), NA, beta),
      mse = MSE
    ))
    
        # Print progress
    if (verbose) {
      cat(sprintf("\nIter %3d | LL_total = %10.2f | LL_OU = %10.2f | sigma = %.6f | alpha = %.6f || LL_NOISE = %10.2f | tau = %.6f | MSE = %.6f\n", iter, ll_total, ll1, params$sigma, params$alpha, ll2, params$tau, MSE))
      cat(paste0(state_names, ': ', round(params$theta, 3) ), '\n')
      flush.console()
    }
      
    if (diagnose){
        diag_result <- diagnose_residuals(X, e_result)
    }
    # Convergence: relative change must be small AND absolute change must be small
    if (iter > 1 & beta >= 1) {
    #if (beta > 1){
      abs_change <- ll_total - ll_prev
      rel_change <- abs(abs_change) / (abs(ll_prev) + 1)
      if (verbose) cat('\nIter', iter, '|| Total Change:', abs_change,'| Relative Change:', rel_change,'\n')
      if (abs_change < tol_abs && rel_change < tol_rel) {
        cat(sprintf("Converged at iteration %d\n", iter))
        cat(sprintf("Final sigma: %.8f\n", params$sigma))  
        break
      }
    }
    ll_prev <- ll_total
  }

  defaults$ntips <- length(phy$tip.label)

  if (save){
      if (is.null(save_pre)){
        save_pre <- paste(sample(c(0:9, letters, LETTERS), 8, replace = TRUE), collapse = "")
      } 
      result <- list(paras = params, e_result = e_result, history = history,  defaults = defaults)
      filename <- paste0(save_pre, '.rds')
      saveRDS(result, filename)
      if (verbose) cat('Results saved to RDS:', filename, '\n')
  }

  return(list(paras = params, e_result = e_result, history = history, defaults = defaults))
}
