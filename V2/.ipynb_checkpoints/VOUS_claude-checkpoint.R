###############################################################################
# VOUS: Variational Ornstein-Uhlenbeck Stochastics
# R Implementation
#
# Models stochastic single-cell gene expression over inferred cell lineage
# trees using an OU process, softplus transformation, and Negative Binomial
# observation model with mean-field variational inference.
###############################################################################

library(ape)       # For phylogenetic tree manipulation
library(Matrix)    # For sparse matrices
library(stats)     # For optimization

###############################################################################
# 1. TREE UTILITIES
###############################################################################

#' Parse a lineage tree and extract structural information
#' 
#' @param tree An object of class "phylo" (from ape package)
#' @param regime_labels Named vector: leaf_name -> regime_label
#' @return List with tree structure information
parse_tree <- function(tree, regime_labels) {
  n_leaves <- length(tree$tip.label)
  n_nodes <- tree$Nnode + n_leaves
  
  # Get distances from root to each node
  root_dists <- node.depth.edgelength(tree)
  
  # Get leaf indices
  leaf_indices <- 1:n_leaves
  
  # Assign regimes to internal nodes using Fitch-Hartigan (parsimony)
  # For simplicity, we use a basic ancestral reconstruction
  node_regimes <- assign_regimes_fitch(tree, regime_labels)
  
  # Compute shared path lengths (root to MRCA) for each pair of leaves
  # This is needed for the OU covariance matrix
  mrca_matrix <- mrca(tree)  # matrix of MRCA node indices
  shared_lengths <- matrix(0, n_leaves, n_leaves)
  for (i in 1:n_leaves) {
    for (j in i:n_leaves) {
      mrca_node <- mrca_matrix[i, j]
      shared_lengths[i, j] <- root_dists[mrca_node]
      shared_lengths[j, i] <- shared_lengths[i, j]
    }
  }
  
  # Root-to-leaf distances
  leaf_dists <- root_dists[leaf_indices]
  
  # Build regime path information for weight matrix W
  regime_paths <- build_regime_paths(tree, node_regimes, root_dists)
  
  # Unique regime names
  unique_regimes <- unique(regime_labels)
  
  list(
    tree = tree,
    n_leaves = n_leaves,
    n_nodes = n_nodes,
    leaf_indices = leaf_indices,
    leaf_dists = leaf_dists,
    shared_lengths = shared_lengths,
    node_regimes = node_regimes,
    regime_paths = regime_paths,
    unique_regimes = unique_regimes,
    regime_labels = regime_labels,
    root_dists = root_dists
  )
}

#' Simple Fitch-Hartigan parsimony assignment of regimes to internal nodes
#' 
#' @param tree phylo object
#' @param tip_regimes Named vector of regime labels for tips
#' @return Named vector of regime labels for all nodes
assign_regimes_fitch <- function(tree, tip_regimes) {
  n_tips <- length(tree$tip.label)
  n_nodes <- tree$Nnode + n_tips
  
  # Initialize
  node_states <- vector("list", n_nodes)
  node_regime <- character(n_nodes)
  
  # Assign tip states
  for (i in 1:n_tips) {
    tip_name <- tree$tip.label[i]
    node_states[[i]] <- tip_regimes[tip_name]
    node_regime[i] <- tip_regimes[tip_name]
  }
  
  # Bottom-up pass (Fitch)
  # Process internal nodes in post-order
  post_order <- postorder(tree)
  
  for (idx in seq_along(post_order)) {
    edge_idx <- post_order[idx]
    parent <- tree$edge[edge_idx, 1]
    child <- tree$edge[edge_idx, 2]
    
    if (is.null(node_states[[parent]])) {
      node_states[[parent]] <- node_states[[child]]
    } else {
      intersection <- intersect(node_states[[parent]], node_states[[child]])
      if (length(intersection) > 0) {
        node_states[[parent]] <- intersection
      } else {
        node_states[[parent]] <- union(node_states[[parent]], node_states[[child]])
      }
    }
  }
  
  # Top-down pass (Hartigan)
  root <- n_tips + 1
  node_regime[root] <- node_states[[root]][1]
  
  pre_order_edges <- rev(postorder(tree))
  for (idx in seq_along(pre_order_edges)) {
    edge_idx <- pre_order_edges[idx]
    parent <- tree$edge[edge_idx, 1]
    child <- tree$edge[edge_idx, 2]
    
    if (node_regime[parent] %in% node_states[[child]]) {
      node_regime[child] <- node_regime[parent]
    } else {
      node_regime[child] <- node_states[[child]][1]
    }
  }
  
  names(node_regime) <- 1:n_nodes
  node_regime
}

#' Build regime path information for the weight matrix W
#' 
#' @param tree phylo object
#' @param node_regimes Regime assignments for all nodes
#' @param root_dists Root-to-node distances
#' @return List of regime path info for each leaf
build_regime_paths <- function(tree, node_regimes, root_dists) {
  n_tips <- length(tree$tip.label)
  root <- n_tips + 1
  
  # For each leaf, find the path from root to leaf and identify regime segments
  paths <- list()
  for (i in 1:n_tips) {
    # Get path from root to tip i
    path_nodes <- get_path_from_root(tree, i)
    
    # Identify regime segments
    segments <- list()
    current_regime <- node_regimes[path_nodes[1]]
    segment_start_dist <- root_dists[path_nodes[1]]
    
    for (k in 2:length(path_nodes)) {
      node <- path_nodes[k]
      if (node_regimes[node] != current_regime) {
        # End current segment
        segment_end_dist <- root_dists[path_nodes[k - 1]]
        segments[[length(segments) + 1]] <- list(
          regime = current_regime,
          start_dist = segment_start_dist,
          end_dist = segment_end_dist
        )
        # Start new segment
        current_regime <- node_regimes[node]
        segment_start_dist <- root_dists[path_nodes[k - 1]]
      }
    }
    # Final segment
    segments[[length(segments) + 1]] <- list(
      regime = current_regime,
      start_dist = segment_start_dist,
      end_dist = root_dists[path_nodes[length(path_nodes)]]
    )
    
    paths[[i]] <- segments
  }
  
  paths
}

#' Get the path from root to a given tip
#' 
#' @param tree phylo object
#' @param tip_index Index of the tip
#' @return Vector of node indices from root to tip
get_path_from_root <- function(tree, tip_index) {
  n_tips <- length(tree$tip.label)
  root <- n_tips + 1
  
  # Build parent lookup
  parent_of <- integer(n_tips + tree$Nnode)
  for (i in 1:nrow(tree$edge)) {
    parent_of[tree$edge[i, 2]] <- tree$edge[i, 1]
  }
  
  # Trace from tip to root
  path <- tip_index
  current <- tip_index
  while (current != root) {
    current <- parent_of[current]
    path <- c(current, path)
  }
  
  path
}

###############################################################################
# 2. OU PROCESS MODEL
###############################################################################

#' Compute the OU weight matrix W
#' 
#' W determines the expected latent state at each leaf as a weighted sum of
#' regime-specific optima theta.
#' Equation 10 in the paper.
#'
#' @param tree_info Output from parse_tree
#' @param alpha OU selection strength parameter
#' @return Weight matrix W (n_leaves x n_regimes)
compute_W_matrix <- function(tree_info, alpha) {
  n_leaves <- tree_info$n_leaves
  regimes <- tree_info$unique_regimes
  n_regimes <- length(regimes)
  
  W <- matrix(0, nrow = n_leaves, ncol = n_regimes)
  colnames(W) <- regimes
  
  for (i in 1:n_leaves) {
    segments <- tree_info$regime_paths[[i]]
    t_total <- tree_info$leaf_dists[i]
    
    for (seg in segments) {
      regime_col <- which(regimes == seg$regime)
      # Contribution from this regime segment
      # e^{-alpha*(t - t_tau)} - e^{-alpha*(t - t_{tau-1})}
      contrib <- exp(-alpha * (t_total - seg$end_dist)) - 
                 exp(-alpha * (t_total - seg$start_dist))
      W[i, regime_col] <- W[i, regime_col] + contrib
    }
  }
  
  W
}

#' Compute the OU covariance matrix Sigma
#' 
#' Equation 11 in the paper.
#' Sigma_ij = (sigma^2 / 2*alpha) * exp(-alpha*(t_i + t_j - 2*s_ij)) * (1 - exp(-2*alpha*s_ij))
#'
#' @param tree_info Output from parse_tree
#' @param alpha OU selection strength
#' @param sigma OU diffusion parameter
#' @return Covariance matrix Sigma (n_leaves x n_leaves)
compute_Sigma <- function(tree_info, alpha, sigma) {
  n <- tree_info$n_leaves
  ti <- tree_info$leaf_dists
  sij <- tree_info$shared_lengths
  
  Sigma <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in i:n) {
      val <- (sigma^2 / (2 * alpha)) * 
             exp(-alpha * (ti[i] + ti[j] - 2 * sij[i, j])) * 
             (1 - exp(-2 * alpha * sij[i, j]))
      Sigma[i, j] <- val
      Sigma[j, i] <- val
    }
  }
  
  # Add small jitter for numerical stability
  Sigma <- Sigma + diag(1e-6, n)
  
  Sigma
}

#' Compute the scaled covariance Sigma_tilde = (1/sigma^2) * Sigma
#' 
#' @param tree_info Output from parse_tree
#' @param alpha OU selection strength
#' @return Scaled covariance matrix
compute_Sigma_tilde <- function(tree_info, alpha) {
  n <- tree_info$n_leaves
  ti <- tree_info$leaf_dists
  sij <- tree_info$shared_lengths
  
  Sigma_tilde <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in i:n) {
      val <- (1 / (2 * alpha)) * 
             exp(-alpha * (ti[i] + ti[j] - 2 * sij[i, j])) * 
             (1 - exp(-2 * alpha * sij[i, j]))
      Sigma_tilde[i, j] <- val
      Sigma_tilde[j, i] <- val
    }
  }
  
  Sigma_tilde + diag(1e-6, n)
}

###############################################################################
# 3. OBSERVATION MODELS
###############################################################################

#' Softplus transformation: lambda = log(1 + exp(z))
#' Maps latent states to positive expression rates. Equation 4.
#'
#' @param z Latent state values
#' @return Positive expression rates lambda
softplus <- function(z) {
  # Numerically stable version
  ifelse(z > 20, z, log(1 + exp(z)))
}

#' Derivative of softplus: sigmoid function
#' @param z Latent state values
#' @return Derivatives
softplus_deriv <- function(z) {
  1 / (1 + exp(-z))
}

#' Log of softplus
#' @param z Latent state values
#' @return log(softplus(z))
log_softplus <- function(z) {
  log(softplus(z))
}

#' Negative Binomial log-likelihood (per cell)
#' 
#' Equation 13 in the paper.
#'
#' @param x Observed read counts
#' @param lambda Expected expression rates
#' @param r Dispersion parameter
#' @return Log-likelihood values
negbin_loglik <- function(x, lambda, r) {
  lgamma(x + r) - lgamma(r) - lgamma(x + 1) +
    x * log(lambda / (r + lambda)) +
    r * log(r / (r + lambda))
}

#' Poisson log-likelihood (per cell)
#' 
#' Equation 12 in the paper.
#'
#' @param x Observed read counts
#' @param lambda Expected expression rates
#' @return Log-likelihood values
poisson_loglik <- function(x, lambda) {
  x * log(lambda) - lambda - lgamma(x + 1)
}

###############################################################################
# 4. VARIATIONAL INFERENCE
###############################################################################

#' Compute expected values under q(z) using Monte Carlo sampling
#' 
#' Approximates E_q[softplus(z)], E_q[log(softplus(z))], 
#' and E_q[log(r + softplus(z))] using MC samples.
#'
#' @param mu_q Variational means (n_leaves)
#' @param sigma_q Variational standard deviations (n_leaves)
#' @param r Dispersion parameter
#' @param n_samples Number of MC samples
#' @return List with expected values
compute_mc_expectations <- function(mu_q, sigma_q, r, n_samples = 100) {
  n <- length(mu_q)
  
  # Sample from q(z) = N(mu_q, sigma_q^2)
  # z_samples: n_samples x n matrix
  z_samples <- matrix(rnorm(n * n_samples), nrow = n_samples, ncol = n)
  z_samples <- sweep(z_samples, 2, sigma_q, "*")
  z_samples <- sweep(z_samples, 2, mu_q, "+")
  
  # Compute softplus for all samples
  lambda_samples <- softplus(z_samples)
  
  # E[lambda]
  E_lambda <- colMeans(lambda_samples)
  
  # E[log(lambda)]
  E_log_lambda <- colMeans(log(lambda_samples + 1e-10))
  
  # E[log(r + lambda)]
  E_log_r_plus_lambda <- colMeans(log(r + lambda_samples))
  
  list(
    E_lambda = E_lambda,
    E_log_lambda = E_log_lambda,
    E_log_r_plus_lambda = E_log_r_plus_lambda
  )
}

#' Compute expected values using Taylor approximations (alternative to MC)
#' 
#' Equations 27-29 in the paper.
#'
#' @param mu_q Variational means
#' @param sigma_q Variational standard deviations
#' @param r Dispersion parameter
#' @return List with expected values
compute_taylor_expectations <- function(mu_q, sigma_q, r) {
  n <- length(mu_q)
  sigma_q2 <- sigma_q^2
  
  # E[softplus(z)] - Equation 27
  E_lambda <- ifelse(
    mu_q < 5,
    softplus(mu_q) + sigma_q2 / 2 * sigmoid(mu_q) * (1 - sigmoid(mu_q)),
    mu_q
  )
  
  # E[log(softplus(z))] - Equation 28
  E_log_lambda <- numeric(n)
  for (i in 1:n) {
    if (mu_q[i] < 2) {
      omega <- 1 / (1 + exp(2 * (mu_q[i] + 2)))
      T0 <- log(log(2)) + mu_q[i] / (2 * log(2)) + 
            (log(2) - 1) * (mu_q[i]^2 + sigma_q2[i]) / (8 * log(2)^2)
      E_log_lambda[i] <- (1 - omega) * T0 + omega * mu_q[i]
    } else if (mu_q[i] < 10) {
      E_log_lambda[i] <- log(mu_q[i] + exp(-mu_q[i])) - 
                          (1 - exp(-mu_q[i]))^2 / (2 * (mu_q[i] + exp(-mu_q[i]))^2) * sigma_q2[i]
    } else {
      E_log_lambda[i] <- log(mu_q[i])
    }
  }
  
  # E[log(r + softplus(z))] via MC (no closed-form Taylor given)
  # Fall back to MC for this term
  z_samples <- matrix(rnorm(n * 100), nrow = 100, ncol = n)
  z_samples <- sweep(z_samples, 2, sigma_q, "*")
  z_samples <- sweep(z_samples, 2, mu_q, "+")
  lambda_samples <- softplus(z_samples)
  E_log_r_plus_lambda <- colMeans(log(r + lambda_samples))
  
  list(
    E_lambda = E_lambda,
    E_log_lambda = E_log_lambda,
    E_log_r_plus_lambda = E_log_r_plus_lambda
  )
}

#' Sigmoid function
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

#' Compute the Evidence Lower Bound (ELBO)
#' 
#' Equation 18/26 in the paper:
#' L = E_q[log p(z)] + E_q[log p(x|lambda)] - E_q[log q(z)]
#'
#' @param x Observed read counts (n_leaves vector)
#' @param mu_q Variational means (n_leaves)
#' @param sigma_q Variational std devs (n_leaves)
#' @param tree_info Tree structure from parse_tree
#' @param alpha OU selection strength
#' @param sigma OU diffusion parameter (unused when using KKT)
#' @param theta Regime optima (n_regimes vector, or NULL to use KKT)
#' @param r NB dispersion parameter
#' @param beta L2 regularization strength for alpha
#' @param use_kkt Whether to use KKT conditions for theta and sigma
#' @param obs_model "negbin" or "poisson"
#' @param mc_samples Number of MC samples for expectations
#' @return ELBO value (scalar)
compute_ELBO <- function(x, mu_q, sigma_q, tree_info, alpha, sigma = NULL,
                          theta = NULL, r, beta = 1.0, use_kkt = TRUE,
                          obs_model = "negbin", mc_samples = 100) {
  
  n <- tree_info$n_leaves
  sigma_q2 <- sigma_q^2
  
  # --- OU prior term: E_q[log p(z)] ---
  W <- compute_W_matrix(tree_info, alpha)
  Sigma_tilde <- compute_Sigma_tilde(tree_info, alpha)
  
  if (use_kkt) {
    # KKT conditions: Equations 22-24
    Sigma_tilde_inv <- solve(Sigma_tilde)
    
    # Estimate theta_hat (Equation 22)
    theta_hat <- solve(t(W) %*% Sigma_tilde_inv %*% W) %*% 
                 t(W) %*% Sigma_tilde_inv %*% mu_q
    theta <- as.vector(theta_hat)
    
    # Estimate sigma_hat^2 (Equation 23)
    residual <- mu_q - W %*% theta
    sigma_hat2 <- (1/n) * (t(residual) %*% Sigma_tilde_inv %*% residual +
                            sum(diag(Sigma_tilde_inv) * sigma_q2))
    sigma_hat2 <- as.numeric(sigma_hat2)
    
    # Simplified E_q[log p(z)] using Equation 24
    log_det_Sigma_tilde <- determinant(Sigma_tilde, logarithm = TRUE)$modulus
    E_log_prior <- -0.5 * (as.numeric(log_det_Sigma_tilde) + n * log(sigma_hat2))
    
  } else {
    # Direct computation
    Sigma <- compute_Sigma(tree_info, alpha, sigma)
    Sigma_inv <- solve(Sigma)
    log_det_Sigma <- determinant(Sigma, logarithm = TRUE)$modulus
    
    mean_vec <- W %*% theta
    residual <- mu_q - mean_vec
    
    E_log_prior <- -0.5 * (as.numeric(log_det_Sigma) +
                            t(residual) %*% Sigma_inv %*% residual +
                            sum(diag(Sigma_inv) * sigma_q2))
    sigma_hat2 <- sigma^2
  }
  
  # --- Observation likelihood term: E_q[log p(x|lambda)] ---
  mc <- compute_mc_expectations(mu_q, sigma_q, r, n_samples = mc_samples)
  
  if (obs_model == "negbin") {
    # Equation 26
    E_log_obs <- sum(
      lgamma(x + r) - lgamma(r) + r * log(r) +
        x * mc$E_log_lambda -
        (x + r) * mc$E_log_r_plus_lambda
    )
  } else {
    # Poisson: Equation 25
    E_log_obs <- sum(x * mc$E_log_lambda - mc$E_lambda)
  }
  
  # --- Entropy of q(z): -E_q[log q(z)] = 0.5 * sum(log(sigma_q^2)) + const ---
  entropy <- 0.5 * sum(log(sigma_q2))
  
  # --- L2 regularization on alpha (Appendix A.2) ---
  reg <- -alpha^2 / (2 * beta)
  
  ELBO <- as.numeric(E_log_prior) + E_log_obs + entropy + reg
  
  attr(ELBO, "theta") <- theta
  attr(ELBO, "sigma2") <- sigma_hat2
  attr(ELBO, "E_log_prior") <- as.numeric(E_log_prior)
  attr(ELBO, "E_log_obs") <- E_log_obs
  attr(ELBO, "entropy") <- entropy
  
  ELBO
}

###############################################################################
# 5. OPTIMIZATION
###############################################################################

#' Fit VOUS model for a single gene
#' 
#' Optimizes variational parameters (mu_q, sigma_q) and model parameters
#' (alpha, r) to maximize the ELBO.
#'
#' @param x Observed read counts (n_leaves vector)
#' @param tree_info Tree structure from parse_tree
#' @param obs_model "negbin" or "poisson"
#' @param max_iter Maximum number of optimization iterations
#' @param lr Learning rate for gradient-based updates
#' @param tol Convergence tolerance (relative change in ELBO)
#' @param beta L2 regularization strength for alpha
#' @param init_alpha Initial alpha value
#' @param init_r Initial dispersion parameter
#' @param use_em Whether to use EM algorithm
#' @param mc_samples Number of MC samples
#' @param verbose Print progress
#' @return List with fitted parameters and ELBO trajectory
fit_vous_gene <- function(x, tree_info, obs_model = "negbin",
                           max_iter = 500, lr = 0.01, tol = 1e-5,
                           beta = 1.0, init_alpha = 1.0, init_r = 1.0,
                           use_em = FALSE, mc_samples = 100,
                           verbose = TRUE) {
  
  n <- tree_info$n_leaves
  
  # Initialize variational parameters
  # mu_q initialized from log(x + 1) as a reasonable starting point
  mu_q <- log(x + 1)
  log_sigma_q <- rep(log(0.5), n)  # log-scale for positivity
  
  # Initialize model parameters (log-scale for positivity)
  log_alpha <- log(init_alpha)
  log_r <- log(init_r)
  
  # ELBO tracking
  elbo_history <- numeric(max_iter)
  best_elbo <- -Inf
  best_params <- NULL
  
  if (use_em) {
    # EM Algorithm
    for (iter in 1:max_iter) {
      alpha <- exp(log_alpha)
      r <- exp(log_r)
      sigma_q <- exp(log_sigma_q)
      
      # ---- E-step: optimize variational parameters ----
      e_obj <- function(par) {
        mq <- par[1:n]
        lsq <- par[(n+1):(2*n)]
        sq <- exp(lsq)
        elbo <- compute_ELBO(x, mq, sq, tree_info, alpha, r = r,
                             beta = beta, obs_model = obs_model,
                             mc_samples = mc_samples)
        -elbo  # minimize negative ELBO
      }
      
      e_init <- c(mu_q, log_sigma_q)
      e_result <- optim(e_init, e_obj, method = "L-BFGS-B",
                        control = list(maxit = 50, factr = 1e7))
      mu_q <- e_result$par[1:n]
      log_sigma_q <- e_result$par[(n+1):(2*n)]
      sigma_q <- exp(log_sigma_q)
      
      # ---- M-step: optimize model parameters ----
      m_obj <- function(par) {
        la <- par[1]
        lr_val <- par[2]
        a <- exp(la)
        rv <- exp(lr_val)
        elbo <- compute_ELBO(x, mu_q, sigma_q, tree_info, a, r = rv,
                             beta = beta, obs_model = obs_model,
                             mc_samples = mc_samples)
        -elbo
      }
      
      m_init <- c(log_alpha, log_r)
      m_result <- optim(m_init, m_obj, method = "L-BFGS-B",
                        control = list(maxit = 50, factr = 1e7))
      log_alpha <- m_result$par[1]
      log_r <- m_result$par[2]
      
      # Compute ELBO
      alpha <- exp(log_alpha)
      r <- exp(log_r)
      sigma_q <- exp(log_sigma_q)
      current_elbo <- compute_ELBO(x, mu_q, sigma_q, tree_info, alpha, r = r,
                                    beta = beta, obs_model = obs_model,
                                    mc_samples = mc_samples)
      elbo_history[iter] <- current_elbo
      
      if (current_elbo > best_elbo) {
        best_elbo <- current_elbo
        best_params <- list(mu_q = mu_q, sigma_q = sigma_q,
                           alpha = alpha, r = r,
                           theta = attr(current_elbo, "theta"),
                           sigma2 = attr(current_elbo, "sigma2"))
      }
      
      if (verbose && iter %% 10 == 0) {
        cat(sprintf("EM Iter %d: ELBO = %.4f, alpha = %.4f, r = %.4f\n",
                    iter, current_elbo, alpha, r))
      }
      
      # Check convergence
      if (iter > 1) {
        rel_change <- abs((elbo_history[iter] - elbo_history[iter-1]) / 
                         (abs(elbo_history[iter-1]) + 1e-10))
        if (rel_change < tol) {
          if (verbose) cat(sprintf("Converged at iteration %d\n", iter))
          elbo_history <- elbo_history[1:iter]
          break
        }
      }
    }
    
  } else {
    # Joint optimization (default, as described in paper)
    joint_obj <- function(par) {
      mq <- par[1:n]
      lsq <- par[(n+1):(2*n)]
      la <- par[2*n + 1]
      lr_val <- par[2*n + 2]
      
      sq <- exp(lsq)
      a <- exp(la)
      rv <- exp(lr_val)
      
      elbo <- compute_ELBO(x, mq, sq, tree_info, a, r = rv,
                           beta = beta, obs_model = obs_model,
                           mc_samples = mc_samples)
      -elbo
    }
    
    init_par <- c(mu_q, log_sigma_q, log_alpha, log_r)
    
    # Use L-BFGS-B for optimization
    result <- optim(init_par, joint_obj, method = "L-BFGS-B",
                    control = list(maxit = max_iter, factr = 1e7, 
                                   trace = ifelse(verbose, 1, 0)))
    
    # Extract results
    mu_q <- result$par[1:n]
    sigma_q <- exp(result$par[(n+1):(2*n)])
    alpha <- exp(result$par[2*n + 1])
    r <- exp(result$par[2*n + 2])
    
    final_elbo <- compute_ELBO(x, mu_q, sigma_q, tree_info, alpha, r = r,
                                beta = beta, obs_model = obs_model,
                                mc_samples = mc_samples)
    
    best_params <- list(
      mu_q = mu_q,
      sigma_q = sigma_q,
      alpha = alpha,
      r = r,
      theta = attr(final_elbo, "theta"),
      sigma2 = attr(final_elbo, "sigma2")
    )
    best_elbo <- as.numeric(final_elbo)
    elbo_history <- best_elbo
  }
  
  list(
    params = best_params,
    elbo = best_elbo,
    elbo_history = elbo_history
  )
}

###############################################################################
# 6. HYPOTHESIS TESTING
###############################################################################

#' Fit null and alternative models and perform likelihood ratio test
#' 
#' Null model (H0): uniform regime (single theta across the tree)
#' Alternative model (H1): distinct regimes (different theta per tissue)
#' Equation 19 in the paper.
#'
#' @param x Observed read counts (n_leaves vector)
#' @param tree_info_null Tree info for null model (uniform regime)
#' @param tree_info_alt Tree info for alternative model (distinct regimes)
#' @param obs_model "negbin" or "poisson"
#' @param max_iter Maximum iterations
#' @param beta Regularization strength
#' @param mc_samples MC samples
#' @param verbose Print progress
#' @return List with LRT statistic, p-value, and fitted models
vous_lrt <- function(x, tree_info_null, tree_info_alt,
                      obs_model = "negbin", max_iter = 500,
                      beta = 1.0, mc_samples = 100, verbose = FALSE) {
  
  if (verbose) cat("Fitting null model...\n")
  fit_null <- fit_vous_gene(x, tree_info_null, obs_model = obs_model,
                             max_iter = max_iter, beta = beta,
                             mc_samples = mc_samples, verbose = verbose)
  
  if (verbose) cat("Fitting alternative model...\n")
  fit_alt <- fit_vous_gene(x, tree_info_alt, obs_model = obs_model,
                            max_iter = max_iter, beta = beta,
                            mc_samples = mc_samples, verbose = verbose)
  
  # Likelihood ratio statistic (Equation 19)
  lambda_LR <- 2 * (fit_alt$elbo - fit_null$elbo)
  lambda_LR <- max(0, lambda_LR)  # Ensure non-negative
  
  # Degrees of freedom = difference in number of theta parameters
  df <- length(tree_info_alt$unique_regimes) - length(tree_info_null$unique_regimes)
  
  # P-value from chi-squared distribution
  p_value <- pchisq(lambda_LR, df = df, lower.tail = FALSE)
  
  list(
    lambda_LR = lambda_LR,
    df = df,
    p_value = p_value,
    fit_null = fit_null,
    fit_alt = fit_alt,
    log2_fc = compute_log2fc(fit_alt$params$theta, tree_info_alt$unique_regimes)
  )
}

#' Compute log2 fold change between regimes using softplus-transformed theta
#' 
#' @param theta Estimated theta values
#' @param regimes Regime names
#' @return Log2 fold change (last regime vs first regime)
compute_log2fc <- function(theta, regimes) {
  if (length(theta) < 2) return(NA)
  lambda_theta <- softplus(theta)
  log2(lambda_theta[length(theta)] / lambda_theta[1])
}

#' Run VOUS on multiple genes with BH correction
#' 
#' @param count_matrix Gene x cell count matrix (genes in rows, cells in columns)
#' @param tree_info_null Null model tree info
#' @param tree_info_alt Alternative model tree info
#' @param gene_names Names of genes to test
#' @param obs_model Observation model
#' @param max_iter Maximum iterations per gene
#' @param beta Regularization strength
#' @param mc_samples MC samples
#' @param verbose Print progress
#' @return Data frame with results for each gene
vous_test_genes <- function(count_matrix, tree_info_null, tree_info_alt,
                             gene_names = NULL, obs_model = "negbin",
                             max_iter = 500, beta = 1.0, mc_samples = 100,
                             verbose = TRUE) {
  
  if (is.null(gene_names)) gene_names <- rownames(count_matrix)
  n_genes <- length(gene_names)
  
  results <- data.frame(
    gene = gene_names,
    lambda_LR = numeric(n_genes),
    p_value = numeric(n_genes),
    q_value = numeric(n_genes),
    log2fc = numeric(n_genes),
    alpha_null = numeric(n_genes),
    r_null = numeric(n_genes),
    alpha_alt = numeric(n_genes),
    r_alt = numeric(n_genes),
    converged = logical(n_genes),
    stringsAsFactors = FALSE
  )
  
  for (g in 1:n_genes) {
    if (verbose && g %% 10 == 0) {
      cat(sprintf("Processing gene %d / %d: %s\n", g, n_genes, gene_names[g]))
    }
    
    x <- as.numeric(count_matrix[gene_names[g], ])
    
    tryCatch({
      lrt_result <- vous_lrt(x, tree_info_null, tree_info_alt,
                              obs_model = obs_model, max_iter = max_iter,
                              beta = beta, mc_samples = mc_samples,
                              verbose = FALSE)
      
      results$lambda_LR[g] <- lrt_result$lambda_LR
      results$p_value[g] <- lrt_result$p_value
      results$log2fc[g] <- lrt_result$log2_fc
      results$alpha_null[g] <- lrt_result$fit_null$params$alpha
      results$r_null[g] <- lrt_result$fit_null$params$r
      results$alpha_alt[g] <- lrt_result$fit_alt$params$alpha
      results$r_alt[g] <- lrt_result$fit_alt$params$r
      results$converged[g] <- TRUE
    }, error = function(e) {
      if (verbose) cat(sprintf("  Gene %s failed: %s\n", gene_names[g], e$message))
      results$p_value[g] <<- NA
      results$converged[g] <<- FALSE
    })
  }
  
  # Benjamini-Hochberg FDR correction
  results$q_value <- p.adjust(results$p_value, method = "BH")
  
  results
}

###############################################################################
# 7. SIMULATION
###############################################################################

#' Simulate gene expression data using the VOUS generative model
#' 
#' Follows the full generative process: OU -> softplus -> NegBin
#'
#' @param tree_info Tree structure from parse_tree
#' @param alpha OU selection strength
#' @param sigma OU diffusion parameter
#' @param theta Regime optima (named vector matching regime names)
#' @param r Negative Binomial dispersion
#' @param n_genes Number of genes to simulate
#' @param obs_model "negbin" or "poisson"
#' @return List with simulated data
simulate_vous <- function(tree_info, alpha, sigma, theta, r, n_genes = 1,
                           obs_model = "negbin") {
  
  n <- tree_info$n_leaves
  n_regimes <- length(tree_info$unique_regimes)
  
  # Ensure theta is properly ordered
  theta_vec <- theta[tree_info$unique_regimes]
  
  # Compute OU mean and covariance
  W <- compute_W_matrix(tree_info, alpha)
  Sigma <- compute_Sigma(tree_info, alpha, sigma)
  
  mean_vec <- as.vector(W %*% theta_vec)
  
  # Cholesky decomposition for sampling
  L <- t(chol(Sigma))
  
  # Simulate latent states and counts for each gene
  z_matrix <- matrix(0, nrow = n_genes, ncol = n)
  lambda_matrix <- matrix(0, nrow = n_genes, ncol = n)
  x_matrix <- matrix(0, nrow = n_genes, ncol = n)
  
  for (g in 1:n_genes) {
    # Sample z from multivariate normal (OU prior)
    eps <- rnorm(n)
    z <- mean_vec + L %*% eps
    z_matrix[g, ] <- z
    
    # Softplus transformation
    lambda <- softplus(z)
    lambda_matrix[g, ] <- lambda
    
    # Sample counts
    if (obs_model == "negbin") {
      x_matrix[g, ] <- rnbinom(n, size = r, mu = lambda)
    } else {
      x_matrix[g, ] <- rpois(n, lambda = lambda)
    }
  }
  
  colnames(z_matrix) <- tree_info$tree$tip.label
  colnames(lambda_matrix) <- tree_info$tree$tip.label
  colnames(x_matrix) <- tree_info$tree$tip.label
  
  list(
    z = z_matrix,
    lambda = lambda_matrix,
    x = x_matrix,
    params = list(alpha = alpha, sigma = sigma, theta = theta, r = r)
  )
}

###############################################################################
# 8. ANCESTRAL STATE RECONSTRUCTION
###############################################################################

#' Reconstruct ancestral cell states on the lineage tree
#' 
#' Uses the fitted variational parameters and OU model to estimate
#' latent states at internal nodes.
#'
#' @param fit Fitted VOUS model (output from fit_vous_gene)
#' @param tree_info Tree structure
#' @return Vector of estimated latent states at all nodes
reconstruct_ancestral <- function(fit, tree_info) {
  params <- fit$params
  mu_q <- params$mu_q
  alpha <- params$alpha
  theta <- params$theta
  sigma2 <- params$sigma2
  
  tree <- tree_info$tree
  n_tips <- length(tree$tip.label)
  n_nodes <- tree$Nnode + n_tips
  root <- n_tips + 1
  
  # Initialize with leaf values
  states <- numeric(n_nodes)
  states[1:n_tips] <- mu_q
  
  # For internal nodes, use weighted average of descendant states
  # based on OU dynamics (simple approach)
  parent_of <- integer(n_nodes)
  children_of <- vector("list", n_nodes)
  edge_lengths <- numeric(n_nodes)
  
  for (i in 1:nrow(tree$edge)) {
    p <- tree$edge[i, 1]
    c <- tree$edge[i, 2]
    parent_of[c] <- p
    children_of[[p]] <- c(children_of[[p]], c)
    edge_lengths[c] <- tree$edge.length[i]
  }
  
  # Process nodes in post-order
  post_order <- postorder(tree)
  for (idx in seq_along(post_order)) {
    edge_idx <- post_order[idx]
    p <- tree$edge[edge_idx, 1]
    c <- tree$edge[edge_idx, 2]
  }
  
  # Simple reconstruction: for each internal node, average descendant
  # states weighted by branch length contribution
  node_order <- rev(unique(tree$edge[postorder(tree), 1]))
  
  for (node in node_order) {
    kids <- children_of[[node]]
    if (length(kids) > 0) {
      # Weight by inverse of branch length (closer children contribute more)
      weights <- exp(-alpha * edge_lengths[kids])
      weights <- weights / sum(weights)
      
      regime <- tree_info$node_regimes[as.character(node)]
      theta_node <- theta[which(tree_info$unique_regimes == regime)]
      
      # OU conditional: state is pulled toward theta
      child_states <- states[kids]
      # Weighted combination of children's states reverted through OU
      states[node] <- sum(weights * (child_states - theta_node) * 
                         exp(-alpha * edge_lengths[kids])) + theta_node
    }
  }
  
  names(states) <- c(tree$tip.label, paste0("Node", (n_tips+1):n_nodes))
  states
}

###############################################################################
# 9. VISUALIZATION UTILITIES
###############################################################################

#' Plot ELBO convergence
#' 
#' @param elbo_history Vector of ELBO values during optimization
#' @param gene_name Gene name for title
plot_elbo <- function(elbo_history, gene_name = "") {
  plot(seq_along(elbo_history), elbo_history, type = "l",
       xlab = "Iteration", ylab = "ELBO",
       main = paste("ELBO Convergence", gene_name),
       col = "steelblue", lwd = 2)
  grid()
}

#' Plot fitted expression on tree
#' 
#' @param tree_info Tree structure
#' @param fit Fitted model
#' @param x Observed counts
#' @param gene_name Gene name
plot_vous_fit <- function(tree_info, fit, x, gene_name = "") {
  par(mfrow = c(1, 2))
  
  # Plot tree with regime colors
  regime_colors <- rainbow(length(tree_info$unique_regimes))
  names(regime_colors) <- tree_info$unique_regimes
  tip_colors <- regime_colors[tree_info$regime_labels[tree_info$tree$tip.label]]
  
  plot(tree_info$tree, show.tip.label = FALSE, main = paste("Lineage Tree:", gene_name))
  tiplabels(pch = 16, col = tip_colors, cex = 0.5)
  legend("bottomleft", legend = tree_info$unique_regimes,
         col = regime_colors, pch = 16, cex = 0.8)
  
  # Plot observed vs fitted expression
  lambda_fit <- softplus(fit$params$mu_q)
  plot(x, lambda_fit, pch = 16, col = tip_colors, cex = 0.8,
       xlab = "Observed Counts", ylab = "Fitted Lambda",
       main = paste("Observed vs Fitted:", gene_name))
  abline(0, 1, lty = 2, col = "gray50")
  grid()
  
  par(mfrow = c(1, 1))
}

#' Volcano plot of VOUS results
#' 
#' @param results Data frame from vous_test_genes
#' @param q_threshold FDR threshold
#' @param fc_threshold Log2 fold change threshold
volcano_plot <- function(results, q_threshold = 0.05, fc_threshold = 1) {
  results <- results[results$converged, ]
  
  sig <- results$q_value < q_threshold & abs(results$log2fc) > fc_threshold
  
  cols <- ifelse(sig & results$log2fc > 0, "red",
          ifelse(sig & results$log2fc < 0, "blue", "gray70"))
  
  plot(results$log2fc, -log10(results$q_value),
       pch = 16, cex = 0.6, col = cols,
       xlab = expression(log[2] ~ "Fold Change"),
       ylab = expression(-log[10] ~ "q-value"),
       main = "VOUS Differential Expression")
  abline(h = -log10(q_threshold), lty = 2, col = "gray40")
  abline(v = c(-fc_threshold, fc_threshold), lty = 2, col = "gray40")
  
  # Label top genes
  top_genes <- results[sig, ]
  top_genes <- top_genes[order(top_genes$q_value), ]
  if (nrow(top_genes) > 10) top_genes <- top_genes[1:10, ]
  text(top_genes$log2fc, -log10(top_genes$q_value),
       labels = top_genes$gene, pos = 3, cex = 0.7)
  
  grid()
}

###############################################################################
# 10. COMPLETE EXAMPLE / DEMO
###############################################################################

#' Run a complete VOUS analysis demo
#' 
#' Creates a simulated tree, generates data, and tests for DE genes.
run_vous_demo <- function() {
  set.seed(42)
  cat("=== VOUS Demo ===\n\n")
  
  # --- Step 1: Create a simulated tree ---
  cat("Step 1: Simulating lineage tree...\n")
  n_tips <- 50
  tree <- rtree(n_tips)
  tree$edge.length <- tree$edge.length / max(node.depth.edgelength(tree))
  
  # Assign regimes: roughly half as "primary", half as "metastatic"
  # Use a subtree to define the metastatic regime
  tip_names <- tree$tip.label
  
  # Pick a random internal node and label its descendants as metastatic
  n_internal <- tree$Nnode
  subtree_node <- n_tips + sample(2:n_internal, 1)
  met_tips <- extract.clade(tree, subtree_node)$tip.label
  
  regime_labels <- setNames(
    ifelse(tip_names %in% met_tips, "metastatic", "primary"),
    tip_names
  )
  
  cat(sprintf("  Tree has %d leaves, %d primary, %d metastatic\n",
              n_tips, sum(regime_labels == "primary"), sum(regime_labels == "metastatic")))
  
  # --- Step 2: Parse tree for null and alternative models ---
  cat("Step 2: Parsing tree structures...\n")
  
  # Alternative model: two regimes
  tree_info_alt <- parse_tree(tree, regime_labels)
  
  # Null model: single regime
  null_labels <- setNames(rep("all", n_tips), tip_names)
  tree_info_null <- parse_tree(tree, null_labels)
  
  # --- Step 3: Simulate gene expression ---
  cat("Step 3: Simulating gene expression...\n")
  
  alpha_true <- 3
  sigma_true <- 1
  r_true <- 0.5
  
  # Simulate positive genes (DE) with theta shift in metastatic
  n_pos <- 20
  theta_pos <- c(primary = 0, metastatic = 1.5)
  sim_pos <- simulate_vous(tree_info_alt, alpha_true, sigma_true, theta_pos, r_true, n_pos)
  
  # Simulate negative genes (non-DE) with uniform theta
  n_neg <- 20
  theta_neg <- c(primary = 0, metastatic = 0)
  sim_neg <- simulate_vous(tree_info_alt, alpha_true, sigma_true, theta_neg, r_true, n_neg)
  
  # Combine into count matrix
  count_matrix <- rbind(sim_pos$x, sim_neg$x)
  gene_names <- c(paste0("DE_gene_", 1:n_pos), paste0("nonDE_gene_", 1:n_neg))
  rownames(count_matrix) <- gene_names
  
  cat(sprintf("  Simulated %d DE genes and %d non-DE genes\n", n_pos, n_neg))
  
  # --- Step 4: Run VOUS differential expression analysis ---
  cat("Step 4: Running VOUS...\n")
  
  results <- vous_test_genes(
    count_matrix, tree_info_null, tree_info_alt,
    gene_names = gene_names,
    obs_model = "negbin",
    max_iter = 200,
    beta = 1.0,
    mc_samples = 50,
    verbose = TRUE
  )
  
  # --- Step 5: Summarize results ---
  cat("\nStep 5: Results summary\n")
  cat(sprintf("  Total genes tested: %d\n", nrow(results)))
  cat(sprintf("  Converged: %d\n", sum(results$converged)))
  
  sig_genes <- results[results$q_value < 0.05 & !is.na(results$q_value), ]
  cat(sprintf("  Significant (q < 0.05): %d\n", nrow(sig_genes)))
  
  # Classification performance
  truth <- grepl("^DE_", results$gene)
  predicted <- results$q_value < 0.05 & !is.na(results$q_value)
  
  tp <- sum(truth & predicted)
  fp <- sum(!truth & predicted)
  fn <- sum(truth & !predicted)
  tn <- sum(!truth & !predicted)
  
  precision <- tp / max(tp + fp, 1)
  recall <- tp / max(tp + fn, 1)
  
  cat(sprintf("  True positives: %d, False positives: %d\n", tp, fp))
  cat(sprintf("  Precision: %.3f, Recall: %.3f\n", precision, recall))
  
  # Print top results
  cat("\nTop significant genes:\n")
  top <- head(results[order(results$p_value), ], 10)
  print(top[, c("gene", "lambda_LR", "p_value", "q_value", "log2fc")])
  
  # --- Step 6: Fit a single gene and show details ---
  cat("\nStep 6: Detailed fit for first DE gene...\n")
  x_example <- as.numeric(count_matrix["DE_gene_1", ])
  fit_example <- fit_vous_gene(x_example, tree_info_alt, obs_model = "negbin",
                                max_iter = 300, use_em = TRUE, verbose = TRUE)
  
  cat(sprintf("  Fitted alpha: %.4f\n", fit_example$params$alpha))
  cat(sprintf("  Fitted r: %.4f\n", fit_example$params$r))
  cat(sprintf("  Fitted theta: %s\n", 
              paste(sprintf("%.4f", fit_example$params$theta), collapse = ", ")))
  cat(sprintf("  Fitted sigma2: %.4f\n", fit_example$params$sigma2))
  cat(sprintf("  Final ELBO: %.4f\n", fit_example$elbo))
  
  # --- Optional: Ancestral reconstruction ---
  cat("\nStep 7: Ancestral state reconstruction...\n")
  ancestral <- reconstruct_ancestral(fit_example, tree_info_alt)
  cat(sprintf("  Root state estimate: %.4f\n", ancestral[n_tips + 1]))
  
  # Return everything
  invisible(list(
    tree_info_alt = tree_info_alt,
    tree_info_null = tree_info_null,
    count_matrix = count_matrix,
    results = results,
    fit_example = fit_example,
    ancestral = ancestral
  ))
}

###############################################################################
# RUN THE DEMO
###############################################################################
# Uncomment to run:
# demo_results <- run_vous_demo()