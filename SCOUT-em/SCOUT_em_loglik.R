##### likelihood functions 
invert_LL <- function(phy, edges, Rate.mat, root.state, y, W) {
# Function to calculate log liklihood using the invert algorithm. 
# Removed the tip.fog requirements. 
	N <- length(y)
   # print(phy)
   # print(head(edges))
   # print(Rate.mat)
    V <- varcov.ou(phy, edges, Rate.mat, 
        root.state=root.state, 
        simmap.tree=FALSE, 
        root.age=NULL, 
        scaleHeight=FALSE, 
        assume.station=FALSE, 
        shift.point=0.5)
   # print(V[0:10, 0:10])
    if (any(is.nan(diag(V))) || any(is.infinite(diag(V)))) {
      #  print('NA in diag of V. Returning max.')
        return(1000000)
    }

    theta <- Inf
    try(theta <- pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)%*%t(W)%*%pseudoinverse(V)%*%y, silent=TRUE)

    if(any(theta==Inf)){
        return(10000000)
    }
    DET <- sum(log(abs(Re(diag(qr(V)$qr)))))

    #However, sometimes this fails (not sure yet why) so I just toggle between this and another approach:
    if(!is.finite(DET)){
        DET <- determinant(V, logarithm=TRUE)
        logl <- -.5*(t(y-W%*%theta)%*%pseudoinverse(V)%*%(y-W%*%theta))-.5*as.numeric(DET$modulus)-.5*(N*log(2*pi))
    }else{
        logl <- -.5*(t(y-W%*%theta)%*%pseudoinverse(V)%*%(y-W%*%theta))-.5*as.numeric(DET)-.5*(N*log(2*pi))
    }

    return(-logl)
}

# NEED TO DOUBLE CHECK 1. get.root.theta as an option and 2. whether tip.fog gets set at some point. 
tp_LL <- function(phy, edges, Rate.mat, y, W){
	pars <- matrix(c(Rate.mat[3,], Rate.mat[2,], Rate.mat[1,]), dim(Rate.mat)[2], 3, 
                dimnames = list(levels(factor((tot.states))), c("opt", "sig", "alp")))
    if(get.root.theta == TRUE){
        root.par.index <- length(p)
        theta0 <- p[root.par.index]
        expected.vals <- colSums(t(W) * c(theta0, pars[,1]))
        names(expected.vals) <- phy$tip.label
    }else{
        expected.vals <- colSums(t(W) * pars[,1])
        names(expected.vals) <- phy$tip.label
    }
              
    comp <- NA
    try(comp <- phylolm::three.point.compute(transformed.tree$tree, y, expected.vals, transformed.tree$diag), silent=TRUE)
    if(is.na(comp[1])){
        return(10000000)
    }else{
        logl <- -as.numeric(Ntip(phy) * log(2 * pi) + comp$logd + (comp$PP - 2 * comp$QP + comp$QQ))/2
    }

    return (-logl)
}


loglikPois <- function(par) {
    lambda <- exp(par[1])
          
    ll <- sum(y * log(lambda) - lambda - lgamma(y + 1L))
          
    return(-ll)
}
loglikNB <- function(par){
    mu <- exp(par[1])
    k <- exp(par[2])

    ll <- sum(lgamma(y + k) - lgamma(k) - lgamma(y + 1L) + k * log(k) + y * log(mu) - (y + k ) + log(mu + k))

    return(-ll)
}

loglikZINB <- function(par) {
  # par[1] = log(mu) ; par[2] = log(k) ; par[3] = logit(pi)
  mu   <- exp(par[1])
  k    <- exp(par[2])
  pi   <- plogis(par[3])               # logistic transform ⇒ 0 < π < 1

  # Log‑likelihood split into zeros and positives
  is0   <- (y == 0)
  ispos <- !is0

  # Zero part
  logp0  <- log( pi + (1 - pi) * (k / (k + mu))^k )

  # Positive part (same NB kernel as before)
  logp_pos <- lgamma(y[ispos] + k) - lgamma(k) - lgamma(y[ispos] + 1) +
              k * log(k) + y[ispos] * log(mu) -
              (y[ispos] + k) * log(mu + k) + log(1 - pi)

  ll <- sum(logp0[is0]) + sum(logp_pos)
  return(-ll) 
}

loglikBP <- function(...){}
# Need to suss this out further. But would be great to be able to estimate kon, koff, and s and somehow use that to translate back to log space? 
# just kidding idk if that makes any sense at all... 
# and I have no idea how we do this without doing to log approach. I am not sure if I see the problem with using gaussian. 


################# Not doing the analytic version.... #################

pois2 <- function(y, lambda){
  ll <- sum(dpois(y, lambda, log=TRUE))
  return(-ll)
}

nb2 <- function(y, size, mu){
  ll <- sum(dnbinom(y, size = size, mu = mu, log = TRUE))
  return(-ll)
}

gaus <- function(y, mu, sd){
    ll <- sum(dnorm(y, mean = mu, sd = sd))
    return(-ll)
}

optim_joint <- function(p, obs, tree, edges, index.mat, Rate.mat, root.age, root.state, get.root.theta, version = NULL, mode = 'ignore') {
    p <- exp(p)
    np <- length(p)

    ix.max <- max(index.mat, na.rm=TRUE)

    if (any(is.na(index.mat))){
        index.mat[is.na(index.mat)] <- np + 1  # fill NA with 1e-10, fills in below. 
    }

    Rate.mat[] <- c(p, 1e-10)[index.mat]
    #print(Rate.mat)
    if(get.root.theta == TRUE){
        W <- weight.mat(tree, edges, Rate.mat, root.state=root.state, simmap.tree=FALSE, root.age=root.age, scaleHeight=FALSE, assume.station=FALSE, shift.point=0.5)
    }else{
        W <- weight.mat(tree, edges, Rate.mat, root.state=root.state, simmap.tree=FALSE, root.age=root.age, scaleHeight=FALSE, assume.station=TRUE, shift.point=0.5)
    }
    
    # Version with Gaussian noise, non-generative. 

    if (is.null(version)){
        obs_norm <- obs
        ll2 <- 0     

    } else if (version == 'NORM') {
        mu <- p[ix.max+1]; sd <- p[ix.max+2]
        if (mode == 'gen') {
            norm <- rnorm(length(obs), mean = mu, sd = sd)
            names(obs) <- tree$tip.label

            obs_sorted <- sort(obs)
            norm_sorted <- sort(norm)

            names(norm_sorted) <- names(obs_sorted)
            obs_norm <- norm_sorted[tree$tip.label]
        } else {
            obs_norm <- obs
        }

        #ll1 <- invert_LL(tree, edges, Rate.mat, root.state, obs_norm, W)
        # Calculate gaussian likelihood
        ll2 <- gaus(obs, mu = mu, sd = sd)
        
    } else if (version == 'NB') {
        mu <- p[ix.max+1]; size <- p[ix.max+2]
        
        if (mode == 'gen') {
            gen <- rnbinom(length(obs), mu = mu, size = size)
            names(obs) <- tree$tip.label
            obs_sorted <- sort(obs)
            gen_sorted <- sort(gen)

            names(gen_sorted) <- names(obs_sorted)
            obs_l <- gen_sorted[tree$tip.label]

            obs_norm <- log1p(obs_l)
        } else {
            obs_norm <- log1p(obs)
        }
        
        #ll1 <- invert_LL(tree, edges, Rate.mat, obs_norm, W)
        ll2 <- nb2(obs, size = size, mu = mu)
 
    } else if (version == 'POIS') {
        lambda <- p[ix.max+1]
        
        if (mode == 'gen') {
            gen <- rpois(length(obs), lambda = lambda)
            names(obs) <- tree$tip.label
            obs_sorted <- sort(obs)
            gen_sorted <- sort(gen)

            names(gen_sorted) <- names(obs_sorted)
            obs_l <- gen_sorted[tree$tip.label]

            obs_norm <- log1p(obs_l)
        } else {
            obs_norm <- log1p(obs)
        }
        
        #ll1 <- invert_LL(tree, edges, Rate.mat, obs_norm, W)
        ll2 <- pois2(obs, lambda = lambda)
        
    }

    # establish obs norm across the board. 
    # Calculate OU-likelihood
    ll1 <- invert_LL(tree, edges, Rate.mat, root.state, obs_norm, W)
  
   # cat('\nOU-likelihood:\n')
   # print(ll1)

    #cat('\nNoise likelihood:\n')
    #print(ll2)
    # if nothing special then ll2 should = 0. 
    if (any(is.na(c(ll1, ll2)))){return(10000000)}
    totll <- ll1 + ll2

    return(totll)
}

