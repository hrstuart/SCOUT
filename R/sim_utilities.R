#############################################################

#Citation: Zhang, X., Xu, C. & Yosef, N. Simulating multiple faceted variability in single cell RNA sequencing. Nat Commun 10, 2611 (2019). https://doi.org/10.1038/s41467-019-10500-w
#https://github.com/YosefLab/SymSim/blob/master/R/simulation_functions.R

#############################################################
# Master Equation Related Functions
#############################################################

#' Getting GeneEffects matrices 
#'
#' This function randomly generates the effect size of each evf on the dynamic expression parameters
#' @param ngenes number of genes
#' @param nevf number of evfs
#' @param randomseed (should produce same result if ngenes, nevf and randseed are all the same)
#' @param prob the probability that the effect size is not 0
#' @param geffect_mean the mean of the normal distribution where the non-zero effect sizes are dropped from 
#' @param geffect_sd the standard deviation of the normal distribution where the non-zero effect sizes are dropped from 
#' @param evf_res the EVFs generated for cells
#' @param is_in_module a vector of length ngenes. 0 means the gene is not in any gene module. In case of non-zero values, genes with the same value in this vector are in the same module.
#' @return a list of 3 matrices, each of dimension ngenes * nevf
GeneEffects <- function(ngenes,nevf,randseed,prob,geffect_mean,geffect_sd, evf_res, is_in_module){
  set.seed(randseed)
  
  gene_effects <- lapply(c('kon','koff','s'),function(param){
    effect <- lapply(c(1:ngenes),function(i){
      nonzero <- sample(size=nevf,x=c(0,1),prob=c((1-prob),prob),replace=T)
      nonzero[nonzero!=0]=rnorm(sum(nonzero),mean=geffect_mean,sd=geffect_sd)
      return(nonzero)
    })
    return(do.call(rbind,effect))
  })
  
  
  mod_strength <- 0.8
  if (sum(is_in_module) > 0){ # some genes need to be assigned as module genes; their gene_effects for s will be updated accordingly
    for (pop in unique(is_in_module[is_in_module>0])){ # for every population, find the 
      #which(evf_res[[2]]$pop == pop)
      for (iparam in c(1,3)){
        nonzero_pos <- order(colMeans(evf_res[[1]][[iparam]][which(evf_res[[2]]$pop==pop),])- colMeans(evf_res[[1]][[iparam]]),
                             decreasing = T)[1:ceiling(nevf*prob)]
        geffects_centric <- rnorm(length(nonzero_pos), mean = 0.5, sd=0.5)
        for (igene in which(is_in_module==pop)){
          for (ipos in 1:length(nonzero_pos)){
            gene_effects[[iparam]][igene,nonzero_pos[ipos]] <- rnorm(1, mean=geffects_centric[ipos], sd=(1-mod_strength))
          }
          gene_effects[[iparam]][igene,setdiff((1:nevf),nonzero_pos)] <- 0
        }
      }
      
      nonzero_pos <- order(colMeans(evf_res[[1]][[2]][which(evf_res[[2]]$pop==pop),])- colMeans(evf_res[[1]][[2]]),
                           decreasing = T)[1:ceiling(nevf*prob)]
      geffects_centric <- rnorm(length(nonzero_pos), mean = -0.5, sd=0.5)
      for (igene in which(is_in_module==pop)){
        for (ipos in 1:length(nonzero_pos)){
          gene_effects[[2]][igene,nonzero_pos[ipos]] <- rnorm(1, mean=geffects_centric[ipos], sd=(1-mod_strength))
        }
        gene_effects[[2]][igene,setdiff((1:nevf),nonzero_pos)] <- 0
      }
    }
  }
  return(gene_effects)
}

#' sample from smoothed density function
#' @param nsample number of samples needed
#' @param den_fun density function estimated from density() from R default
SampleDen <- function(nsample,den_fun){
  probs <- den_fun$y/sum(den_fun$y)
  bw <- den_fun$x[2]-den_fun$x[1]
  bin_id <- sample(size=nsample,x=c(1:length(probs)),prob=probs,replace=T)
  counts <- table(bin_id)
  sampled_bins <- as.numeric(names(counts))
  samples <- lapply(c(1:length(counts)),function(j){
    runif(n=counts[j],min=(den_fun$x[sampled_bins[j]]-0.5*bw),max=(den_fun$x[sampled_bins[j]]+0.5*bw))
  })
  samples <- do.call(c,samples)
  return(samples)
}

#' Getting the parameters for simulating gene expression from EVf and gene effects
#'
#' This function takes gene_effect and EVF, take their dot product and scale the product to the correct range 
#' by using first a logistic function and then adding/dividing by constant to their correct range
#' @param gene_effects a list of three matrices (generated using the GeneEffects function), 
#' each corresponding to one kinetic parameter. Each matrix has nevf columns, and ngenes rows. 
#' @param evf a vector of length nevf, the cell specific extrinsic variation factor
#' @param match_param_den the fitted parameter distribution density to sample from 
#' @param bimod the bimodality constant
#' @param scale_s a factor to scale the s parameter, which is used to tune the size of the actual cell (small cells have less number of transcripts in total)
#' @return params a matrix of ngenes * 3
#' @examples 
#' Get_params()
Get_params.OLD <- function(gene_effects,evf,match_param_den,bimod,scale_s){
  params <- lapply(1:3, function(iparam){evf[[iparam]] %*% t(gene_effects[[iparam]])})
  scaled_params <- lapply(c(1:3),function(i){
    X <- params[[i]]
    # X=matrix(data=c(1:10),ncol=2) 
    # this line is to check that the row and columns did not flip
    temp <- plyr::alply(X, 1, function(Y){Y})
    values <- do.call(c,temp)
    ranks <- rank(values)
    sorted <- sort(SampleDen(nsample=max(ranks),den_fun=match_param_den[[i]]))
    temp3 <- matrix(data=sorted[ranks],ncol=length(X[1,]),byrow=T)
    return(temp3)
  })
  
  bimod_perc <- 1
  ngenes <- dim(scaled_params[[1]])[2]; bimod_vec <- numeric(ngenes)
  bimod_vec[1:ceiling(ngenes*bimod_perc)] <- bimod
  bimod_vec <- c(rep(bimod, ngenes/2), rep(0, ngenes/2))
  scaled_params[[1]] <- apply(t(scaled_params[[1]]),2,function(x){x <- 10^(x - bimod_vec)})
  scaled_params[[2]] <- apply(t(scaled_params[[2]]),2,function(x){x <- 10^(x - bimod_vec)})
  scaled_params[[3]] <- t(apply(scaled_params[[3]],2,function(x){x<-10^x}))*scale_s
  
  return(scaled_params)
}

#' Kinetic parameters from EVFs and gene effects.
#'
#' @param param_mode 'rank' (original SymSim behaviour) or 'affine'.
#' @param signal_gain length-3 over c(kon, koff, s): per-gene dynamic range in log10 units. NA leaves
#'   that parameter on the rank transform. Only used when param_mode = 'affine'.
#' @param burst_scale multiplies kon and koff. Beta mean kon/(kon+koff) is unchanged, variance falls
#'   ~1/burst_scale. 1 is the calibrated SymSim setting.
#' @param support_clip clamp the affine-mapped log10 parameters to the same support the rank path
#'   draws from, i.e. range(match_param_den[[i]]$x). Only used when param_mode = 'affine'; the rank
#'   path cannot leave that support in the first place, so this is inert there. Without it the
#'   affine map is unbounded and a single gene-cell count can reach ~1e6, far outside any real
#'   scRNA-seq range. The bound on the s parameter is 10^4.043, so the count ceiling is
#'   10^4.043 * scale_s.
#' @return list of three ngenes x ncells matrices: kon, koff, s.
Get_params <- function(gene_effects, evf, match_param_den, bimod, scale_s,
                          param_mode = c('rank', 'affine'),
                          signal_gain = c(NA, NA, 0.3),
                          burst_scale = 1,
                          support_clip = TRUE) {

  param_mode <- match.arg(param_mode)
  if (length(signal_gain) == 1) signal_gain <- rep(signal_gain, 3)
  if (length(signal_gain) != 3) stop('signal_gain must be length 1 or 3 (kon, koff, s).')

  params <- lapply(1:3, function(iparam) { evf[[iparam]] %*% t(gene_effects[[iparam]]) })

  scaled_params <- lapply(c(1:3), function(i) {
    X <- params[[i]]                                   # ncells x ngenes

    if (!(param_mode == 'affine' && !is.na(signal_gain[i]))) {
      # ---- original SymSim path, verbatim (same plyr round-trip, same SampleDen RNG draw) ----
      temp <- plyr::alply(X, 1, function(Y){Y})
      values <- do.call(c, temp)
      ranks <- rank(values)
      sorted <- sort(SampleDen(nsample = max(ranks), den_fun = match_param_den[[i]]))
      temp3 <- matrix(data = sorted[ranks], ncol = length(X[1, ]), byrow = T)
      return(temp3)
    }

    # ---- affine path: preserve the latent trait's spread ----
    ngenes  <- ncol(X)
    Xc      <- sweep(X, 2, colMeans(X), '-')
    # ONE pooled scale, not per-gene: genes genuinely differ in effect size (GeneEffects weights sum
    # to mean 0, sd 2.2, so ~1 gene in 6 has a cancelled signal) and per-gene standardisation would
    # flatten that real heterogeneity away.
    sd_pool <- stats::median(apply(X, 2, stats::sd))
    if (!is.finite(sd_pool) || sd_pool == 0) sd_pool <- 1

    # Realistic per-gene baseline levels, assigned in the order the gene means already imply, so
    # between-gene expression heterogeneity is unchanged and only the within-gene range moves.
    mu_g <- sort(SampleDen(nsample = ngenes, den_fun = match_param_den[[i]]))[
              rank(colMeans(X), ties.method = 'first')]

    out <- sweep(Xc / sd_pool * signal_gain[i], 2, mu_g, '+')

    # Clamp to the rank path's own support, so the two modes differ ONLY in how the latent trait is
    # mapped inside that support and not in how far outside it they are allowed to go. SampleDen
    # draws from match_param_den[[i]]$x, so this range is exactly what rank can produce.
    if (support_clip) {
      lim <- range(match_param_den[[i]]$x)
      out[] <- pmin(pmax(out, lim[1]), lim[2])
    }
    out
  })

  bimod_perc <- 1
  ngenes <- dim(scaled_params[[1]])[2]; bimod_vec <- numeric(ngenes)
  bimod_vec[1:ceiling(ngenes*bimod_perc)] <- bimod
  bimod_vec <- c(rep(bimod, ngenes/2), rep(0, ngenes/2))
  scaled_params[[1]] <- apply(t(scaled_params[[1]]),2,function(x){x <- 10^(x - bimod_vec)})
  scaled_params[[2]] <- apply(t(scaled_params[[2]]),2,function(x){x <- 10^(x - bimod_vec)})
  scaled_params[[3]] <- t(apply(scaled_params[[3]],2,function(x){x<-10^x}))*scale_s

  if (burst_scale != 1) {
    scaled_params[[1]] <- scaled_params[[1]] * burst_scale
    scaled_params[[2]] <- scaled_params[[2]] * burst_scale
  }

  scaled_params
}

#' Getting the parameters for simulating gene expression from EVf and gene effects
#'
#' This function takes gene_effect and EVF, take their dot product and scale the product to the correct range 
#' by using first a logistic function and then adding/dividing by constant to their correct range
#' @param evf a vector of length nevf, the cell specific extrinsic variation factor
#' @param gene_effects a list of three matrices (generated using the GeneEffects function), 
#' each corresponding to one kinetic parameter. Each matrix has nevf columns, and ngenes rows. 
#' @param param_realdata the fitted parameter distribution to sample from 
#' @param bimod the bimodality constant
#' @return params a matrix of ngenes * 3
#' @examples 
#' Get_params()
Get_params2 <- function(gene_effects,evf,bimod,ranges){
  params <- lapply(gene_effects,function(X){evf %*% t(X)})
  scaled_params <- lapply(c(1:3),function(i){
    X <- params[[i]]
    temp <- apply(X,2,function(x){1/(1+exp(-x))})
    temp2 <- temp*(ranges[[i]][2]-ranges[[i]][1])+ranges[[i]][1]
    return(temp2)
  })
  scaled_params[[1]]<-apply(scaled_params[[1]],2,function(x){x <- 10^(x - bimod)})
  scaled_params[[2]]<-apply(scaled_params[[2]],2,function(x){x <- 10^(x - bimod)})
  scaled_params[[3]]<-apply(scaled_params[[3]],2,function(x){x<-abs(x)})
  scaled_params <- lapply(scaled_params,t)
  return(scaled_params)
}

get_prob <- function(glength){
  if (glength >= 1000){prob <- 0.7} else{
    if (glength >= 100 & glength < 1000){prob <- 0.78}
    else if (glength < 100) {prob <- 0}
  }
  return(prob)
}