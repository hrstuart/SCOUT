#' Calculate fit metrics and weights from model summaries
#' @param ngenes Integer describing number of genes to simuate for the tedsim pipeline. 
#' @param ncells Integer describing number of cells to simulate.
#' @param tree Newick string describing cell state tree. Controls how many final states are simulated.
#' @param outdir Path to output directory
#' @param out_prefix Prefix for saving files. 
#' @import TedSim
#' @import ape
#' @return List of results. The same set are also saved to outdir. 
#' @export
simulate_lineage <- function(ngenes, ncells, tree, outdir, out_prefix){
    max_walk <- 6
    p_d <- 0
    n_cif <- 30
    n_diff <- 20
    
    N_nodes <- 2*ncells-2
    N_char <- 64
    mu<-0.1
    p_a <- 0.4
    cif_step <- 0.4
    phyla <- ape::read.tree(text=tree)
    
    create_directory_if_not_exists(outdir)
    
    returnlist <- SIFGenerate(phyla, n_diff,step = cif_step)
     
    cifs <- SimulateCIFs(ncells,phyla ,p_a = p_a, n_CIF = n_cif, n_diff = n_diff,
                     step = cif_step, p_d = p_d, 
                     Sigma = 0.5, N_char = N_char, max_walk = max_walk, 
                     SIF_res = returnlist, unif_on = FALSE)
    
    cif_leaves <- lapply(c(1:3),function(parami){
    cif_leaves_all <- cifs[[1]][[parami]][c(1:ncells),]
          return(cif_leaves_all)
        })

    cif_res <- list(cif_leaves,cifs[[2]])
    states <- cifs[[2]]
    states <- states[1:N_nodes,]
    states_leaves <- states[1:ncells,]
    muts <- cifs[[7]]
    rownames(muts) <- paste("cell",states[,4],sep = "_")
    muts_leaves <- muts[1:ncells,]
    
    true_counts_res <- CIF2Truecounts(ngenes = ngenes, ncif = n_cif, ge_prob = 0.3, ncells = ncells, cif_res = cif_res)
    
    data(gene_len_pool)
    gene_len <- sample(gene_len_pool, ngenes, replace = FALSE)
    observed_counts <- True2ObservedCounts(true_counts=true_counts_res[[1]], meta_cell=true_counts_res[[3]], protocol="UMI", alpha_mean=0.2, alpha_sd=0.05, gene_len=gene_len, depth_mean=1e5, depth_sd=3e3)
    
    make_path <- paste0(outdir, '/', out_prefix)
    gene_expression_dir <- paste0(make_path, '_counts.csv')
    cell_meta_dir <- paste0(make_path, '_cell_meta.csv')
    character_matrix_dir <-paste0(make_path, '_character_matrix.csv')
    tree_gt_dir <- paste0(make_path, '_newick.nwk')
    internal_dir <- paste0(make_path, '_internal_node_annot.csv')

    states_leaves <- as.data.frame(states_leaves)
    states_leaves$nodeID <- paste0('t', states_leaves$cellID)

    ape::write.tree(cifs[[4]], tree_gt_dir)
    write.csv(observed_counts[[1]], gene_expression_dir, row.names = FALSE)
    write.csv(states_leaves, cell_meta_dir)
    write.table(muts_leaves, character_matrix_dir)
    write.csv(cifs[[2]], internal_dir)

    return(list(tree = cifs[[4]], observed_counts = observed_counts[[1]], 
                meta = states_leaves, states=cifs[[2]]))
}

#' Utility function to generate a vector of thetas. 
generate_thetas <- function(x, nstats, step) {
  # Calculate starting point so that vector is centered around x
  start <- x - step * floor((nstats - 1) / 2)
  if (start < 0) {
  	stop('Adjust theta parameters. The floor to simulate the set of theta values for OUM is less than 0.')
  }
  seq(from = start, by = step, length.out = nstats)
}

#'@import OUwie
generate_EvoEVFs <- function(real_tree, metadata, nevfs, a, s, t0, theta_step, n_states){
	#oum_thetas <- generate_thetas(t0, n_states, theta_step)

    bm_evf <- lapply(1:3, function(iparam){
        evfs <- sapply(1:nevfs, function(i){
            d <- OUwie.sim.edited(real_tree, metadata,  alpha=rep(1e-10, n_states), sigma.sq=rep(s, n_states), theta0=t0, theta=rep(t0, n_states), simmap.tree = FALSE)
            return(d$X)
        })
    return(evfs)
    })

    ou1_evf <- lapply(1:3, function(iparam){
        evfs <- sapply(1:nevfs, function(i){
            d <- OUwie.sim.edited(real_tree, metadata, alpha=rep(a, n_states), sigma.sq=rep(s, n_states), theta0=t0, theta=rep(t0+theta_step, n_states), simmap.tree = FALSE)
            return(d$X)
        })
        return(evfs)
    })

    #thetas <- generate_thetas(t0, n_states, theta_step)
    thetas <- seq(from = t0, by=theta_step, length.out = n_states)
    oum_evf <- lapply(1:3, function(iparam){
        evfs <- sapply(1:nevfs, function(i){
            d <- OUwie.sim.edited(real_tree, metadata,  alpha=rep(a, n_states), sigma.sq=rep(s, n_states), theta0=t0, theta=thetas, simmap.tree = FALSE)
            return(d$X)
        })
        return(evfs)
    })

    eevfs <- list(bm_evf, ou1_evf, oum_evf)
    
    return(eevfs)
}

generate_EvoCounts <- function(ncells, metadata, eevfs, nevfs,  ngenes){
    ecounts <- lapply(eevfs, function(imod){
        gene_effects <- GeneEffects(ngenes = ngenes, nevf = nevfs, randseed = 123, prob = 0.3, 
                                    geffect_mean = 0, geffect_sd = 1, is_in_module=0)

        data(match_params)
        match_params[,1]=log(base=10,match_params[,1])
                match_params[,2]=log(base=10,match_params[,2])
                match_params[,3]=log(base=10,match_params[,3])
                match_params_den <- lapply(c(1:3),function(i){
                    density(match_params[,i],n=2000)
                })

        params <- Get_params(gene_effects,imod,match_params_den,bimod=0,scale_s=10)

        counts <- lapply(c(1:ngenes),function(i){
            count <- sapply(c(1:ncells),function(j){
                y <- rbeta(1,params[[1]][i,j],params[[2]][i,j])
                x <- rpois(1,y*params[[3]][i,j])
                return(x)
            })
        })

        counts <- do.call(rbind, counts)
        return(counts)
    })
}

simulate_OU_genes <- function(tree_, metadata, a, s, t0, theta_step, nstates, nevfs, ngenes, ncells){ 
    # real_tree, metadata, nevfs, a, s, ngenes, n_states

    eevfs <- generate_EvoEVFs(tree_, metadata, nevfs, a, s, t0, theta_step, nstates)
    res <- generate_EvoCounts(ncells, metadata, eevfs, nevfs,  ngenes)

    #### GENERATE COUNTS DATA 
    counts <- do.call(rbind, res)
    counts <- t(counts)

    states <- metadata$cluster; names(states) <- metadata$cellID

    #### CLEAN FINAL DATASETS ####
    counts_df <- as.data.frame(counts)
    colnames(counts_df) <- c(paste0('BM1_', 1:ngenes), paste0('OU1_', 1:ngenes), paste0('OUM_', 1:ngenes))
    counts_df$OUM <- states[tree_$tip.label]
    counts_df$species <- tree_$tip.label

    ou_evfs <- cbind( eevfs[[1]][[1]], eevfs[[2]][[1]], eevfs[[3]][[1]] ) 
    ou_evf_df <- as.data.frame(ou_evfs)
    colnames(ou_evf_df) <- c(paste0('BM1_', 1:nevfs), paste0('OU1_', 1:nevfs), paste0('OUM_', 1:nevfs))
    ou_evf_df$OUM <- states[tree_$tip.label]
    ou_evf_df$species <- tree_$tip.label

    return(list(counts = counts_df, evfs = ou_evf_df))
}


#' Calculate fit metrics and weights from model summaries
#' @param ngenes Integer describing number of genes to simuate for the tedsim pipeline. 
#' @param ncells Integer describing number of cells to simulate. Note ncells < ngenes. 
#' @param tree Newick string describing cell state tree. Controls how many final states are simulated.
#' @param outdir Path to output directory
#' @param out_prefix Prefix for saving files. 
#' @param a Integer or list: Alpha to simulate values. 
#' @param s Integer of list: Sigma to simulate values. 
#' @param t0 Integer: Root theta values. 
#' @param theta_step Integer: Used to generate optimal values depending on the model. 
#' @param nevfs Integer: Number of EVFs to generate. 
#' @return List of results. The same set are also saved to outdir. 
#' @export
simulate_test_data <- function(ngenes, ncells, tree, outdir, out_prefix, a, s, t0, theta_step, nevfs = 20) {
	log_message(sprintf('Generating scLT data for %d cells.', ncells), verbose =TRUE)
	log_message(sprintf('Files will be saved to %s.', outdir), verbose =TRUE)

	test_dataset <- simulate_lineage(250, ncells, tree, outdir, out_prefix) # hardcoding this because we don't actually use these genes. 
    #states <- test_dataset[[3]][,2]; names(states) <- paste0('t',test_dataset[[3]][,4])
    metadata_f <- data.frame(test_dataset[[3]])
    metadata_f$cellID <- paste0('t', metadata_f$cellID)

    metadata <- metadata_f[, c('cellID', 'cluster')]
    metadata$cluster <- as.factor( metadata$cluster)
    
    #tree_$edge
    full_annot <- as.data.frame(test_dataset[[4]])
    full_annot <- full_annot[, c('cellID', 'cluster')]
    full_annot$cluster <- as.factor( full_annot$cluster)
    
    tree_ <- test_dataset[[1]]

    value_lookup <- setNames(full_annot$cluster, full_annot$cellID)
    for (nid in unique(tree_$edge[, 1])) {
        if (nid %in% names(value_lookup)) {
            tree_$node.label[nid] <- as.character(value_lookup[[as.character(nid)]])
        }
      }
    lower.bound <- ncells+1
    tree_$node.label <- tree_$node.label[lower.bound:nrow(full_annot)]
    
    nstates <- length(unique(metadata$cluster))
    combos = expand.grid(a, s)
    simulated_res = list('counts' = list(), 'evfs' = list())
    log_message(sprintf('Simulating OU genes for %d parameter combinations.', nrow(combos)), verbose =TRUE)
	
	for (i in 1:nrow(combos)){
		a_i <- combos[i, 1]
		s_i <- combos[i, 2]

		OU_sim <- simulate_OU_genes(tree_, metadata, a_i, s_i, t0, theta_step, nstates, nevfs, ngenes , ncells)

		simulated_res$counts[[sprintf('a%.2f_s%.2f', a_i, s_i)]] <- OU_sim[['counts']]
		simulated_res$evfs[[sprintf('a%.2f_s%.2f', a_i, s_i)]] <- OU_sim[['evfs']]

		write.csv(OU_sim$counts, paste0(outdir, '/', out_prefix, '_alpha_', a_i, '_sigma_', s_i, '_OUEVF_counts.csv'))
        write.csv(OU_sim$evfs, paste0(outdir, '/', out_prefix,'_alpha_', a_i, '_sigma_', s_i, '_OUEVF_evf_counts.csv'))
	 }

	 log_message('Data generation complete', verbose =TRUE)

	return(list('scLT' = test_dataset, 'simulated_OU' = simulated_res)) 

}




##OUwie Simulator##

#written by Jeremy M. Beaulieu

#Simulates the Hansen model of continuous characters evolving under discrete selective
#regimes. The input is a tree of class "phylo" that has the regimes as internal node labels
#and a data file the contains the regime states for each species. The trait file must be in
#the following order: Species names then Regime. The user must specify the parameters values
#for each simulation (i.e. alpha, sigma.sq, theta0, theta).

##The following examples assume 2 selective regimes and different models can be specified:
#single rate Brownian motion BM1: alpha=c(0,0); sigma.sq=c(0.9); theta0=0; theta=c(0,0)
#two rate Brownian motion BMS: alpha=c(0,0); sigma.sq=c(0.45,.9); theta0=0; theta=c(0,0)
#global OU (OU1): alpha=c(0.1,0.1); sigma.sq=c(0.9,0.9); theta0=1; theta=c(1,1)
#normal OU (OUSM): alpha=c(0.1,0.1); sigma.sq=c(0.9,0.9); theta0=0; theta=c(1,2)
#multiple sigmas (OUSMV): alpha=c(0.1,0.1); sigma.sq=c(0.45,0.9);
#multiple alphas (OUSMA): alpha=c(0.5,0.1); sigma.sq=c(0.9,0.9); theta0=0; theta=c(1,2)
#multiple alphas and sigmas (OUSMVA): alpha=c(0.5,0.1); sigma.sq=c(0.45,0.9); theta0=0; theta=c(1,2)

OUwie.sim.edited <- function(phy=NULL, data=NULL, simmap.tree=FALSE, root.age=NULL, scaleHeight=FALSE, 
    alpha=NULL, sigma.sq=NULL, theta0=NULL, theta=NULL, mserr="none", shift.point=0.5, fitted.object=NULL, get.all=FALSE){
    mserr_vector <- NA
    if(!is.null(fitted.object)) {
        if(grepl("BM", fitted.object$model) | grepl("OU1", fitted.object$model)) {
            stop(paste("not implemented yet for ", fitted.object$model))
        }
        if(!is.null(alpha) | !is.null(theta0) | !is.null(theta)) {
            stop("You're passing in parameters to simulate from AND a fitted object to simulate under. You can do one or the other")
        }
        phy <- fitted.object$phy
        data <- cbind(phy$tip.label, fitted.object$data)
        alpha <- fitted.object$solution['alpha',]
        alpha[which(is.na(alpha))] <- 0
        sigma.sq <- fitted.object$solution['sigma.sq',]
        
        if(mserr != "none"){
            warning("measurement error is not yet handled for simulations from fitted.object")
        }
        
        if (fitted.object$root.station == TRUE | fitted.object$root.station==FALSE){
            if (fitted.object$model == "OU1"){
                theta <- matrix(t(fitted.object$theta[1,]), 2, length(levels(fitted.object$tot.states)))[1,]
                theta0 <- theta[phy$node.label[1]]
            }
        }
        if (fitted.object$root.station == TRUE | !grepl("OU", fitted.object$model)){ # BM1 or BMS as well
            if (fitted.object$model != "OU1"){
                theta <- matrix(t(fitted.object$theta), 2, length(levels(fitted.object$tot.states)))[1,]
                theta0 <- theta[phy$node.label[1]]
            }
        }
        if (fitted.object$root.station == FALSE & grepl("OU", fitted.object$model)){
            if (fitted.object$model != "OU1"){
                if(fitted.object$get.root.theta == TRUE){
                    theta.all <- matrix(t(fitted.object$theta), 2, 1:length(levels(fitted.object$tot.states))+1)[1,]
                    theta <- theta.all[2:length(theta.all)]
                    theta0 <- theta.all[1]
                }else{
                    int.states <- factor(phy$node.label)
                    phy.tmp <- phy
                    phy.tmp$node.label <- as.numeric(int.states)
                    theta <- matrix(t(fitted.object$theta), 2, length(levels(fitted.object$tot.states)))[1,]
                    theta0 <- theta[phy.tmp$node.label[1]]
                }
            }
        }
    }

    bt <- branching.times(phy)
    if(is.null(root.age)){
        if(any(bt<0)){
            cat("Warning: Looks like your tree is producing negative branching times. Adjusting to make branches positive." )
            bt <- bt - min(bt)
        }
    }

    #Makes sure the data is in the same order as the tip labels
    if(simmap.tree == FALSE){
        #This is annoying, but the second column has to be in there twice otherwise, error.
        if(mserr == "none"){
            data <- data.frame(data[,2], data[,2], row.names=data[,1])
        }
        if(mserr == "known"){
            data <- data.frame(data[,2], data[,3], row.names=data[,1])
            mserr_vector <- data[,2] #because we've shifted things over
        }
        if(is.numeric(mserr)){
            if(length(mserr) == length(phy$tip.label)){
                data <- data.frame(data[,2], mserr, row.names=data[,1])
                mserr_vector <- mserr
            }   
            if(length(mserr)==1){
                data <- data.frame(data[,2], rep(mserr, length(phy$tip.label)), row.names=data[,1])
                mserr_vector <- rep(mserr, length(phy$tip.label))
            }
            mserr <- "known"
        }

        data <- data[phy$tip.label,]

        n <- max(phy$edge[,1])
        ntips <- length(phy$tip.label)

        int.states <- factor(phy$node.label)
        phy$node.label <- as.numeric(int.states)
        tip.states <- factor(data[,1])
        data[,1] <- as.numeric(tip.states)
        tot.states <- factor(c(phy$node.label,as.character(data[,1])))
        k <- length(levels(tot.states))

        regime <- matrix(rep(0,(n-1)*k), n-1, k)

        #Obtain root state and internal node labels
        root.state <- phy$node.label[1]
        int.state <- phy$node.label[-1]

        #New tree matrix to be used for subsetting regimes
        edges <- cbind(c(1:(n-1)),phy$edge,MakeAgeTable(phy, root.age=root.age))
        if(scaleHeight == TRUE){
            edges[,4:5] <- edges[,4:5]/max(MakeAgeTable(phy, root.age=root.age))
            root.age <- 1
        }
        
        edges <- edges[sort.list(edges[,3]),]
        mm <- c(data[,1],int.state)

        regime <- matrix(0,nrow=length(mm),ncol=length(unique(mm)))
        #Generates an indicator matrix from the regime vector
        for (i in 1:length(mm)) {
            regime[i,mm[i]] <- 1
        }
        #Finishes the edges matrix
        edges <- cbind(edges,regime)

        #Resort the edge matrix so that it looks like the original matrix order
        edges <- edges[sort.list(edges[,1]),]

        oldregime <- root.state

        alpha <- alpha
        alpha[alpha==0] <- 1e-10
        sigma <- sqrt(sigma.sq)
        theta  <- theta

        x <- matrix(0, n, 1)
        TIPS <- 1:ntips
        ROOT <- ntips + 1L
        x[ROOT,] <- theta0

        for(i in 1:length(edges[,1])){
            anc <- edges[i,2]
            desc <- edges[i,3]
            oldtime <- edges[i,4]
            newtime <- edges[i,5]
            if(anc%in%edges[,3]){
                start <- which(edges[,3]==anc)
                oldregime <- which(edges[start,6:(k+5)]==1)
            }else{
                #For the root:
                oldregime <- root.state
            }
            newregime=which(edges[i,6:(k+5)]==1)

            if(oldregime==newregime){
                x[edges[i,3],] <- (x[edges[i,2],]*exp(-alpha[oldregime]*(newtime-oldtime))) + (theta[oldregime]*(1-exp(-alpha[oldregime]*(newtime-oldtime)))) + (sigma[oldregime]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[oldregime]*(newtime-oldtime)))/(2*alpha[oldregime])))
            }else{
                shifttime <- newtime-((newtime-oldtime) * shift.point)
                epoch1 <- (x[edges[i,2],]*exp(-alpha[oldregime]*(shifttime-oldtime))) + (theta[oldregime]*(1-exp(-alpha[oldregime]*(shifttime-oldtime)))) + (sigma[oldregime]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[oldregime]*(shifttime-oldtime)))/(2*alpha[oldregime])))
                oldtime <- shifttime
                newtime <- newtime
                x[edges[i,3],] <- (epoch1*exp(-alpha[newregime]*(newtime-oldtime))) + (theta[newregime]*(1-exp(-alpha[newregime]*(newtime-oldtime)))) + (sigma[newregime]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[newregime]*(newtime-oldtime)))/(2*alpha[newregime])))
            }
        }

        if(get.all == TRUE){
            sim.dat <- matrix(,length(x),3)
            sim.dat <- data.frame(sim.dat)

            sim.dat[,1] <- NA
            sim.dat[TIPS,1] <- phy$tip.label
            sim.dat[,2] <- c(data[,1], phy$node.label)
            sim.dat[,3] <- x
            
            if(mserr == "known"){
                for(i in TIPS){
                    sim.dat[i,3] <- rnorm(1,sim.dat[i,3],data[i,2])
                }
            }
        }else{
            sim.dat <- matrix(,ntips,3)
            sim.dat <- data.frame(sim.dat)

            sim.dat[,1] <- phy$tip.label
            sim.dat[,2] <- data[,1]
            sim.dat[,3] <- x[TIPS]
            
            if(mserr == "known"){
                for(i in TIPS){
                    sim.dat[i,3] <- rnorm(1,sim.dat[i,3],data[i,2])
                }
            }
        }

        colnames(sim.dat)<-c("Genus_species","Reg","X")
        if(mserr=="known"){
            sim.dat$mserr <- mserr_vector   
        }
    }
    sim.dat
}
