## For one iter 

####### SUPPORT FUNCTIONS #######
# This assumes a stationary theta. -- double check that this makes sense
# Might need to make appropriate changes for BM1. 

transformPhy.new <- function(phy, map, Alpha_j, Sigma_j, tip.paths = NULL) {
    # This is the new transformPhy from OUwie but reduced so that it only needs a single alpha and sigma command. 
    nTip <- length(phy$tip.label)

    #ensure we traverse edges from root toward tips so Ato[ancestor] is ready
    phy <- ape::reorder.phylo(phy, "cladewise")

    nEdge <- nrow(phy$edge)
    nNodeTot <- nTip + phy$Nnode

    Ato <- numeric(nNodeTot)
    root <- nTip + 1
    Ato[root] <- 0
    D <- V_Tilde <- numeric(nEdge)
    #cat('Alpha_j:', Alpha_j, 'Sigma_j:', Sigma_j, '\t')
    for (i in 1:nEdge) {
        anc <- phy$edge[i, 1]
        des <- phy$edge[i, 2]
        Map_i <- map[[i]]
        Acur <- Ato[anc] 
        w <- 0
        v <- 0
        for (j in 1:length(Map_i)) {
            dt <- as.numeric(Map_i[j])
            tmp.w <- Alpha_j * dt
            w <- w + tmp.w
            if (Alpha_j == 0) {
                tmp.v <- Sigma_j * exp(2 * Acur) * dt
            }else{
                Aend <- Acur + Alpha_j * dt
                tmp.v <- Sigma_j * (exp(2 * Aend) - exp(2 * Acur)) / (2 * Alpha_j)
            }
            v <- v + tmp.v
            Acur <- Acur + Alpha_j * dt
        }
        V_Tilde[i] <- v
        D[i] <- w
        Ato[des] <- Acur
    }
    phy$edge.length <- V_Tilde
    #calculates the diagonal matrix for each tip i
    DiagWt <- numeric(nTip)
    names(DiagWt) <- phy$tip.label

    if (is.null(tip.paths)) {
        for (i in 1:nTip) {
            DiagWt[i] <- exp(-sum(D[getPathToRoot(phy, i)]))
        }
    }else{
        for (i in 1:nTip) {
            DiagWt[i] <- exp(-sum(D[tip.paths[[i]]]))
        }
    }
    obj <- list(tree = phy, diagWeight = DiagWt)
    return(obj)
}


### WEIGHT MAT FUNCTION
#Weight matrix generator taken from Butler and King (2004) and modified to allow multiple alpha parameters
#written by Jeremy M. Beaulieu
#Taken from Beaulieu 2012 and unmodified to only need one alpha parameter. 
weight.mat <- function(phy, edges, alpha, root.state, root.age = NULL, scaleHeight = FALSE, assume.station = TRUE, shift.point = 0.5){
    # THis is the updated version of weight.mat from the OUwie package. Might change this out for the code from Claude but for now just want to see if this makes sense. 
    age.table <- MakeAgeTable(phy, root.age = root.age)
    Tmax <- max(age.table)
    n <- max(phy$edge[, 1])
    ntips <- length(phy$tip.label)
    
    if(is.null(root.state)) {
        root.state <- which(edges[dim(edges)[1], ] == 1) - 5
        edges <- edges[-1 * dim(edges)[1], ]
    }

    mm <- dim(edges)
    k <- length(6:mm[2])
    pp <- prop.part(phy)

    root_node <- setdiff(unique(edges[, 2]), unique(edges[, 3]))[1]
    max_node_id <- max(c(edges[, 2], edges[, 3]))
    Ato <- numeric(max_node_id)
    Ato[root_node] <- 0

    nodevar.root.tot <- rep(0, max(edges[, 3]))
    nodevar.k <- rep(0, max(edges[, 3]))

    W <- matrix(0, ntips, k)

    for(j in 1:k){
        n.cov.root.tot <- matrix(0, n, 1)
        n.cov.k <- matrix(0, n, 1)
        Ato[] <- 0
        Ato[root_node] <- 0

        for(i in 1:nrow(edges)){
            anc <- edges[i, 2]
            desc <- edges[i, 3]

            #cumulative A at start of the edge
            Acur <- Ato[anc]

            #per edge increments
            nodevar.root.tot[i] <- 0
            nodevar.k[i] <- 0

            oldtime <- edges[i, 4]
            newtime <- edges[i, 5]
            if(anc %in% edges[, 3]){
                start <- which(edges[, 3] == anc)
                oldregime <- which(edges[start, 6:(k + 5)] == 1)
            }else{
                oldregime <- root.state
            }
            newregime <- which(edges[i, 6:(k + 5)] == 1)
            if(oldregime == newregime){
                dt <- newtime - oldtime
                nodevar.root.tot[i] <- nodevar.root.tot[i] - alpha * dt
                if(newregime == j){
                    nodevar.k[i] <- nodevar.k[i] + (exp(Acur + alpha * dt) - exp(Acur))
                }
                Acur <- Acur + alpha * dt

            }else{
                shifttime <- newtime - ((newtime - oldtime) * shift.point)
                # epoch 1 - oldregime
                dt1 <- shifttime - oldtime
                nodevar.root.tot[i] <- nodevar.root.tot[i] - alpha * dt1
                if(oldregime == j){
                    nodevar.k[i] <- nodevar.k[i] + (exp(Acur + alpha * dt1) - exp(Acur))
                }
                Acur <- Acur + alpha * dt1
                #epoch 2 - newregime
                dt2 <- newtime - shifttime
                nodevar.root.tot[i] <- nodevar.root.tot[i] - alpha * dt2
                if(newregime == j){
                    nodevar.k[i] <- nodevar.k[i] + (exp(Acur + alpha * dt2) - exp(Acur))
                }
                Acur <- Acur + alpha * dt2
            }
            n.cov.k[desc, ] <- nodevar.k[i]
            n.cov.root.tot[desc, ] <- nodevar.root.tot[i]
    
            Ato[desc] <- Acur
        }
        w.k <- mat.gen(phy, n.cov.k, pp)
        w.root.tot <- mat.gen(phy, n.cov.root.tot, pp)
        W[, j] <- exp(diag(w.root.tot)) * diag(w.k)
    }
    w.root.tot <- mat.gen(phy, n.cov.root.tot, pp)
    w_root <- exp(diag(w.root.tot))
    
    if (assume.station == TRUE) {
        W[, root.state] <- W[, root.state] + w_root
    }else{
        W <- cbind(w_root, W)
    }
        #Restandardizes W so that the rows sum to 1 -- Generalized. Will reduce to the simpler model if assuming 1 alpha parameter, but when alpha varies by regime they will sum to 1 (though proportionally should be ok).
    W <- W / rowSums(W)
    return(W)
}

HoAne_transform_tree <- function(phy, model, add.root, alpha, sigma){
  if (model == 'BM1') {
          mm <- 'BM'
        } else if (add.root == TRUE) {
          mm <- 'OUfixedRoot'
        } else {
          mm <- 'OUrandomRoot'
        }

  transformed.tree <- phylolm::transf.branch.lengths(phy, model = mm, parameters = list(alpha = alpha, sigma = sigma)) # alpha should be ignored 
  return(transformed.tree)
}

##Matrix generating function taken from vcv.phylo in ape:
mat.gen<-function(phy,piece.wise,pp){
    phy <- reorder(phy, "pruningwise")
    n <- length(phy$tip.label)
    anc <- phy$edge[,1]
    des <- phy$edge[,2]
    ep <- piece.wise[,1]
    comp <- numeric(n + phy$Nnode)
    mat <- matrix(0, n, n)
    
    for (i in length(anc):1) {
        focal <- comp[anc[i]]
        comp[des[i]] <- focal + ep[des[i]]
        j <- i - 1L
        while (anc[j] == anc[i] && j > 0) {
            left <- if (des[j] > n) pp[[des[j] - n]] else des[j]
            right <- if (des[i] > n) pp[[des[i] - n]] else des[i]
            mat[left, right] <- mat[right, left] <- focal
            j <- j - 1L
        }
    }
    diag.elts <- 1 + 0:(n - 1)*(n + 1)
    mat[diag.elts] <- comp[1:n]
    
    mat
}


##### OU VAR-COVAR FUNCTIONS 
#OU variance-covariance matrix generator - modified to include alpha and sigma independently. 

#written by Jeremy M. Beaulieu
varcov.ou <- function(phy, edges, alpha, sigma, root.state, scaleHeight=FALSE, assume.station=TRUE){
    root.age=NULL
    shift.point=.5
    if(assume.station == TRUE){
        alpha=alpha
        sigma=sigma 
        vcv <- quickVCV(phy=phy, alpha=alpha, sigma.sq=sigma, scaleHeight=scaleHeight)
    }else{
        if(is.null(root.state)) {
            root.state <- which(edges[dim(edges)[1],]==1)-5
            edges <- edges[-1*dim(edges)[1],]
        }
        n=max(phy$edge[,1])
        ntips=length(phy$tip.label)
        
        mm<-dim(edges)
        k<-length(6:mm[2])
        
        pp <- prop.part(phy)
        oldregime=root.state
        nodevar1=rep(0,max(edges[,3]))
        nodevar2=rep(0,max(edges[,3]))
        n.cov1=matrix(rep(0,n), n, 1)
        n.cov2=matrix(rep(0,n), n, 1)
        
        for(i in 1:length(edges[,1])){
            anc <- edges[i,2]
            oldtime <- edges[i,4]
            newtime <- edges[i,5]
            if(anc%in%edges[,3]){
                start <- which(edges[,3]==anc)
                oldregime <- which(edges[start,6:(k+5)]==1)
            }
            else{
                    #For the root:
                oldregime=root.state
            }
            newregime=which(edges[i,6:(k+5)]==1)

            # Simplify since only one alpha and sigma value. 
            nodevar1[i] <- alpha*(newtime-oldtime)
            nodevar2[i] <- sigma*((exp(2*alpha*newtime)-exp(2*alpha*oldtime))/(2*alpha))
           
            oldregime <- newregime
            n.cov1[edges[i,3],] <- nodevar1[i]
            n.cov2[edges[i,3],] <- nodevar2[i]
        }

        vcv1 <- mat.gen(phy,n.cov1,pp)
        vcv2 <- mat.gen(phy,n.cov2,pp)
        if(any(abs(diff(alpha)) > 0)){
            species.variances <- diag(vcv1)
            species.total.variances <- matrix(0, dim(vcv1)[2], dim(vcv1)[2])
            count=0
            for(i in 1:dim(vcv1)[2]) {
                for(j in 1:dim(vcv1)[2]){
                    species.total.variances[i,j] <- exp(-(species.variances[i] + species.variances[j]))
                    count=count+1 # the count is always watching
                }
            }
            vcv <- species.total.variances * vcv2
        }else{
            if(is.null(root.age)){
                root.age <- max(branching.times(phy))
            }
            vcv <- exp(-2*alpha*max(root.age)) * vcv2
        }
    }
    return(vcv)
}


## Quick VCV maker of OU1 and OUM -- since alpha and sigma.sq are constants, and since the regime does not matter, we can just do a simple plug and chug. It's a tad slower, but mostly for testing purposes.
quickVCV <- function(phy, alpha, sigma.sq, scaleHeight){
    phy$node.label <- NULL
    vcv <- matrix(0, Ntip(phy), Ntip(phy))
    split.times <- branching.times(phy)
    if(scaleHeight == TRUE){
        split.times <- split.times/max(split.times)
    }
    tot.time <- max(split.times)
    for(i in 1:Ntip(phy)){
        for(j in 1:Ntip(phy)){
            if(i == j){
                dij <- 0
                tij <- tot.time
                vcv[i,j] <- (sigma.sq /(2 * alpha)) * exp(-2 * alpha * dij)
            }else{
                split <- getMRCA(phy, tip=c(phy$tip.label[i], phy$tip.label[j]))
                dij <- split.times[which(names(split.times)==split)]
                tij <- tot.time - dij
                vcv[i,j] <- (sigma.sq /(2 * alpha)) * exp(-2 * alpha * dij)
            }
        }
    }
    return(vcv)
}

######## NEW UTILITY FUNCTIONS -- Claude support 
compute_W_matrix <- function(tree_info, alpha, normalize=FALSE, add.root=TRUE) {
  n_leaves <- tree_info$n_leaves
  regimes <- tree_info$unique_regimes
  n_regimes <- length(regimes)
  root.state <- tree_info$root.state 
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

  root <- sapply(1:n_leaves, function(i){exp(-alpha * tree_info$leaf_dists[i])})
  W[, root.state] <- W[,root.state]+root # changing so that the root is already a column in here! 
  # just kidding thats a larger change because need to change the root.state variable and I am hesitant to do that. 
  # also adding the cbind ensures that we get two columns for BM1. 
  #if (!add.root){
  #      W[, root.state] <- W[,root.state]+root
  #  } else {
  #      W <- cbind(root, W)
   # }

  if (normalize){
    W <- W/rowSums(W)
  }
  
  return(W)
}

compute_VCV<- function(tree_info, alpha, sigma, add.root=TRUE) {
  n <- tree_info$n_leaves
  t <- tree_info$leaf_dists
  s <- tree_info$shared_lengths
  
  V <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in i:n) {
      sij <- s[i, j] # shared time 
      dij <- t[i] + t[j] - 2*sij  # distance for tips i, j to root minus 2*shared time (aka root --> MRCA) 
      if (i == j & dij != 0 ){stop('Error: distance to self does not equal 0.')}
      # we are assuming that sigma is already squared. 
      # Ho and Ane 2014 -- form that conditions on theta0. 
      if (add.root){
        val <- (sigma / (2 * alpha)) * exp(-alpha * dij) * (1 - exp(-2 * alpha * sij)) # --> fixed root; like ouwie assume.station = FALSE
      } else {
        val <- (sigma / (2 * alpha)) * exp(-alpha * dij) # --> stationary root; like ouwie assume.station = TRUE 
      }
      V[i, j] <- val
      V[j, i] <- val
    }
  }
  # Add small jitter for numerical stability
  #V <- V + diag(1e-6, n) # ridge added in Estep. 
  
  return(V)
}

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
  
  return(paths)
}

#' Get the path from root to a given tip
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
  
  return(path)
}


######## OTHER UTILITY FUNCTIONS --- cite Beaulieu et al / OUwie 
#### THREE POINT SUPPORT FUNCTIONS 
#OU functions for using three-point algorithm of Ho and Ane 2013

#written by James D. Boyko
MakeAgeTable <- function(phy, root.age=NULL){
    if(is.null(root.age)){
        node.ages <- dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))
    }else{
        node.ages <- dateNodes(phy, rootAge=root.age)
    }
    max.age <- max(node.ages)
    table.ages <- matrix(0, dim(phy$edge)[1], 2)
    for(row.index in 1:dim(phy$edge)[1]){
        table.ages[row.index,1] <- max.age - node.ages[phy$edge[row.index,1]]
        table.ages[row.index,2] <- max.age - node.ages[phy$edge[row.index,2]]
    }
    return(table.ages)
}

## vcvbypass.functions.R
getMapFromNode <- function(phy, tipstates, nodestates, shift.point){
  Map <- vector("list", dim(phy$edge)[1])
  Data <- c(tipstates, nodestates)
  NodeStates <- data.frame(Data[phy$edge[,1]], Data[phy$edge[,2]])

  #NodeStates <- cbind(Data[phy$edge[,1]], Data[phy$edge[,2]])
  for(i in 1:dim(phy$edge)[1]){
    from <- as.character(NodeStates[i,1])
    to <- as.character(NodeStates[i,2])
    #from <- NodeStates[i, 1]
    #to <- NodeStates[i, 2]

    if(from == to){
      tmp <- phy$edge.length[i]
      names(tmp) <- from
      Map[[i]] <- tmp
    }else{
      shift.time <- shift.point * phy$edge.length[i]
      tmp <- c(phy$edge.length[i] - shift.time, shift.time)
      names(tmp) <- c(from, to)
      Map[[i]] <- tmp
    }
  }
  return(Map)
}

getPathToRoot <- function(phy, tip){
  nTip <- length(phy$tip.label)
  root <- nTip + 1
  path <- 0
  count <- 1
  while(tip != root){
    tip.ind <- which(phy$edge[,2] == tip)
    path <- c(path, tip.ind)
    count <- count + 1
    tip <- phy$edge[tip.ind,1]
  }
  path <- path[-1]
  
  return(path)
}

# Need to be careful because this needs to be of length 1 if BM1 invert. 
# k = 1 for OU1 and BM1. n = nstates for OUM. This SHOULD work out to be the adjustment needed for the three.point algorithm.  
establishBounds <- function(model, algorithm, k, get.root.theta, em=NULL, lb=1e-9, ub=1e5){
    if (model == 'BM1'){
        baseline <- 1
    } else {
        baseline <- 2
    }
    if (is.null(em)){
        adjust <- 0
    } else {
        if (em %in% c('NORM', 'NB')){
            adjust <- 2
        } else if (em == 'POIS') {
            adjust <- 1
        } else {
            adjust <- 0
        }
    }

    np <- baseline + adjust
    if (algorithm == 'three.point'){
        if (model == 'BM1' & get.root.theta == FALSE){
            np <- np
        } else {
         np <- np + k
        }
    }

    lower <- rep(log(lb), np )
    upper <- rep(log(ub), np )

    return(list(lower, upper))
}

establishIP <- function(model, data, algorithm, phy, k, Tmax, em, starting.vals, est.mean=TRUE){
    x <- data[,2]
    if (ncol(data) > 2){
        xobs <- data[,3]
    } else {
        xobs <- data[,2]
    }

    # POST scale-height Tmax 
    if(is.null(starting.vals)){
        if (model == 'BM1'){
            ip <-  c(mean(pic(x, phy)^2))
        } else {
            ip <- c(log(2)/Tmax, mean(pic(x, phy)^2))
        }
    } else {
        ip <- starting.vals
    }

    if(algorithm == "three.point"){
        means.by.regime <- with(data, tapply(data[,2], data[,1], mean))
        names(means.by.regime) <- NULL
        ip <- c(ip, means.by.regime)
    }

    if (!is.null(em)){
        if (em == "NORM"){
            emip <- c( sd(xobs))
        } else if (em == 'NB'){
            emip <- c(1)
        } else if (em == 'POIS'){
            emip <- c()
        } else {
        stop('Specified model not supported. Try from list: NORM, NB, POIS')
        }

        if (est.mean == TRUE){
            mu <- mean(xobs)
            emip <- c(mu, emip)
        }

        ip <- c(ip, emip)
        
    }

    return(ip)
}


######### preprocessing functions 
makeEdges <- function(phy, data, model, k, Tmax_init, int.state, root.age, scaleHeight){
    #New tree matrix to be used for subsetting regimes
    n <- max(phy$edge[,1])
    edges=cbind(c(1:(n-1)), phy$edge, MakeAgeTable(phy, root.age=root.age))
    if(scaleHeight==TRUE){
        edges[,4:5]<-edges[,4:5]/Tmax_init
    }

    edges=edges[sort.list(edges[,3]),]

    if (model %in% c('BM1', 'OU1')) {
        regime <- matrix(0,nrow=length(edges[,1]),ncol=k)
        regime[,1] <- 1
        if (k>1) {
            regime[,2:k] <- 0
        }
    } else {
        mm <- c(data[,1],int.state)
        regime <- matrix(0,nrow=length(mm),ncol=k)
        #Generates an indicator matrix from the regime vector
        for (i in 1:length(mm)) {
            regime[i,mm[i]] <- 1
        }
    }
    edges=cbind(edges,regime)

    #Resort the edge matrix so that it looks like the original matrix order
    edges <- edges[sort.list(edges[,1]),]

    return (edges)
}

# here we want a matrix of n-ou-param by n-states with the corresponding init parameters index in the correct spot. 
establishIndex <- function(model, k, algorithm, get.root.theta){
    index.mat <- matrix(0, 2, k)
    if (model == "BM1"){
        init.index <- c(NA, 1)
    } else {
        init.index <- c(1, 2)
    }

    index.mat[] <- init.index
    if(algorithm == "three.point"){
        index.mat <- rbind(index.mat, rep(NA,k))
        max.par.so.far <- max(index.mat, na.rm = TRUE)
        if (model == 'OU1') {
            index.mat[is.na(index.mat)] <- 3
        } else if (model == 'OUM') {
            index.mat[is.na(index.mat)] <- (max.par.so.far + 1):(max.par.so.far + k)
        }
    }
    return(index.mat)
}


####=======#### 
check_bounds <- function(params, lower_bounds, upper_bounds, verbose=FALSE){
    below_lower <- params < lower_bounds
    if (any(below_lower)) {
        if (verbose) {
            cat("Parameters below lower bound:\n")
            print(names(params)[below_lower])
            cat("  Original: ", params[below_lower], "\n")
            cat("  Clamped to: ", lower_bounds[below_lower], "\n\n")
        }
        params[below_lower] <- lower_bounds[below_lower]
    }

    above_upper <- params > upper_bounds
    if (any(above_upper)) {
        if (verbose) {
            cat("Parameters above upper bound:\n")
            print(names(params)[above_upper])
            cat("  Original: ", params[above_upper], "\n")
            cat("  Clamped to: ", upper_bounds[above_upper], "\n\n")
        }
        params[above_upper] <- upper_bounds[above_upper]
    }

    if (!any(below_lower | above_upper) & verbose) {
        cat("All parameters within bounds.\n")
    }

    return(params)
}


lineage_smooth <- function(tree_, lognorm__, ka, s=NULL) {
    adj_mat <- as.matrix(cophenetic(tree_), sparse = TRUE)
    
    tip_order <- row.names(adj_mat)
    ln_ordered <- as.matrix(lognorm__[tip_order, ])
    
    m <- nrow(adj_mat)
    sigmas <- rep(0, m)
    windows <- matrix(nrow = nrow(adj_mat), ncol = ncol(adj_mat))
    
    for (i in 1:m){
        dists <- adj_mat[i, ]
        max_dist <- sort(dists)[ka+1]
        windows[i,] <- dists > max_dist # store the boolean of whether the distance is greater or less than the max dist. TRUE = make 0. 
        sigmas[i] <- max_dist
    }
    
    # Overwrite sigmas if there is an established s, otherwise do the dynamic option. 
    if (!is.null(s)){
        sigmas <- rep(s, m)}
    
    sigmas[sigmas == 0] <- 1 # replace zeros with 1 to remove division by zero errors, not sure if this is the way to go though. 

    affinity_matrix <- exp(- (adj_mat^2 / ( sigmas) ))
    affinity_matrix <- affinity_matrix + t(affinity_matrix) # make symmetric
    affinity_matrix[windows] <- 0
    
    normalized_affinity <- as.matrix(affinity_matrix / rowSums(affinity_matrix))
    diag(normalized_affinity) <- 1 
    
    smoothed_data <- normalized_affinity %*% ln_ordered
    
    if (sum(rownames(smoothed_data) == tree_$tip.label) != length(tip_order)) {
        stop('Data not in the same order at tree. Check results before moving on.')
    }

    return(smoothed_data)
}

####=======#### 

preprocess.OLD <- function(inputs, reg, gene_name, 
    scaleHeight = FALSE, 
    root.station = FALSE, 
    get.root.theta = FALSE, 
    root.age = NULL ){

    phy <- inputs$inputs[[reg]]$tree
    model <- inputs$inputs[[reg]]$model # model now the actual model being tested NOT the regime. so in OU1, BM1, or OUM. 
    reg_col <- inputs$inputs[[reg]]$rcol # this is the name of the column in the object with the regime. 
    
    if (!gene_name %in% colnames(inputs$meta_data)){
        stop('Gene name not in data. Try again.')
    }

    if (!reg_col %in% colnames(inputs$meta_data)){
        stop('Model not annotated in metadata. Add to data.')
    }
        
    data.orig <- inputs$meta_data[, c('species', reg_col, gene_name)]
    phy <- reorder.phylo(phy, "cladewise")
        
    if(model == "BM1" & root.station == FALSE) { get.root.theta = TRUE } ## Need to check whats up with this. 
        
    data <- data.frame(data.orig[, c(reg_col, gene_name)], row.names = data.orig[, 'species'])
    #print(head(data))
        
    # Initialize constants. 
    shift.point=0.5
    data <- data[phy$tip.label,]
    tip.states.cp <- factor(data[,1]) # not sure if this is used ever. 
    ntips <- length(phy$tip.label)
    Tmax.i <- max(MakeAgeTable(phy, root.age=root.age)) # with root.age == NULL 
        
        #### Establish a bunch of values used later. 
    if (model %in% c('BM1', 'OU1')){
          #Begins the construction of the edges matrix -- similar to the ouch format##
          #Makes a vector of absolute times in proportion of the total length of the tree
          k <- 1
          phy$node.label <- rep(1, Nnode(phy))
          tip.states <- factor(rep(1, length(data[,1])))
          int.states <- factor(phy$node.label)
          tot.states <- factor(c(phy$node.label,tip.states))
          
          #Obtain root state -- for both models assume the root state to be 1 since no other state is used even if provided in the tree
          root.state <- 1
    } else {
          #Obtain a a list of all the regime states. This is a solution for instances when tip states and
          #the internal nodes are not of equal length:
          tot.states <- factor(c(phy$node.label,as.character(data[,1])))
          k <- length(levels(tot.states))
          int.states <- factor(phy$node.label)
          phy$node.label <- int.states# as.numeric(int.states)
          tip.states <- factor(data[,1])
          data[,1] <- as.numeric(tip.states)
          #Obtain root state and internal node labels
          root.state <- phy$node.label[1]
          int.state <- phy$node.label[-1]
    }
        
    edges <- makeEdges(phy, data, model, k, Tmax.i, int.state, root.age, scaleHeight)
        
    map <- getMapFromNode(phy, tip.states, int.states, shift.point)
    if(scaleHeight==TRUE){
        map <- lapply(map, function(x) x/Tmax.i)
    }
        
    if(scaleHeight==TRUE){
        phy$edge.length <- phy$edge.length/Tmax.i
        Tmax <- 1
        root.age <- 1
    } else {
        Tmax <- Tmax.i
    }
        
    phy$states <- factor(data[phy$tip.label, 1])

    X <- inputs$meta_data[, gene_name] #res$X
    names(X) <- inputs$meta_data[, 'species'] #row.names(res) #inputs$meta_data[, 'species']

    map <- getMapFromNode(phy, tip.states, int.states, shift.point)
    if(scaleHeight==TRUE){
        phy$edge.length <- phy$edge.length/Tmax.i
        map <- lapply(map, function(x) x/Tmax.i)
        Tmax <- 1
        root.age <- 1
    } else {
        Tmax <- Tmax.i
    }

    tip.paths <- lapply(1:length(X), function(x) getPathToRoot(phy, x))

    alpha.i <- log(2)/Tmax
    tau.i <- sd(as.vector(X))
    sigma.i <- mean(pic(X, phy)^2)
    means.by.reg.init <-  tapply(X, data[,1], mean)

    if (get.root.theta == TRUE){
        # For BM1: Assumed to be the same as the expected value, so pulled once in the optimization and for the Wmat calculation 
        # For OU+: NOT assumed to be the same so included in the parameter opimization twice. 
        # This will solve the issues with calculating the W matrix but requires some reformatting for the optimzation problem. 
        assume.station = FALSE 
        # build theta0 into theta if asking to get root theta. 
        if (model == 'BM1'){
            theta <- c( means.by.reg.init) 
            names(theta) <- c('root')
        } else {
            theta <- c(means.by.reg.init[root.state], means.by.reg.init) 
            names(theta) <- c('root', names(means.by.reg.init))
        }
    } else {
        assume.station = TRUE
        theta <-c(means.by.reg.init)
    }

    defaults <- list('root.state' = root.state,  'assume.station' = assume.station, 
                    'scaleHeight' = scaleHeight , 
                    'model' = model, 
                    tp_const = list( 'tippath' = tip.paths, 'map' = map))

    cat('Running with defaults:\n')
    cat(sprintf('\troot.state = %s | assume.station = %s | scaleHeight = %s | model = %s\n', root.state, assume.station, scaleHeight, model))

    #names(means.by.reg.init) <- NULL
    if (model == 'BM1'){
        # Need to add a constant tiny alpha for BM1 to all the above code. 
        par0 <- list('tau' = tau.i, alpha = 1e-10, 'sigma' = sigma.i, 'theta' = theta) # tau, alpha, sigma2, means.by.regime. 
    } else {
        par0 <- list('tau' = tau.i, 'alpha' = alpha.i, 'sigma' = sigma.i, 'theta' = theta) # tau, alpha, sigma2, means.by.regime. 
    }
   
    return(list(param_init = par0, defaults = defaults, X = X, phy = phy, edges = edges))

}

##################################################################################

preprocessTree <- function(inputs, reg, 
    gene_cols = NULL, 
    scaleHeight = FALSE, 
    root.fixed = FALSE, # root.fixed == FALSE --> assume.station --> true ---> add.root = FALSE 
    root.age = NULL, 
    weight_type = 'alt', # whether to use the new or old way to calculate the varcovar and weights matrices 
    algorithm = 'three.point',  # which algorithm to use (inverse or threepoint )
    ttransf = 'ouwie', # ha2014 toggles the other option. 
    skipTau = FALSE
    ){

    phy <- inputs$inputs[[reg]]$tree
    model <- inputs$inputs[[reg]]$model # model now the actual model being tested NOT the regime. so in OU1, BM1, or OUM. 
    reg_col <- inputs$inputs[[reg]]$rcol # this is the name of the column in the object with the regime. 
    
    if (!reg_col %in% colnames(inputs$meta_data)){
        stop('Model not annotated in metadata. Add to data.')
    }

    if (is.null(gene_cols)){
        gene_cols <- inputs$gene_cols
    } else if (length(intersect(gene_cols, colnames(inputs$meta_data))) == 0){
        stop('No genes found. Check names. Or inputs$gene_cols.')
    }
        
    data.orig <- inputs$meta_data[, c('species', reg_col, gene_cols)]
    phy <- reorder.phylo(phy, "cladewise")
    
    if(model == "BM1" & root.fixed == FALSE) { root.fixed = TRUE } # change the argument if BM1. 
    
    data <- data.frame(data.orig[, c(reg_col, gene_cols)], row.names = data.orig[, 'species'])
        
    # Initialize constants. 
    shift.point=0.5
    data <- data[phy$tip.label,]
    ntips <- length(phy$tip.label)
    Tmax.i <- max(MakeAgeTable(phy, root.age=root.age)) # with root.age == NULL 

        #### Establish a bunch of values used later. 
    if (model %in% c('BM1', 'OU1')){
          #Begins the construction of the edges matrix -- similar to the ouch format##
          #Makes a vector of absolute times in proportion of the total length of the tree
          k <- 1
          phy$node.label <- rep(1, Nnode(phy))
          tip.states <- factor(rep(1, length(data[,1])))
          int.states <- factor(phy$node.label)
          tot.states <- factor(c(phy$node.label,tip.states))
          
          #Obtain root state -- for both models assume the root state to be 1 since no other state is used even if provided in the tree
          root.state <- 1
    } else {
          #Obtain a a list of all the regime states. This is a solution for instances when tip states and
          #the internal nodes are not of equal length:
          tot.states <- factor(c(phy$node.label,as.character(data[,1])))
          k <- length(levels(tot.states))
          int.states <- factor(phy$node.label)
          phy$node.label <- int.states# as.numeric(int.states)
          tip.states <- factor(data[,1])
          data[,1] <- as.numeric(tip.states)
          #Obtain root state and internal node labels
          root.state <- phy$node.label[1]
          int.state <- phy$node.label[-1]
    }
        
    edges <- makeEdges(phy, data, model, k, Tmax.i, int.state, root.age, scaleHeight)
        
    map <- getMapFromNode(phy, tip.states, int.states, shift.point)
    if(scaleHeight==TRUE){
        map <- lapply(map, function(x) x/Tmax.i)
    }
        
    if(scaleHeight==TRUE){
        phy$edge.length <- phy$edge.length/Tmax.i
        Tmax <- 1
        root.age <- 1
    } else {
        Tmax <- Tmax.i
    }
        
    phy$states <- tip.states

    map <- getMapFromNode(phy, tip.states, int.states, shift.point)
    if(scaleHeight==TRUE){
        phy$edge.length <- phy$edge.length/Tmax.i
        map <- lapply(map, function(x) x/Tmax.i)
        Tmax <- 1
        root.age <- 1
    } else {
        Tmax <- Tmax.i
    }

    # legacy parameters for OUwie functions. 
    if (root.fixed == TRUE){
        assume.station = FALSE 
    } else {
        assume.station = TRUE
    }

    #cat('Running with defaults:\n')
    #cat(sprintf('\troot.state = %s | assume.station = %s | scaleHeight = %s | model = %s\n', root.state, assume.station, scaleHeight, model))

    ### GATHER ALTERNATIVE INPUTS FOR THE NEW VERSION OF THE WEIGHTS MATRIX. 
    n_leaves = length(phy$tip.label)
    n_nodes = phy$Nnode + n_leaves
    root_dists <- node.depth.edgelength(phy)
    leaf_indices <- 1:n_leaves
    leaf_dists <- root_dists[leaf_indices]
    node_regimes <- c(tip.states, phy$node.label)
    names(node_regimes) <- 1:n_nodes
    regime_paths <- build_regime_paths(phy, node_regimes, root_dists)
    unique_regimes <- unique(node_regimes)
    alt.root.state <- if (root.fixed) {'root' } else {phy$node.label[1]}
    mrca_matrix <- mrca(phy)  # matrix of MRCA node indices
    shared_lengths <- matrix(0, n_leaves, n_leaves)
    for (i in 1:n_leaves) {
        for (j in i:n_leaves) {
          mrca_node <- mrca_matrix[i, j]
          shared_lengths[i, j] <- root_dists[mrca_node]
          shared_lengths[j, i] <- shared_lengths[i, j]
        }
    }

    alt_weight_mat = list(n_leaves = n_leaves, 
                        regime_paths = regime_paths,
                        unique_regimes = unique_regimes,
                        leaf_dists = leaf_dists,
                        shared_lengths = shared_lengths,
                        root.state = alt.root.state)

    tip.paths <- lapply(1:nrow(data), function(x) getPathToRoot(phy, x))

    defaults_list <- list(root.state = root.state,  
        add.root = root.fixed, 
        assume.station = assume.station, 
        scaleHeight = scaleHeight , 
        model = model, 
        weight_type = weight_type, 
        transf = ttransf, # this controls whether to use the phylolm:transf.branch.lengths or the ouwie version. 
        algorithm = algorithm, 
        skipTau = skipTau, 
        tp_const = list( tippath = tip.paths, map = map), 
        parsed_alt_tree = alt_weight_mat)

    return(list(defaults = defaults_list, tree_const = list(tree = phy, edges = edges, Tmax = Tmax)))
}

preprocessGene <- function(inputs, gene_name, tree_const, root.state, root.fixed, model, skipTau){    
    # Isolating what needs a specific gene to run for parallelization. 
    X <- inputs$meta_data[, gene_name] #res$X
    names(X) <- inputs$meta_data[, 'species'] #row.names(res) #inputs$meta_data[, 'species']
    phy <- tree_const$tree
    Tmax <- tree_const$Tmax

    alpha.i <- log(2)/Tmax
    tau.i <- sd(as.vector(X))
    sigma.i <- tryCatch({mean(pic(X, phy)^2)}, error = function(err){
        stop(err$message)
    })

    means.by.reg.init <-  tapply(X[phy$tip.label], phy$states, mean)

    if (root.fixed){
        # If root.fixed = TRUE it means we are ADDING a root. 
        # For BM1: Assumed to be the same as the expected value, so pulled once in the optimization and for the Wmat calculation 
        # For OU+: NOT assumed to be the same so included in the parameter opimization twice. 
        # This will solve the issues with calculating the W matrix but requires some reformatting for the optimzation problem. 
        
        # build theta0 into theta if asking to get root theta. 

        if (model == 'BM1'){
            theta <- c( means.by.reg.init) 
          #  names(theta) <- c('root')
        } else {
            theta <- c(means.by.reg.init[root.state], means.by.reg.init) 
            names(theta) <- c('root', names(means.by.reg.init))
        }
    } else {
        theta <-c(means.by.reg.init)
    }
    
    #names(means.by.reg.init) <- NULL
    if (model == 'BM1'){
         alpha.i = 1e-10
        # Need to add a constant tiny alpha for BM1 to all the above code. 
    } 

    if (skipTau){
        tau.i <- 0
    }

    par0 <- list('tau' = tau.i, 'alpha' = alpha.i, 'sigma' = sigma.i, 'theta' = theta) # tau, alpha, sigma2, means.by.regime. 
   
    return(list(param_init = par0, X = X))

}


#### OTHER SCOUT UTILS 
log_message <- function(message, log_file = NULL, verbose = FALSE) {
  timestamped_msg <- paste0(Sys.time(), " | ", message)
  
  if (!is.null(log_file)) {
    cat(timestamped_msg, file = log_file, append = TRUE, sep = "\n")
  }
  
  if (verbose) {
    message(timestamped_msg)
  }
}

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
######### FROM runSCOUT.R, prepare_data 
formatSCOUT <- function(tree_path, metadata_path, outpath, 
    species_key = NULL, 
    quant_traits = NULL,
    regimes = "BM1", 
    anc_infer = "ape",
    normalize = FALSE, 
    floor = 0, 
    smoothing_k = NULL, 
    blacklist = NULL, 
    logfile=NULL) {

    create_directory_if_not_exists(outpath)
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
    if (length(intersect(reg_oums, colnames(meta))) != length(reg_oums)){
        stop(sprintf('Specified regimes not found in metadata columns: %s. Please correct. Exiting...', paste(reg_oums, collapse=', ')))
    }

    # check format -- make the states characters just to make this easy. 
    for (col in reg_oums) {
        if (!is.character(meta[[col]])) {
            meta[[col]] <- paste0('r', meta[[col]]) # pasting r in front of any numeric regime. 
        }
    }

    # clean up reg names 
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

    ################### Data Preprocessing ###################

    if (normalize){
        log_message(paste(Sys.time(), 'Log normalizing gene expression counts.'), log_file=logfile,  verbose=TRUE)
        meta[, gene_cols] <- log1p(meta[, gene_cols])
    }

    if (!is.null(smoothing_k)) {
        log_message(paste('Smoothing counts at k =', smoothing_k), log_file=logfile,  verbose=TRUE)
        smooth_meta <- lineage_smooth(phy, meta[, gene_cols], smoothing_k)
        meta[phy$tip.label, gene_cols] <- smooth_meta[phy$tip.label, ]
    } 

    if (floor > 0){
        meta[, gene_cols][meta[, gene_cols] == 0] <- floor
    }

    cln_meta <- meta[phy$tip.label, ]
    cln_meta$OU1 <- 'global'#, length(phy$tip.label))

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
            ### DONT CHANGE THING FROM ANC things fail without any regime shifts for some reason.
            phy1$node.label <- rep('ANC', phy$Nnode)
            rcol <- 'OU1'
            model <- rname 
        }
        log_message(paste0('Regime: ', rname, ' Rcol: ', rcol, ' Model: ', model), verbose = TRUE)
        model_input[[rname]] <- list(model = model, tree = phy1, rcol = rcol)
        phy$node.label <- NULL # reset labels 
    }

    rcols <- unique(unlist(lapply(model_input, function(x) {x$rcol})))
    cln_meta <- cln_meta[, c('species', rcols, gene_cols)]
    rownames(cln_meta) <- NULL

    return(list(inputs= model_input, 
        meta_data = cln_meta, 
        gene_cols = gene_cols, 
        regimes = regimes,
        anc_infer = infer_anc))
}


runSCOUT<- function(idata, 
    lambda1, 
    lambda2, 
    fixed.root, 
    weight_type = 'alt', 
    algorithm = 'three.point', 
    tau_prior_mean = 0.2, # N(sd(X)*0.2, 0.1)
    tau_prior_sd = 0.1, # N(sd(X)*0.2, 0.1)
    alpha_shape = 2, # control on the alpha prior. 
    alpha_rate = 1, 
    M_only = FALSE, 
    skipTau = FALSE, 
    scaleHeight = FALSE, 
    runid=NULL, 
    cores=1, 
    logfile=NULL,
    verbose=FALSE){

    genes <- idata$gene_cols 
    regimes <- idata$regimes 
    combos <- expand.grid(genes, regimes)
    combos <- apply(combos, 2, as.character)
    colnames(combos) <- c('gene_name', 'regime')#, 'burnin', 'lambda1', 'tau_prior_mean', 'tau_prior_sd', 'lambda2')
    code <- if (is.null(runid)) {paste(sample(c(0:9, letters, LETTERS), 8, replace = TRUE), collapse = "") } else {runid}

    cat('Running', nrow(combos), 'combinations. Results will be saved at run ID =', code, '\n')
    if (cores > nrow(combos)){
        cores <- nrow(combos)
    }

    all_tree_info <- list() 
    for (reg in regimes){
        all_tree_info[[reg]] <- preprocessTree(inputs=idata, reg = reg, scaleHeight=scaleHeight, root.fixed=fixed.root, root.age=NULL, skipTau = skipTau)
    }

    total_time <- 0
    plan(multisession, workers = cores)

    start_time <- Sys.time()
            
    log_message(paste0('Start time = ', format(start_time, "%H:%M:%S")[[1]], '\n'), logfile, verbose=TRUE)
    with_progress({
        p <- progressor(along = seq_len(nrow(combos)))
        results <- future_lapply(seq_len(nrow(combos)), function(j){
            # Unpack combos 
            gene_name <- combos[j, 1]
            reg <- combos[j, 2]

            tree_info <- all_tree_info[[reg]]
            defaults <- tree_info$defaults
            tree_constants <- tree_info$tree_const
            gene_info <- preprocessGene(idata, gene_name, tree_constants, defaults$root.state, defaults$add.root, defaults$model, defaults$skipTau) # I can tuck this into run_em later. 
            gene_info$param_init <- verify_initial_params(gene_info$param_init, defaults, FALSE )   

            # Cleaning up naming for the regimes so that the Weights are calculated correctly. 
            if (any(names(gene_info$param_init$theta) != defaults$parsed_alt_tree$unique_regimes)){
                defaults$parsed_alt_tree$unique_regimes <- names(gene_info$param_init$theta)
            }
            if (!'root' %in% defaults$parsed_alt_tree$unique_regimes & defaults$add.root){
                defaults$parsed_alt_tree$unique_regimes <- c('root', defaults$parsed_alt_tree$unique_regimes)
            }

            res1 <- run_em(tree_constants$tree, gene_info$X, gene_info$param_init, tree_constants$edges, defaults, 
                            max_iter = 100, # default
                            tol_rel = 1e-4, # default 
                            tol_abs = 0.01, # smaller than default.  
                            lambda1 = lambda1,
                            lambda2 = lambda2, 
                            diagnose=FALSE , 
                            tau_prior_mean = tau_prior_mean, # both defaults for now. 
                            tau_prior_sd = tau_prior_sd,
                            alpha_shape = alpha_shape, 
                            alpha_rate = alpha_rate, 
                            save=FALSE, # default 
                            save_pre=NULL, 
                            skipE = M_only, # if M-only then skipE == TRUE. 
                            verbose=verbose)

            res1$settings <- as.vector(combos[j, ])
            p()
            return(res1)
            }, future.seed=TRUE)
    })
    end_time <- Sys.time()
    log_message(paste0('End time = ', format(end_time, "%H:%M:%S")[[1]], '\n'), logfile, verbose=TRUE)
    
    return(results)
}


runSCOUT.batches <- function(idata, 
    lambda1, 
    lambda2, 
    fixed.root, 
    weight_type = 'alt', 
    algorithm = 'three.point', 
    tau_prior_mean = 0.2, # N(sd(X)*0.2, 0.1)
    tau_prior_sd = 0.1, # N(sd(X)*0.2, 0.1)
    M_only = FALSE, 
    skipTau = FALSE, 
    runid=NULL, 
    cores=1, 
    batch_size = 50, 
    results_dir = 'results_tmp',
    logfile=NULL,
    verbose=FALSE){

    genes <- idata$gene_cols 
    if (batch_size > length(genes)){
        batch_size = length(genes)
        log_message('WARNING: Batch size bigger than the number of genes. Resetting to number of genes.')
    }

    regimes <- idata$regimes 
    combos <- expand.grid(genes, regimes)
    combos <- apply(combos, 2, as.character)
    colnames(combos) <- c('gene_name', 'regime')#, 'burnin', 'lambda1', 'tau_prior_mean', 'tau_prior_sd', 'lambda2')
    code <- if (is.null(runid)) {paste(sample(c(0:9, letters, LETTERS), 8, replace = TRUE), collapse = "") } else {runid}

    ## SET UP TMP RESULTS 
    if (!dir.exists(results_dir)) dir.create(results_dir)
    results_file <- file.path(results_dir, "combined_results.rds")

    
    all_tree_info <- list() 
    for (reg in regimes){
        all_tree_info[[reg]] <- preprocessTree(inputs=idata, reg = reg, scaleHeight=FALSE, root.fixed=fixed.root, root.age=NULL, skipTau = skipTau)
    }

   
    # Pre-compute batches
    n_combos <- nrow(combos)
    batch_indices <- split(seq_len(n_combos), ceiling(seq_len(n_combos) / batch_size))
    cat('Running', nrow(combos), 'combinations in', length(batch_indices), 'batches. Results will be saved at run ID =', code, '\n')
    if (cores > nrow(combos)){
        cores <- nrow(combos)
    }

    total_time <- 0
    plan(multisession, workers = cores)

    start_time <- Sys.time()
            
    log_message(paste0('Start time = ', format(start_time, "%H:%M:%S")[[1]], '\n'), logfile, verbose=TRUE)


    with_progress({
      p <- progressor(along = batch_indices)
      
      future_lapply(batch_indices, function(batch_idx) {
        
        # Process all combos in this batch
        batch_results <- list()
        
        for (j in batch_idx) {
          # Unpack combos
          gene_name <- combos[j, 1]
          reg <- combos[j, 2]
          tree_info <- all_tree_info[[reg]]
          defaults <- tree_info$defaults
          tree_constants <- tree_info$tree_const
          
          # Preprocess gene
          gene_info <- preprocessGene(idata, gene_name, tree_constants,
                                       defaults$root.state, defaults$add.root,
                                       defaults$model, defaults$skipTau)
          gene_info$param_init <- verify_initial_params(gene_info$param_init,
                                                         defaults, FALSE)
          
          # Clean up regime names
          if (any(names(gene_info$param_init$theta) != defaults$parsed_alt_tree$unique_regimes)){
            defaults$parsed_alt_tree$unique_regimes <- names(gene_info$param_init$theta)
          }
          if (!'root' %in% defaults$parsed_alt_tree$unique_regimes & defaults$add.root){
            defaults$parsed_alt_tree$unique_regimes <- c('root', defaults$parsed_alt_tree$unique_regimes)
          }
          
          # Run EM
          res1 <- run_em(tree_constants$tree, gene_info$X, gene_info$param_init,
                         tree_constants$edges, defaults,
                         max_iter = 100,
                         tol_rel = 1e-4,
                         tol_abs = 0.01,
                         lambda1 = lambda1,
                         lambda2 = lambda2,
                         diagnose = FALSE,
                         tau_prior_mean = tau_prior_mean,
                         tau_prior_sd = tau_prior_sd,
                         save = FALSE,
                         save_pre = NULL,
                         skipE = M_only,
                         verbose = verbose)
          
          res1$settings <- as.vector(combos[j, ])
          batch_results[[length(batch_results) + 1]] <- res1
        }
        
        # Write batch to disk
        batch_file <- file.path(results_dir, sprintf("batch_%03d.rds", match(batch_idx[1], unlist(batch_indices))))
        saveRDS(batch_results, batch_file)
        
        p()
        
        # Return nothing (results are on disk)
        NULL
        
      }, future.seed = TRUE)
    })

    end_time <- Sys.time()
    log_message(paste0('End time = ', format(end_time, "%H:%M:%S")[[1]], '\n'), logfile, verbose=TRUE)
    


}

####################
# Helpers 

combine_disk_batch <- function(results_dir, tid, batch_size=10, tmpdir='./tmp', keeptemp=FALSE) {
  batch_files <- list.files(results_dir, pattern = "^batch_.*\\.rds$",
                            full.names = TRUE)
  batch_files <- sort(batch_files)
  nbatch_files <- length(batch_files)
  batch_indices <- split(seq_len(nbatch_files), ceiling(seq_len(nbatch_files) / batch_size))

  create_directory_if_not_exists(tmpdir)

  for (batch_ix in batch_indices){
    batch_results <- list()
  
    for (b in batch_ix) {
        batch_data <- readRDS(batch_files[b])
        batch_results <- c(batch_results, batch_data)
    }
    batch_results_df <- data.frame(
          gene = sapply(batch_results, function(x) x$settings[1]),
          regime = sapply(batch_results, function(x) x$settings[2]),
          loglik = sapply(batch_results, function(x) x$final$ll_total),
          alpha = sapply(batch_results, function(x) x$final$alpha),
          sigma = sapply(batch_results, function(x) x$final$sigma),
          objective = sapply(batch_results, function(x) x$final$objective), 
          tau = sapply(batch_results, function(x) x$final$tau), 
          iter = sapply(batch_results, function(x) x$final$iter), 
          converged = sapply(batch_results, function(x) x$final$converge),
          stringsAsFactors = FALSE
    )
    write.csv(batch_results_df, paste0(tmpdir, '/tmp_res_batch_', batch_ix[1], '_', tid, '.csv'))
    saveRDS(batch_results, paste0(tmpdir, '/tmp_res_batch_', batch_ix[1], '_', tid, '.rds'))
  }
  
  tmpfiles <- list.files(tmpdir, pattern = 'tmp_res_batch_.*.csv', full.names=TRUE)
  combined_batches <- do.call(rbind, lapply(tmpfiles, function(i){read.csv(i, row.names=1)}))
  write.csv(combined_batches, sprintf('%s/%s_%s_results_condensed.csv', results_dir, tid, format(Sys.Date(), "%Y%m%d")))

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

create_directory_if_not_exists <- function(directory_path, message = TRUE) {
  if (!dir.exists(directory_path)) {
    dir.create(directory_path, recursive = TRUE)
    if (message) cat("Directory created:", directory_path, "\n")
  } else {
    if (message) cat("Directory already exists:", directory_path, "\n")
  }
}

`%nin%` = Negate(`%in%`)




######## DIAGNOSTICS -- Claude 


check_tree_health <- function(phy) {
  cat(sprintf("N tips: %d\n", length(phy$tip.label)))
  cat(sprintf("Branch lengths: min=%.6f, max=%.6f, mean=%.6f\n",
              min(phy$edge.length), max(phy$edge.length), 
              mean(phy$edge.length)))
  
  V <- ape::vcv.phylo(phy)
  eig <- eigen(V)$values
  
  cat(sprintf("V eigenvalues: min=%.2e, max=%.2e\n", min(eig), max(eig)))
  cat(sprintf("Condition number: %.2e\n", max(eig) / min(eig)))
  
  if (max(eig) / min(eig) > 1e10) {
    cat("WARNING: Tree is ill-conditioned (very large condition number)\n")
    cat("This will cause numerical issues.\n")
  }
}

# In your EM loop, before running EM:
verify_initial_params <- function(paras, defaults, verbose=TRUE) {
  
  if (verbose){
      cat("=== PARAMETER VERIFICATION ===\n")
      cat(sprintf("Alpha: %.8f\n", paras$alpha))
      cat(sprintf("Sigma: %.8f\n", paras$sigma))
      cat(sprintf("Theta mean: %.6f\n", mean(paras$theta)))
  }
  # Check: are they reasonable?
  if (paras$sigma < 0.05) {
    if (verbose) cat("WARNING: Sigma smaller than lower bound, resetting...\n") # want to change this to be an input parameter that is flexible! Should also try with better lower bounds. 
    paras$sigma <- 0.1
  }
  
  if (paras$alpha < 1e-6 && defaults$model != 'BM1') {
    if (verbose) cat("WARNING: Alpha suspiciously small, resetting...\n")
    paras$alpha <- 0.1
  }

  if (any(paras$theta < 1e-10)) {
    paras$theta <- pmax(pmin(paras$theta, 1000), 1e-7)
  }
  
  return(paras)
}


########## TROUBLESHOOTING 
diagnose_residuals <- function(X, e_result) {
  
  residuals <- X - e_result$E_Z
  
  cat("=== RESIDUAL DIAGNOSTICS ===\n")
  cat(sprintf("Mean residual:      %.6f\n", mean(residuals)))
  cat(sprintf("SD of residuals:    %.6f\n", sd(residuals)))
  cat(sprintf("Max abs residual:   %.6f\n", max(abs(residuals))))
  cat(sprintf("Mean |residual|:    %.6f\n", mean(abs(residuals))))
  
  cat("\n=== VARIANCE DECOMPOSITION ===\n")
  resid_var <- mean(residuals^2)
  post_var  <- mean(diag(e_result$Var_Z_given_X))
  total     <- resid_var + post_var
  
  cat(sprintf("Residual variance:   %.6f  (%5.1f%% of total)\n", 
              resid_var, 100 * resid_var / total))
  cat(sprintf("Posterior variance:  %.6f  (%5.1f%% of total)\n", 
              post_var, 100 * post_var / total))
  cat(sprintf("Total (tau^2):       %.6f\n", total))
  cat(sprintf("Implied tau:         %.6f\n", sqrt(total)))
  
  # Check Kalman gain
  cat("\n=== KALMAN GAIN ===\n")
  cat("If K ~ I, E[Z] ~ X and residuals will be near zero\n")
  
  return(invisible(list(
    resid_var = resid_var,
    post_var  = post_var,
    ratio     = resid_var / total
  )))
}

check_kalman_gain <- function(K) {

  
  # Eigenvalues of K tell you how close to identity it is
  eig_K <- eigen(K, symmetric = TRUE)$values
  
  cat(sprintf("Kalman gain eigenvalues: min = %.4f, max = %.4f\n",
              min(eig_K), max(eig_K)))
  
  # If max eigenvalue close to 1: K ~ I, residuals will be tiny
  if (max(eig_K) > 0.95) {
    cat("WARNING: K close to identity matrix\n")
    cat("E[Z] will be very close to X\n")
    cat("Residuals will be near zero\n")
  }
  
  return(eig_K)
}

# FUnction to check whether the likelihoods are oscillating. For now , will implement an early stop if this is the case. 
is_oscillating <- function(x, min_cycles = 3) {
  if (length(x) < (2 * min_cycles + 1)) return(FALSE)
  dx <- diff(x)
  dx <- dx[abs(dx) > 1e-10]
  if (length(dx) < 2) return(FALSE)
  sum(diff(sign(dx)) != 0) >= (2 * min_cycles - 1)
}

