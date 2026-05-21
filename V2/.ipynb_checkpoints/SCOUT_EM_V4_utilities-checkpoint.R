## For one iter 

####### SUPPORT FUNCTIONS #######
# This assumes a stationary theta. -- double check that this makes sense
# Might need to make appropriate changes for BM1. 
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


transformPhy.new <- function(phy, map, Alpha_j, Sigma_j, tip.paths=NULL){
  nTip <- length(phy$tip.label)
  RootAge <- max(branching.times(phy))
  NodeAges <- branching.times(phy)[phy$edge[,1] - nTip]
  # ModMap <- Map <- map
  D <- V_Tilde <- numeric(dim(phy$edge)[1])
  for(i in 1:dim(phy$edge)[1]){
    # evaluate the map for this particular edge and calculate the tipward variance
    NodeAge_i <- NodeAges[i]
    DistRoot_i <- RootAge - NodeAge_i
    Map_i <- map[[i]]
    # the age of epoch j starts at the node age
    Dist_rootward <- DistRoot_i
    z <- w <- v <- 0
    for(j in 1:length(Map_i)){
      # distance the root of epoch j starts at the node distance and ends at node dist + epoch length
      Dist_tipward <- Dist_rootward + Map_i[j]
      # the length of the epoch is scaled by the alpha parameter of that epoch -- alpha and sigma defined in function.  
      # calculate the descendant distance from the root based on a fixed root distribution
      tmp.w <- Alpha_j * (Dist_tipward - Dist_rootward)
      tmp.v <- Sigma_j * (exp(2 * Alpha_j * Dist_tipward) - exp(2 * Alpha_j * Dist_rootward))/2/Alpha_j
      v <- v + tmp.v
      w <- w + tmp.w
      # ModMap[[i]][j] <- tmp.v
      # The new distance from nodes
      Dist_rootward <- Dist_tipward
    }
    V_Tilde[i] <- v 
    D[i] <- w
  }
  phy$edge.length <- V_Tilde
  # calculates the diagonal matrix for each tip i
  DiagWt <- numeric(nTip)
  names(DiagWt) <- phy$tip.label
  if(is.null(tip.paths)){
    for(i in 1:nTip){
      DiagWt[i] <- exp(-sum(D[getPathToRoot(phy, i)]))
    }
  }else{
    for(i in 1:nTip){
      DiagWt[i] <- exp(-sum(D[tip.paths[[i]]]))
    }
  }
  obj <- list(tree = phy, diag = DiagWt)
  return(obj)
}

### WEIGHT MAT FUNCTION
#Weight matrix generator taken from Butler and King (2004) and modified to allow multiple alpha parameters
#written by Jeremy M. Beaulieu
#Taken from Beaulieu 2012 and unmodified to only need one alpha parameter. 

weight.mat<-function(phy, edges, alpha, root.state,  root.age=NULL,  assume.station=TRUE, shift.point=0.5){
    age.table <- MakeAgeTable(phy, root.age=root.age)
    Tmax <- max(age.table)
    n <- max(phy$edge[,1])
    ntips <- length(phy$tip.label)
    
    if(is.null(root.state)) {
        root.state<-which(edges[dim(edges)[1],]==1)-5
        edges <- edges[-1*dim(edges)[1],]
    }

    mm <- dim(edges)
    k <- length(6:mm[2])

    pp <- prop.part(phy)
    edges <- edges
    nodevar.root.tot <- rep(0,max(edges[,3]))
    nodevar.k <- rep(0,max(edges[,3]))
    W <- matrix(0, ntips, k)

    for(j in 1:k){
        oldregime <- root.state
        n.cov.root.tot = matrix(0, n, 1)
        n.cov.k <- matrix(0, n, 1)
        #Weights for each species per regime
        for(i in 1:length(edges[,1])){
            anc <- edges[i, 2]
            oldtime <- edges[i,4]
            newtime <- edges[i,5]
            
            if(anc%in%edges[,3]){
                start <- which(edges[,3]==anc)
                oldregime <- which(edges[start,6:(k+5)]==1)
            }
            else{
                #For the root:
                oldregime <- root.state
            }
            newregime <- which(edges[i,6:(k+5)]==1)
            if(oldregime==newregime){
                if(oldregime == j){
                    nodevar.root.tot[i] <- -alpha*(newtime-oldtime)
                    nodevar.k[i] <- exp(alpha*newtime)-exp(alpha*oldtime)
                }
                else{
                    nodevar.root.tot[i] <- -alpha*(newtime-oldtime)
                    nodevar.k[i] <- 0
                }
            }
            else{
                shifttime <- newtime-((newtime-oldtime) * shift.point)
                epoch1a <- -alpha*(shifttime-oldtime)
                epoch1b <- exp(alpha*shifttime)-exp(alpha*oldtime)
                oldtime <- shifttime
                epoch2a <- -alpha*(newtime-oldtime)
                epoch2b <- exp(alpha*newtime)-exp(alpha*oldtime)
                nodevar.root.tot[i] <- epoch1a + epoch2a
                    
                if(oldregime==j){
                    nodevar.k[i] <- epoch1b
                }
                if(newregime==j){
                    nodevar.k[i] <- epoch2b
                }
                if(!newregime==j && !oldregime==j){
                    nodevar.k[i] <- 0
                }
            }
        }
        n.cov.k[edges[i,3],] <- nodevar.k[i]
        n.cov.root.tot[edges[i,3],] <- nodevar.root.tot[i]
    }
    w.k <- mat.gen(phy, n.cov.k, pp)
    w.root.tot <- mat.gen(phy, n.cov.root.tot, pp)

    W[1:(ntips),j] <- exp(diag(w.root.tot)) * diag(w.k)


    if(assume.station == TRUE){
        w.root.tot <- mat.gen(phy, n.cov.root.tot, pp)
        W[,root.state] <- W[,root.state] + exp(diag(w.root.tot))
    }else{
        w.root.tot <- mat.gen(phy, n.cov.root.tot, pp)
        W <- cbind(exp(diag(w.root.tot)), W)
    }

    #Restandardizes W so that the rows sum to 1 -- Generalized. Will reduce to the simpler model if assuming 1 alpha parameter, but when alpha varies by regime they will sum to 1 (though proportionally should be ok).
    W <- W/rowSums(W)
    return(W)
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

#### THREE POINT SUPPORT FUNCTIONS 
#OU functions for using three-point algorithm of Ho and Ane 2013

#written by James D. Boyko


## takes a node based reconstruction and returns a map (identical to a map from simmap)
getMapFromNode <- function(phy, tipstates, nodestates, shift.point){
  Map <- vector("list", dim(phy$edge)[1])
  Data <- c(tipstates, nodestates)
  NodeStates <- cbind(Data[phy$edge[,1]], Data[phy$edge[,2]])
  for(i in 1:dim(phy$edge)[1]){
    from <- as.character(NodeStates[i,1])
    to <- as.character(NodeStates[i,2])
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

# data(tworegime)
# phy <- tree
# tipstates <- trait[,2]
# nodestates <- phy$node.label
# shift.point <- 0.5
# getMapFromNode(phy, round(runif(length(phy$tip.label), 1, 2)), nodestates, 0.5)
# getMapFromNode(phy, tipstates, nodestates, 0)


# gets the path from a vertex to the root as an index of the edge matrix
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


# transforms the phylogeny based on a set of parameters and a simmap
transformPhy <- function(phy, map, pars, tip.paths = NULL) {
  #phy must be of class simmap
  nTip <- length(phy$tip.label)

  #ensure we traverse edges from root toward tips so Ato[ancestor] is ready
  phy <- ape::reorder.phylo(phy, "cladewise")

  nEdge <- nrow(phy$edge)
  nNodeTot <- nTip + phy$Nnode

  Ato <- numeric(nNodeTot)
  root <- nTip + 1
  Ato[root] <- 0
  D <- V_Tilde <- numeric(nEdge)

  for (i in 1:nEdge) {
    anc <- phy$edge[i, 1]
    des <- phy$edge[i, 2]
    Map_i <- map[[i]]
    Acur <- Ato[anc] 
    w <- 0
    v <- 0
    for (j in 1:length(Map_i)) {
      dt <- as.numeric(Map_i[j])
      Sigma_j <- pars[, 2][match(names(Map_i)[j], rownames(pars))]
      Alpha_j <- pars[, 3][match(names(Map_i)[j], rownames(pars))]
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
  obj <- list(tree = phy, diag = DiagWt)
  return(obj)
}


# transforms the phylogeny based on a set of parameters and a simmap
transformPhy.old <- function(phy, map, pars, tip.paths=NULL){
  # phy must be of class simmap
  nTip <- length(phy$tip.label)
  RootAge <- max(branching.times(phy))
  NodeAges <- branching.times(phy)[phy$edge[,1] - nTip]
  # ModMap <- Map <- map
  D <- V_Tilde <- numeric(dim(phy$edge)[1])
  for(i in 1:dim(phy$edge)[1]){
  # evaluate the map for this particular edge and calculate the tipward variance
  NodeAge_i <- NodeAges[i]
  DistRoot_i <- RootAge - NodeAge_i
  Map_i <- map[[i]]
  # the age of epoch j starts at the node age
  Dist_rootward <- DistRoot_i
  z <- w <- v <- 0
  for(j in 1:length(Map_i)){
    # distance the root of epoch j starts at the node distance and ends at node dist + epoch length
    Dist_tipward <- Dist_rootward + Map_i[j]
    # the length of the epoch is scaled by the alpha parameter of that epoch
    Sigma_j <- pars[,2][match(names(Map_i)[j], rownames(pars))]
    Alpha_j <- pars[,3][match(names(Map_i)[j], rownames(pars))]
    # calculate the descendant distance from the root based on a fixed root distribution
    tmp.w <- Alpha_j * (Dist_tipward - Dist_rootward)
    tmp.v <- Sigma_j * (exp(2 * Alpha_j * Dist_tipward) - exp(2 * Alpha_j * Dist_rootward))/2/Alpha_j
    v <- v + tmp.v
    w <- w + tmp.w
    # ModMap[[i]][j] <- tmp.v
    # The new distance from nodes
    Dist_rootward <- Dist_tipward
  }
  V_Tilde[i] <- v
  D[i] <- w
  }
  phy$edge.length <- V_Tilde
  # calculates the diagonal matrix for each tip i
  DiagWt <- numeric(nTip)
  names(DiagWt) <- phy$tip.label
  if(is.null(tip.paths)){
  for(i in 1:nTip){
    DiagWt[i] <- exp(-sum(D[getPathToRoot(phy, i)]))
  }
  }else{
  for(i in 1:nTip){
    DiagWt[i] <- exp(-sum(D[tip.paths[[i]]]))
  }
  }
  obj <- list(tree = phy, diag = DiagWt)
  return(obj)
}



getOULik <- function(phy, y, X, pars){
  # transform the phylogeny based on params
  tre <- transformPhy(phy, pars)
  # use the transformed phylogeny for the three point algorithm
  comp <- three.point.compute(tre$tree, y, X, tre$diag)
  # calculate the likelihood
  lik <- -as.numeric(Ntip(phy) * log(2 * pi) + comp$logd + (comp$PP - 2 * comp$QP + comp$QQ))/2
  return(lik)
}

##### OU VAR-COVAR FUNCTIONS 
#OU variance-covariance matrix generator

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
            vcv <- exp(-2*alpha[1]*max(root.age)) * vcv2
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


##Matrix generating function taken from vcv.phylo in ape:
mat.gen <- function(phy,piece.wise,pp){
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


######## OTHER UTILITY FUNCTIONS 

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

# transforms the phylogeny based on a set of parameters and a simmap
transformPhy <- function(phy, map, pars, tip.paths=NULL){
  # phy must be of class simmap
  nTip <- length(phy$tip.label)
  RootAge <- max(branching.times(phy))
  NodeAges <- branching.times(phy)[phy$edge[,1] - nTip]
  # ModMap <- Map <- map
  D <- V_Tilde <- numeric(dim(phy$edge)[1])
  for(i in 1:dim(phy$edge)[1]){
    # evaluate the map for this particular edge and calculate the tipward variance
    NodeAge_i <- NodeAges[i]
    DistRoot_i <- RootAge - NodeAge_i
    Map_i <- map[[i]]
    # the age of epoch j starts at the node age
    Dist_rootward <- DistRoot_i
    z <- w <- v <- 0
    for(j in 1:length(Map_i)){
      # distance the root of epoch j starts at the node distance and ends at node dist + epoch length
      Dist_tipward <- Dist_rootward + Map_i[j]
      # the length of the epoch is scaled by the alpha parameter of that epoch
      Sigma_j <- pars[,2][match(names(Map_i)[j], rownames(pars))]
      Alpha_j <- pars[,3][match(names(Map_i)[j], rownames(pars))]
      # calculate the descendant distance from the root based on a fixed root distribution
      tmp.w <- Alpha_j * (Dist_tipward - Dist_rootward)
      tmp.v <- Sigma_j * (exp(2 * Alpha_j * Dist_tipward) - exp(2 * Alpha_j * Dist_rootward))/2/Alpha_j
      v <- v + tmp.v
      w <- w + tmp.w
      # ModMap[[i]][j] <- tmp.v
      # The new distance from nodes
      Dist_rootward <- Dist_tipward
    }
    V_Tilde[i] <- v 
    D[i] <- w
  }
  phy$edge.length <- V_Tilde
  # calculates the diagonal matrix for each tip i
  DiagWt <- numeric(nTip)
  names(DiagWt) <- phy$tip.label
  if(is.null(tip.paths)){
    for(i in 1:nTip){
      DiagWt[i] <- exp(-sum(D[getPathToRoot(phy, i)]))
    }
  }else{
    for(i in 1:nTip){
      DiagWt[i] <- exp(-sum(D[tip.paths[[i]]]))
    }
  }
  obj <- list(tree = phy, diag = DiagWt)
  return(obj)
}

####=======#### 

fix_upper <- function(ub, x0) {
  if(any(ub<x0)) {
    warning("Some upper bounds are less than initial values, adjusting")
    ub[which(ub<x0)] <- x0[which(ub<x0)]
  } 
  return(ub)
}

fix_lower <- function(lb, x0) {
  if(any(lb>x0)) {
    warning("Some lower bounds are greater than initial values, adjusting")
    lb[which(lb>x0)] <- x0[which(lb>x0)]
  } 
  return(lb)
}


preprocess <- function(inputs, reg, gene_name, 
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

######## DIAGNOSTICS
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
verify_initial_params <- function(paras, phy, edges, defaults) {
  
  cat("=== PARAMETER VERIFICATION ===\n")
  cat(sprintf("Alpha: %.8f\n", paras$alpha))
  cat(sprintf("Sigma: %.8f\n", paras$sigma))
  cat(sprintf("Theta mean: %.6f\n", mean(paras$theta)))
  
  # Check: are they reasonable?
  if (paras$sigma < 0.001) {
    cat("WARNING: Sigma suspiciously small, resetting...\n")
    paras$sigma <- 1.0
  }
  
  if (paras$alpha < 1e-6 && defaults$model != 'BM1') {
    cat("WARNING: Alpha suspiciously small, resetting...\n")
    paras$alpha <- 0.1
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