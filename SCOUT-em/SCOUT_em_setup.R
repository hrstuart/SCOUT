######### preprocessing functions 
collectEdges <- function(phy, data, model, Tmax_init, root.age, scaleHeight){
    #A boolean for whether the root theta should be estimated -- default is that it should be.
    if (model %in% c('BM1', 'OU1')){
        ##Begins the construction of the edges matrix -- similar to the ouch format##
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
        phy$node.label <- as.numeric(int.states)
        tip.states <- factor(data[,1])
        data[,1] <- as.numeric(tip.states)
        #Obtain root state and internal node labels
        root.state <- phy$node.label[1]
        int.state <- phy$node.label[-1]
    }
    #New tree matrix to be used for subsetting regimes
    n <- max(phy$edge[,1])
    edges=cbind(c(1:(n-1)), phy$edge, MakeAgeTable(phy, root.age=root.age))
    if(scaleHeight==TRUE){
        edges[,4:5]<-edges[,4:5]/Tmax_init
        root.age <- 1
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

    return (list('phy' = phy,'edges' = edges, 'k' = k, 'root' = root.state, 
        'three.point.params' = list('tot.states' = tot.states, 'int.states' = int.states, 'tip.states' = tip.states)))
}

createIndexMat <- function(model, k, algorithm, get.root.theta){
    index.mat <- matrix(0, 2,k)
    
    if(algorithm == "three.point"){
        Rate.mat <- matrix(1, 3, k)
    }else{
        Rate.mat <- matrix(1, 2, k)
    }

    if (model == "BM1"){
        np <- 1
        index.mat[1,1:k] <- NA
        index.mat[2,1:k] <- 1
                
        param.count <- np+1

        if(algorithm == "three.point"){
            index.mat <- rbind(index.mat, rep(NA,k))
        }

        # Moving to the runner function ... 

        #    index.mat[is.na(index.mat)] <- param.count + 1
        #} else {
         #   index.mat[is.na(index.mat)] <- param.count
        #}
    }
        
    if (model %in% c( "OU1", "OUM") ){
        np <- 2

        index.mat[1,1:k] <- 1
        index.mat[2,1:k] <- 2
        
        if(algorithm == "three.point"){
            max.par.so.far <- max(index.mat)
            if (model == 'OU1') {
                index.mat <- rbind(index.mat, rep(3,k))
            } else {
                index.mat <- rbind(index.mat, (max.par.so.far + 1):(max.par.so.far + k))
            }
        }
        param.count <- np + k # k = 1 if OU1 
                    
        if(get.root.theta == TRUE){
            param.count <- param.count + 1
        }

        index.mat[is.na(index.mat)] <- param.count + 1
    }

    return(list('index.mat' = index.mat, 'param.count' = param.count, 'rate.mat' = Rate.mat))
}

##### Function to initialize starting values. 
initialize_bounds<- function(index.mat, np, algorithm='invert', lb=NULL, ub=NULL){
    if(algorithm == "three.point"){
        if(any(is.na(index.mat[1,]))){
            k.alpha <- 0
            k.sigma <- length(unique(index.mat[2,]))
        } else {
            k.alpha <- length(unique(index.mat[1,]))
            k.sigma <- length(unique(index.mat[2,]))
        }
        if(is.null(lb)){
            lb.alpha <- 1e-9
            lb.sigma <- 1e-9
            lower = c(rep(log(lb.alpha), k.alpha), rep(log(lb.sigma), k.sigma))
            lb <- 1e-9
        }else{
            lower = c(rep(log(lb[1]), k.alpha), rep(log(lb[2]), k.sigma))
            lb <- lb[3] # theta's are added later
        }
        if(is.null(ub)){
            ub.alpha <- 1000
            ub.sigma <- 1000
            upper = c(rep(log(ub.alpha), k.alpha), rep(log(ub.sigma), k.sigma))
            ub <- 1000
        } else {
            upper = c(rep(log(ub[1]), k.alpha), rep(log(ub[2]), k.sigma))
            ub <- ub[3] # theta's are added later
        }
    } else {
        if(is.null(lb)){
            lb <- 1e-9
        }
        if(is.null(ub)){
            ub <- 1000
        }
    lower = rep(log(lb), np)
    upper = rep(log(ub), np)  

    }
    return(list('lower' = lower, 'lb'=lb, 'upper' = upper, 'ub' = ub))
}

initialize_starting <- function(phy, model, x, index.mat, k, bounds, Tmax_init, algorithm, get.root.theta, scaleHeight, starting.vals=NULL, em=NULL, tp.params=NULL){ 
    if(algorithm == "three.point"){
        tot.states <- tp.params$tot.states
        tip.states <- tp.params$tip.states
        int.states <- tp.params$int.states
        shift.point <- tp.params$shiftpoint
        map <- getMapFromNode(phy, tip.states, int.states, shift.point)
        if(scaleHeight==TRUE){
            map <- lapply(map, function(x) x/Tmax_init)
        } 
    }

    if(scaleHeight==TRUE){
        phy$edge.length <- phy$edge.length/Tmax_init
        Tmax <- 1
    } else {
        Tmax <- Tmax_init
    }

    lower <- bounds$lower; upper <- bounds$upper
    lb <- bounds$lb; ub <- bounds$ub
    if(model %in% c("OU1", 'OUM')){
        #Initial value for alpha is just the half life based on the entire length of the tree:
        if(is.null(starting.vals)){
            init.ip <- c(log(2)/Tmax, mean(pic(x, phy)^2))
        } else {
            init.ip <- starting.vals
        }

        if(model=="OU1"){
            ip <- init.ip
        } else {
            ip <- c(rep(init.ip[1],length(unique(index.mat[1,]))), rep(init.ip[2],length(unique(index.mat[2,]))))
        }

        if(algorithm == "three.point"){
            if(model == "OU1"){
                ip <- c(ip, mean(x))
            } else {
                means.by.regime <- with(data, tapply(data[,2], data[,1], mean))
                if(length(means.by.regime) < length(levels(tot.states))){
                    means.by.regime <- rep(mean(data[,2]), length(levels(tot.states)))
                }
                names(means.by.regime) <- NULL
                ip <- c(ip, means.by.regime)
            }

            # NEED TO FIX THIS! 
            lower <- c(lower, rep(log(lb), k))
            upper <- c(upper, rep(log(ub), k))

           # never used this with the 'true' setting so might remove. 
           # if(get.root.theta == TRUE){
           #     ip <- c(ip, means.by.regime[root.state])
           #     lower <- c(lower, log(lb))
           #     upper <- c(upper, log(ub))
           # }
        }
    } else {
        if(is.null(starting.vals)){
            sig <- mean(pic(x, phy)^2)
        }else{
            sig <- starting.vals
        }
        
        ip <- sig

        if(get.root.theta == TRUE){
            if(algorithm == "three.point"){
                ip <- c(ip, mean(x))
                lower <- c(lower, log(lb))
                upper <- c(upper, log(ub))
            }
        }
    }
    return(list('ip' = ip, 'lb' = lower, 'ub' = upper))
}

emitBounds <- function(em, bounds, x){
    lb <- bounds$lb; ub <- bounds$ub

     # Not sure how to handle this in an OUM situation... I wonder if we keep it this way without multiple parameters for multiple thetas? That is what chat gpt said... 
    # The starting values for OU need to be in the log-space but if using the NB or POIS model, 
    # the starting values need to be in the counts space...  and the counts need to be passed to the nloptr function. 
    if (em == 'NORM'){
        # Add mean and sd to init & add x2 params to lower and upper bounds. 
        ip <- c(mean(x), sd(x))
        lower <- c(log(lb), log(lb))
        upper <- c(log(ub), log(ub))

    } else if (em == 'NB'){
        # Add mu and size to init & add x2 params to lower and upper bounds. 
        ip <- c( mean(x), 1)
        lower <- c(log(lb), log(lb))
        upper <- c(log(ub), log(ub))
    } else if (em == 'POIS'){
        # Add lambda to init & add x1 params to lower and upper bounds. 
        ip <- mean(x)
        lower <- log(lb)
        upper <-  log(ub)
    } else {
        stop('Specified model not supported. Try from list: NORM, NB, POIS')
    }

    return(list(ip, lower, upper))
}
  
fitSCOUT.EM <- function(phy, data, 
    model=c("BM1","OU1","OUM"),
    em = NULL, # NORM, NB, POIS
    mode = 'ignore', # on/off for generative mode or not. 
    # I don't think generative mode makes a whole lot of sense at the moment. Probably need to switch to estimating in the latent space. 
    root.age=NULL, 
    scaleHeight=FALSE, 
    root.station=FALSE, 
    get.root.theta=FALSE, 
    shift.point=0.5, 
    algorithm=c("invert", "three.point"), 
    diagn=FALSE, # havent decided what to do with this.
    quiet=FALSE, 
    warn=TRUE, 
    starting.vals=NULL,
    lb = NULL, 
    ub = NULL,
    verbose = FALSE) 
{
    phy <- reorder.phylo(phy, "cladewise")

    if(model == "BM1" & root.station == FALSE) { get.root.theta = TRUE }

    ### Start of true prep: 
    if (!is.null(em) & (em == 'NB' || em == 'POIS')){
        # If a counts model, add an extra column with the raw counts for initial values. 
        data <- data.frame(data[,2], log1p(data[,3]),  data[,3], row.names=data[,1]) # reformat to be state, gex, with tip label as row names 
    } else {
        # assumes data is already log transformed. 
        data <- data.frame(data[,2], data[,3], row.names=data[,1]) # reformat to be state, gex, with tip label as row names 
    }
    
    data <- data[phy$tip.label,]
    tip.states.cp <- factor(data[,1]) # not sure if this is used ever. 

    ntips <- length(phy$tip.label)

    # initialize Tmax 
    Tmax.i <- max(MakeAgeTable(phy, root.age=root.age))

    edges0 <- collectEdges(phy, data, model, Tmax.i, root.age, scaleHeight) # this only uses the states data for the edges, not the phylogenetic data.... 
    index0 <- createIndexMat(model, edges0$k, algorithm, get.root.theta) 

    nump <- index0$param.count - edges0$k
    bounds0 <- initialize_bounds(index0$index.mat, nump, algorithm, lb, ub)

    #This is the first line that requires data. 
    if(algorithm == "three.point"){
        x <- data[,2]
        names(x) <- rownames(data)
        tip.paths <- lapply(1:length(data[,2]), function(x) getPathToRoot(phy, x))
        tpparam <- edges0[['three.point.params']]
        tpparam$shiftpoint <- shift.point
        init_vals <- initialize_starting(phy, model, x, index0$index.mat, edges0$k, bounds0, 
            Tmax.i, algorithm, get.root.theta, scaleHeight, starting.vals, em, tpparam)

    } else {
        #x <- as.matrix(data[,2])
        x <- data[,2]; names(x) <- row.names(data)
        init_vals <- initialize_starting(phy, model, x, index0$index.mat, edges0$k, bounds0, Tmax.i, algorithm, get.root.theta, scaleHeight, starting.vals, em)

    }

    #######################
    # Gather inputs 
    ip <- init_vals$ip
    lower <- init_vals$lb
    upper <- init_vals$ub

    index.mat <- index0$index.mat
    edges <- edges0$edges
    Rate.mat <- index0$rate.mat
    root.state <- edges0$root
    phy <- edges0$phy

    #######################
    np_adjust <- 0 # establish as 0. 
    if (!is.null(em)){
        if(em %in% c('POIS', 'NB')) { i <- data[,3] } else {i <- data[,2] }# if NORM then give the log-transformed otherwise, give the counts. 
        adjust <- emitBounds(em, bounds0, i)

        ip <- c(ip, adjust[[1]])
        lower <- c(lower, adjust[[2]])
        upper <- c(upper, adjust[[3]])

        np_adjust <- length(adjust[[1]])

        xobs <- i; names(xobs) <- row.names(data)
    } else {
        xobs <- x
    }

    #print(head(xobs))

   # print(index.mat)
    if (verbose == TRUE){
        cat('\nStarting Values and Bounds:\n')
        
        init_vals_mat <- matrix(c(ip, lower, upper), nrow=3, byrow=TRUE)
        row.names(init_vals_mat) <- c('init_vals', 'lower', 'upper')
        print(init_vals_mat)

        cat('\nStarting fit...\n')
        print(paste0('Mode = ', mode, ' and Version = ', em))

        print(phy)
    }


    # RUN! 
    res_em<- nloptr(
        x0 = log(ip),
        eval_f = optim_joint,
        opts = list('algorithm' = "NLOPT_LN_SBPLX", 'maxeval' = '1000', 'ftol_rel'=.Machine$double.eps^0.5),
        obs = xobs[phy$tip.label], ### fix this. 
        tree = phy,
        edges = edges, 
        index.mat = index.mat,
        Rate.mat = Rate.mat, 
        version = em, 
        mode = mode, 
        root.age = root.age,
        root.state=root.state,
        get.root.theta = get.root.theta, 
        lb = fix_lower(lower, log(ip)), 
        ub = fix_upper(upper, log(ip))
    )

    loglik = -res_em$objective
    param.count = index0$param.count + np_adjust 

    result <- list(loglik = loglik, 
        AIC = -2*loglik+2*param.count, 
        AICc=-2*loglik+(2*param.count*(ntips/(ntips-param.count-1))), 
        solution = exp(res_em$solution), # alpha, sigma, .... 
        emit.model = em,
        algorithm = algorithm)

    if (verbose == TRUE){
        cat('\nSolution:\n')
        print(exp(res_em$solution))

        cat('\nLog-Likelihood:\n')
        print(result$loglik)

        cat('\nAIC:\n')
        print(result$AIC)
    }

    return(result)
}


  #  opts = list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000", "ftol_rel"=.Machine$double.eps^0.5) )

    