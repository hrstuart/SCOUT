## Pagel's lambda screen: a pre-fit filter for genes with no phylogenetic signal.
##
## A gene whose tip values carry no signal on the tree (lambda ~ 0) is i.i.d. noise as far as any
## OU/BM model is concerned, and fitting SCOUT to it only adds false positives (see the ~8-10% OUM
## floor on pure i.i.d. counts). lambda_screen() scores every gene, plot_lambda_screen() shows the
## distribution so a cutoff can be chosen by eye, and select_lambda_genes() returns the genes to
## pass to SCOUT(genes = ...). SCOUT(lambda_filter = TRUE, ...) runs the same screen internally.
##
## TWO ENGINES, ONE ESTIMATOR. Both maximise phytools' likelihoodLambda over the same 10
## sub-intervals of [0, maxLambda] with optimize(), and test against lambda = 0 with a 1-df LRT,
## i.e. phytools::phylosig(method = 'lambda', test = TRUE).
##
##   fast      For an ULTRAMETRIC tree diag(C) is a constant h, so phytools' transform
##             Cl = lambda*(C - diag(diag(C))) + diag(diag(C)) is Cl = lambda*C + (1 - lambda)*h*I,
##             which shares C's eigenvectors. One eigen(C) per tree makes every likelihood
##             evaluation O(ntips) instead of O(ntips^3). Agrees with phytools to ~1e-11.
##   phylosig  phytools::phylosig per gene. Rebuilds vcv.phylo and solves inside every likelihood
##             evaluation (~300 per gene): ~0.5 s/gene at 142 tips, >600 s/gene at 1719 tips.
##             Used when the tree is not ultrametric, where the fast identity does not hold.
##
## Ported from revision_v2/4a_LUAD_CML/scripts/260920_luad_cml_phylosig.R.

#' Is this trait usable? A constant or non-finite trait has no defined lambda, and must come out as
#' NA rather than as a signal of zero.
#' @noRd
lambda_usable <- function(v) all(is.finite(v)) && is.finite(stats::var(v)) && stats::var(v) > 0

NA_LAMBDA <- c(lambda = NA_real_, lambda_logL = NA_real_, lambda_logL0 = NA_real_, lambda_P = NA_real_)

#' Tip depths of a tree (the diagonal of its BM vcv), without building the vcv.
#' @noRd
tip_depths <- function(phy) ape::node.depth.edgelength(phy)[seq_len(ape::Ntip(phy))]

#' Is the tree ultrametric to a relative tolerance on tip depth?
#' @noRd
is_ultrametric_rel <- function(phy, tol = 1e-6) {
    d <- tip_depths(phy)
    diff(range(d)) <= tol * max(d)
}

#' phytools' maxLambda without phytools: height over the deepest internal node for an ultrametric
#' tree (its `max(nodeHeights[,2]) / max(nodeHeights[,1])`), and 1 otherwise.
#' @noRd
max_lambda <- function(phy, ultrametric) {
    if (!ultrametric) return(1)
    d <- ape::node.depth.edgelength(phy)
    nt <- ape::Ntip(phy)
    max(d) / max(d[nt + seq_len(phy$Nnode)])
}

#' Everything about one tree that no gene changes. Built once per screen.
#' `niter` is phytools' own default number of sub-intervals of [0, maxlam].
#' @noRd
lambda_tree_cache <- function(phy, niter = 10) {
    C <- ape::vcv.phylo(phy)
    n <- nrow(C)
    e <- eigen(C, symmetric = TRUE)
    maxlam <- max_lambda(phy, ultrametric = TRUE)
    list(tips = rownames(C), n = n, h = mean(diag(C)),
         U = e$vectors, d = e$values, u1 = as.numeric(crossprod(e$vectors, rep(1, n))),
         maxlam = maxlam, niter = niter,
         interval = cbind(seq(0, maxlam - maxlam / niter, length.out = niter),
                          seq(maxlam / niter, maxlam, length.out = niter)))
}

#' phytools' likelihoodLambda, evaluated in the eigenbasis of C. `y` is the trait projected onto
#' the eigenvectors. determinant()$modulus is log|det|, hence abs(): this reproduces phytools'
#' arithmetic even where Cl is not positive definite, so the optimiser sees the same surface.
#' @noRd
loglik_lambda <- function(lambda, cc, y) {
    w  <- lambda * cc$d + (1 - lambda) * cc$h          # eigenvalues of Cl
    a  <- sum(cc$u1 * y / w) / sum(cc$u1 * cc$u1 / w)  # GLS root state
    z  <- y - a * cc$u1
    s2 <- sum(z * z / w) / cc$n
    -0.5 * cc$n - 0.5 * cc$n * log(2 * pi) - 0.5 * (cc$n * log(abs(s2)) + sum(log(abs(w))))
}

#' Pagel's lambda + its LRT, same sub-intervals and optimiser as phytools.
#' @noRd
lambda_fast <- function(y, cc) {
    fits <- lapply(seq_len(nrow(cc$interval)), function(i)
        stats::optimize(loglik_lambda, interval = cc$interval[i, ], cc = cc, y = y, maximum = TRUE))
    liks <- vapply(fits, function(f) f$objective, 0)
    res  <- fits[[which(liks == max(liks))[1]]]
    l0   <- loglik_lambda(0, cc, y)
    c(lambda = res$maximum, lambda_logL = res$objective, lambda_logL0 = l0,
      lambda_P = stats::pchisq(2 * (res$objective - l0), df = 1, lower.tail = FALSE))
}

#' Reference engine: phytools::phylosig on one gene.
#' @noRd
lambda_phylosig <- function(phy, v) {
    l <- phytools::phylosig(phy, v, method = 'lambda', test = TRUE)
    c(lambda = unname(l$lambda), lambda_logL = unname(l$logL),
      lambda_logL0 = unname(l$logL0), lambda_P = unname(l$P))
}

#' Read counts and tree the way SCOUT()/formatSCOUT() do, and return the tree plus a
#' tips x genes matrix on the scale SCOUT will fit.
#' @noRd
lambda_inputs <- function(counts, tree, species_key, regimes, blacklist, genes, normalize,
                          logfile, verbose) {
    meta <- if (is.data.frame(counts)) counts else utils::read.csv(counts, row.names = 1)
    if (!'species' %in% colnames(meta)) {
        if (!is.null(species_key) && species_key %in% colnames(meta)) {
            names(meta)[names(meta) == species_key] <- 'species'
        } else {
            stop('Please provide a column named `species` or use species_key to identify it.')
        }
    }
    phy <- if (inherits(tree, 'phylo')) tree else ape::read.tree(tree)

    miss <- setdiff(phy$tip.label, as.character(meta$species))
    if (length(miss)) {
        stop(sprintf('Missing counts for %d leaves in the tree, e.g. %s. Please prune the tree.',
                     length(miss), paste(utils::head(miss, 3), collapse = ', ')))
    }
    # Same edge-length fixes as formatSCOUT, so the screen sees the tree SCOUT fits.
    if (is.null(phy$edge.length)) {
        phy$edge.length <- rep(1, ape::Nedge(phy))
        log_message('lambda_screen: no edge lengths, replacing with 1s.', logfile, verbose)
    } else if (any(phy$edge.length == 0)) {
        phy$edge.length[phy$edge.length == 0] <- 1e-7
    }

    if (!is.null(genes)) {
        gene_cols <- gsub('-', '.', unlist(genes))       # formatSCOUT's name cleaning
        absent <- setdiff(gene_cols, colnames(meta))
        if (length(absent)) {
            log_message(sprintf('lambda_screen: %d requested gene(s) not in the counts, e.g. %s',
                                length(absent), paste(utils::head(absent, 3), collapse = ', ')),
                        logfile, verbose)
        }
        gene_cols <- intersect(gene_cols, colnames(meta))
    } else {
        gene_cols <- setdiff(colnames(meta), c('species', species_key, regimes, 'OU1', blacklist))
    }
    numeric_col <- vapply(meta[gene_cols], is.numeric, logical(1))
    if (any(!numeric_col)) {
        log_message(sprintf('lambda_screen: skipping %d non-numeric column(s): %s',
                            sum(!numeric_col), paste(gene_cols[!numeric_col], collapse = ', ')),
                    logfile, verbose)
        gene_cols <- gene_cols[numeric_col]
    }
    if (!length(gene_cols)) stop('lambda_screen: no numeric gene columns to screen.')

    rownames(meta) <- as.character(meta$species)
    X <- as.matrix(meta[phy$tip.label, gene_cols, drop = FALSE])
    storage.mode(X) <- 'double'
    if (normalize) {
        neg <- colSums(X < 0, na.rm = TRUE) > 0
        if (any(neg)) {
            warning(sprintf(paste0('lambda_screen: %d gene(s) have negative values, so log1p gives NaN ',
                                   'and their lambda will be NA. Already log-normalised? Use normalize = FALSE.'),
                            sum(neg)), call. = FALSE)
        }
        X <- suppressWarnings(log1p(X))
    }
    list(phy = phy, X = X)
}

#' Pagel's lambda screen for phylogenetic signal
#'
#' Scores every gene's Pagel's lambda on the tree, as a filter to run before \code{SCOUT()}: genes
#' with lambda near zero carry no phylogenetic signal and are not appropriate for OU/BM model
#' fitting. Use \code{plot_lambda_screen()} to look at the distribution and choose a cutoff, then
#' \code{select_lambda_genes()} to get the gene list for \code{SCOUT(genes = ...)}. Or set
#' \code{SCOUT(lambda_filter = TRUE, lambda_min = ...)} to screen inside SCOUT.
#'
#' Counts and tree are read exactly as \code{SCOUT()} reads them (same name mangling, same edge
#' fixes, \code{log1p} when \code{normalize = TRUE}), so lambda is measured on the scale SCOUT fits.
#'
#' The estimator is \code{phytools::phylosig(method = 'lambda', test = TRUE)}. On an ultrametric
#' tree an exact eigen-decomposition engine is used (one \code{eigen()} per tree, O(ntips) per
#' likelihood evaluation). If the tree is NOT ultrametric that identity fails, a message says so,
#' and the slow engine (\code{phytools::phylosig} per gene, ~O(ntips^3) per likelihood evaluation)
#' is used; phytools must then be installed. For a non-ultrametric tree maxLambda is 1, as in
#' phytools.
#'
#' \code{lambda_P} is the 1-df chi-square LRT of lambda vs 0. Zero is the boundary of the
#' parameter space, where that test is conservative; read it as an upper bound on the p-value.
#' Estimates that land at the bottom of the search interval (~1e-4) mean "zero".
#'
#' CAUTION: lambda measures BM-like covariance, and a strong OU pull erases it. On the bundled
#' example (alpha = 3, 256 tips) true OU1 genes score lambda ~ 0 and most OUM genes < 0.1, so a
#' lambda cutoff can remove exactly the OU genes SCOUT is meant to find. Inspect the distribution
#' with \code{plot_lambda_screen()} and check the cutoff against the alpha*H range of interest.
#'
#' @param counts path to a counts CSV (read with \code{read.csv(row.names = 1)}, as SCOUT does) or
#'   a data.frame.
#' @param tree path to a newick file or a \code{phylo} object. Every tip must be in the counts.
#' @param species_key name of the column holding tip labels, if it is not \code{species}.
#' @param regimes regime columns (as passed to SCOUT); never screened as genes.
#' @param blacklist other non-gene columns to ignore.
#' @param genes optional subset of genes to screen; default is every numeric non-reserved column.
#' @param normalize apply \code{log1p} before scoring; use the same value as in \code{SCOUT()}.
#' @param engine \code{'auto'} (fast if ultrametric, else phylosig), \code{'fast'} (error if the
#'   tree is not ultrametric) or \code{'phylosig'} (always the reference engine).
#' @param cores number of cores (\code{parallel::mclapply}) for the per-gene loop.
#' @param lambda_min,p_max optional cutoffs; when either is given a logical \code{pass} column is
#'   added (lambda >= lambda_min and lambda_P <= p_max; NA never passes).
#' @param ultrametric_tol relative tolerance on tip depths for the ultrametric check.
#' @param logfile,verbose logging, as in \code{SCOUT()}.
#' @return A data.frame of class \code{scout_lambda_screen}, one row per gene: gene_name, lambda,
#'   lambda_P, lambda_logL, lambda_logL0, mean_expr, var_expr, frac_zero, n_nonzero, engine, ntips
#'   (and pass). Attributes \code{engine}, \code{ultrametric}, \code{normalize}, \code{max_lambda}.
#' @examples
#' \dontrun{
#' scr <- lambda_screen('counts.csv', 'tree.nwk', species_key = 'cellBC', regimes = 'OU4')
#' plot_lambda_screen(scr, lambda_min = 0.2)
#' keep <- select_lambda_genes(scr, lambda_min = 0.2)
#' SCOUT('counts.csv', 'tree.nwk', 'out', regimes = c('BM1', 'OU1', 'OU4'),
#'       species_key = 'cellBC', genes = keep)
#' }
#' @export
lambda_screen <- function(counts, tree,
    species_key = 'species',
    regimes = NULL,
    blacklist = NULL,
    genes = NULL,
    normalize = TRUE,
    engine = c('auto', 'fast', 'phylosig'),
    cores = 1,
    lambda_min = NULL,
    p_max = NULL,
    ultrametric_tol = 1e-6,
    logfile = NULL,
    verbose = TRUE) {

    engine <- match.arg(engine)
    inp <- lambda_inputs(counts, tree, species_key, regimes, blacklist, genes, normalize,
                         logfile, verbose)
    phy <- inp$phy; X <- inp$X
    n <- ape::Ntip(phy)
    ultra <- is_ultrametric_rel(phy, ultrametric_tol)

    if (engine == 'auto') engine <- if (ultra) 'fast' else 'phylosig'
    if (engine == 'fast' && !ultra) {
        stop("Tree is not ultrametric; the fast engine assumes constant tip depth. ",
             "Use engine = 'auto' or 'phylosig'.")
    }
    if (engine == 'phylosig') {
        if (!requireNamespace('phytools', quietly = TRUE)) {
            stop("The slow lambda engine needs the phytools package: install.packages('phytools')")
        }
        if (!ultra) {
            # ~1.2e-7 * n^3 s/gene, from 0.51 s at 142 tips and ~600 s at 1719 tips.
            est <- 1.2e-7 * n^3 * ncol(X) / max(1, cores) / 60
            log_message(sprintf(paste0('lambda_screen: tree is not ultrametric (tip depths %.4g-%.4g); ',
                                       'using the SLOW phytools::phylosig engine. %d genes x %d tips ',
                                       'on %d core(s): roughly %.1f min.'),
                                min(tip_depths(phy)), max(tip_depths(phy)), ncol(X), n, cores, est),
                        logfile, verbose = TRUE)
        }
    }
    log_message(sprintf('lambda_screen: %d genes, %d tips, engine %s', ncol(X), n, engine),
                logfile, verbose)

    if (engine == 'fast') {
        cc <- lambda_tree_cache(phy)
        Y  <- crossprod(cc$U, X[cc$tips, , drop = FALSE])     # every gene projected at once
        one <- function(j) {
            if (!lambda_usable(X[, j])) return(NA_LAMBDA)
            tryCatch(lambda_fast(Y[, j], cc), error = function(e) NA_LAMBDA)
        }
        maxlam <- cc$maxlam
    } else {
        one <- function(j) {
            v <- stats::setNames(X[, j], phy$tip.label)
            if (!lambda_usable(v)) return(NA_LAMBDA)
            tryCatch(lambda_phylosig(phy, v), error = function(e) NA_LAMBDA)
        }
        maxlam <- if (ultra) max_lambda(phy, TRUE) else 1
    }

    res <- if (cores > 1) parallel::mclapply(seq_len(ncol(X)), one, mc.cores = cores)
           else lapply(seq_len(ncol(X)), one)
    # A worker killed mid-job returns an error/NULL; never let that drop a gene silently.
    bad <- which(!vapply(res, function(r) is.numeric(r) && length(r) == 4, logical(1)))
    if (length(bad)) res[bad] <- lapply(bad, one)
    S <- do.call(rbind, res)

    out <- data.frame(gene_name = colnames(X),
                      lambda = S[, 'lambda'], lambda_P = S[, 'lambda_P'],
                      lambda_logL = S[, 'lambda_logL'], lambda_logL0 = S[, 'lambda_logL0'],
                      mean_expr = colMeans(X), var_expr = apply(X, 2, stats::var),
                      frac_zero = colMeans(X == 0), n_nonzero = colSums(X != 0),
                      engine = engine, ntips = n,
                      stringsAsFactors = FALSE, row.names = NULL)
    if (!is.null(lambda_min) || !is.null(p_max)) out$pass <- lambda_pass(out, lambda_min, p_max)

    attr(out, 'engine') <- engine
    attr(out, 'ultrametric') <- ultra
    attr(out, 'normalize') <- normalize
    attr(out, 'max_lambda') <- maxlam
    class(out) <- c('scout_lambda_screen', 'data.frame')

    log_message(sprintf('lambda_screen: median lambda %.3f | %d NA | %.1f%% with lambda_P < 0.05',
                        stats::median(out$lambda, na.rm = TRUE), sum(is.na(out$lambda)),
                        100 * mean(out$lambda_P < 0.05, na.rm = TRUE)), logfile, verbose)
    out
}

#' @noRd
lambda_pass <- function(screen, lambda_min = NULL, p_max = NULL) {
    keep <- rep(TRUE, nrow(screen))
    if (!is.null(lambda_min)) keep <- keep & screen$lambda >= lambda_min
    if (!is.null(p_max))      keep <- keep & screen$lambda_P <= p_max
    keep & !is.na(keep)
}

#' Genes that pass a lambda cutoff
#'
#' @param screen output of \code{lambda_screen()}.
#' @param lambda_min keep genes with lambda >= lambda_min.
#' @param p_max keep genes with lambda_P <= p_max (LRT of lambda vs 0).
#' @return Character vector of gene names, for \code{SCOUT(genes = ...)}. Genes with NA lambda
#'   (constant or non-finite) never pass.
#' @export
select_lambda_genes <- function(screen, lambda_min = NULL, p_max = NULL) {
    if (is.null(lambda_min) && is.null(p_max)) stop('Give lambda_min and/or p_max.')
    screen$gene_name[lambda_pass(screen, lambda_min, p_max)]
}

#' Plot the distribution of Pagel's lambda from a screen
#'
#' Histogram of per-gene lambda with the chosen cutoff marked and the number of genes it keeps in
#' the title, to pick a cutoff before running \code{SCOUT()}.
#'
#' @param screen output of \code{lambda_screen()}.
#' @param lambda_min,p_max candidate cutoffs, as in \code{select_lambda_genes()}.
#' @param breaks number of histogram bins.
#' @param main optional title; by default reports the genes kept.
#' @param ... passed to \code{graphics::hist()}.
#' @return Invisibly, a list with n_total, n_na, n_kept and the kept genes.
#' @export
plot_lambda_screen <- function(screen, lambda_min = NULL, p_max = NULL, breaks = 50,
                               main = NULL, ...) {
    lam <- screen$lambda[!is.na(screen$lambda)]
    if (!length(lam)) stop('No finite lambda values to plot.')
    top <- max(attr(screen, 'max_lambda') %||% 1, lam)
    kept <- if (is.null(lambda_min) && is.null(p_max)) screen$gene_name[!is.na(screen$lambda)]
            else select_lambda_genes(screen, lambda_min, p_max)
    if (is.null(main)) {
        crit <- c(if (!is.null(lambda_min)) sprintf('lambda >= %g', lambda_min),
                  if (!is.null(p_max)) sprintf('P <= %g', p_max))
        main <- if (length(crit)) sprintf('%d / %d genes kept (%s)', length(kept), nrow(screen),
                                          paste(crit, collapse = ', '))
                else sprintf("Pagel's lambda, %d genes", nrow(screen))
    }
    graphics::hist(lam, breaks = seq(0, top, length.out = breaks + 1), main = main,
                   xlab = "Pagel's lambda", ylab = 'Genes', col = 'grey80', border = 'white', ...)
    if (!is.null(lambda_min)) graphics::abline(v = lambda_min, col = 'firebrick', lty = 2, lwd = 2)
    n_na <- sum(is.na(screen$lambda))
    if (n_na) graphics::mtext(sprintf('%d gene(s) with NA lambda not shown', n_na), side = 3,
                              line = 0.2, cex = 0.8)
    invisible(list(n_total = nrow(screen), n_na = n_na, n_kept = length(kept), genes = kept))
}
