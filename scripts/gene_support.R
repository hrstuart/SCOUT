# gene_support.R -- gene-level support for SCOUT model calls on real trees: fits SCOUT on a set of
# perturbed trees (make_perturbed_trees, run_gene_support) and summarises how stable each gene's
# call is across perturbations (gene_support, support_by_condition, gene_confidence).
#
# Usage:
#   source('scout_perturb_lib.R'); source('gene_support.R')
#   res  <- run_gene_support(obs_cfg, tree_set, layer = 'expression', cores = 48)
#   SUPP <- gene_support(res)

gene_set_tag <- function(genes) {
    if (requireNamespace('digest', quietly = TRUE)) return(substr(digest::digest(sort(genes)), 1, 8))
    h <- 5381
    for (ch in utf8ToInt(paste(sort(genes), collapse = ','))) h <- (h * 31 + ch) %% 1000000007
    sprintf('h%08d', h)
}

### Tree perturbation baseline 
make_perturbed_trees <- function(tr, cfg, outdir) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    phy <- tr$phy
    ARM_BASE <- c(shuffle = 10000, nni = 20000, collapse = 30000)
    has_pg <- requireNamespace('phangorn', quietly = TRUE)

    grid <- rbind(
        data.frame(arm = 'reference', intensity = 0, replicate = 1L, stringsAsFactors = FALSE),
        expand.grid(arm = cfg$arms, intensity = cfg$intensities,
                    replicate = seq_len(cfg$n_rep), stringsAsFactors = FALSE),
        data.frame(arm = 'shuffle', intensity = cfg$shuffle_at,
                   replicate = seq_len(cfg$n_rep), stringsAsFactors = FALSE))

    trees <- vector('list', nrow(grid))
    rows  <- vector('list', nrow(grid))
    for (i in seq_len(nrow(grid))) {
        g <- grid[i, ]
        if (g$arm == 'reference') {
            pt <- resolve_and_floor(phy); k <- 0L; seed <- NA_integer_
        } else {
            seed <- unname(ARM_BASE[g$arm]) + cfg$seed_tree +
                    as.integer(round(g$intensity * 100)) * 10L + g$replicate
            pt <- perturb_tree(phy, type = g$arm, intensity = g$intensity, randseed = seed)
            k  <- attr(pt, 'perturb_k')
        }
        stem <- sprintf('%s_%s_i%03d_rep%03d', tr$name, g$arm,
                        round(g$intensity * 100), g$replicate)
        f <- file.path(outdir, paste0(stem, '.nwk'))
        ape::write.tree(pt, f)
        d <- ape::node.depth.edgelength(pt)[seq_along(pt$tip.label)]
        trees[[i]] <- pt
        rows[[i]] <- data.frame(
            tree = tr$name, arm = g$arm, intensity = g$intensity, replicate = g$replicate,
            k = k, seed = seed, stem = stem, tree_file = f,
            rf_dist = if (has_pg) suppressWarnings(as.numeric(phangorn::RF.dist(phy, pt, normalize = TRUE))) else NA_real_,
            height = max(d), depth_cv = stats::sd(d) / mean(d),
            n_internal = n_internal_edges(pt), stringsAsFactors = FALSE)
    }
    man <- do.call(rbind, rows)
    stopifnot(!any(duplicated(man$stem)), all(file.exists(man$tree_file)))
    # every arm with intensity > 0 must actually have moved the tree
    moved <- man$arm != 'reference' & man$rf_dist == 0
    if (any(moved & !is.na(man$rf_dist)))
        warning(sprintf('%d perturbed tree(s) have RF = 0 to the reference -- a no-op arm.', sum(moved)))
    list(manifest = man, trees = trees)
}


### Gene Support for Read Trees
run_gene_support <- function(tr, tree_set, cores = 1, layer = 'counts',
                             force = FALSE, verbose = TRUE) {
    
    dat <- tr$meta 
    stopifnot(layer %in% c('counts', 'expression'))
    normalize <- layer == 'counts'      # 'expression' is assumed already on a fittable scale
    tag  <- gene_set_tag(tr$gene_set)
    ckpt <- file.path(tr$outdir, sprintf('support_%s_%s', layer, tag))
    dir.create(ckpt, recursive = TRUE, showWarnings = FALSE)

    man <- tree_set$manifest
    stopifnot(nrow(man) == length(tree_set$trees))
    n <- nrow(man); t_start <- Sys.time()

    out <- vector('list', n)
    for (i in seq_len(n)) {
        row <- man[i, ]
        f <- file.path(ckpt, sprintf('%s.rds', row$stem))
        if (!force && file.exists(f)) {
            out[[i]] <- readRDS(f)
            if (verbose) cat(sprintf('[%2d/%2d] %-28s cached  %d genes\n', i, n, row$stem,
                                     nrow(out[[i]]$selected)))
            next
        }
        phy <- tree_set$trees[[i]]
        stopifnot(setequal(phy$tip.label, dat$species))

        dataset_id <- sprintf('%s|%s|%s|%g|%03d', tr$name, layer, row$arm, row$intensity, row$replicate)
        t0 <- Sys.time()
        res <- score_classification(phy, dat, dataset_id = dataset_id, regimes = tr$regimes,
                                    normalize = normalize, cores = cores, outdir = tempdir())
        el <- as.numeric(difftime(Sys.time(), t0, units = 'mins'))

        sel <- res$selected
        sel$truth <- NULL                       # there is none; see 5.5
        sel <- cbind(row[rep(1, nrow(sel)), c('tree', 'arm', 'intensity', 'replicate', 'k','rf_dist', 'height', 'depth_cv')],
                     sel, row.names = NULL)
        sel$layer <- layer
        rec <- list(meta = cbind(row[, c('tree', 'arm', 'intensity', 'replicate', 'k', 'rf_dist',
                                         'height', 'depth_cv')], layer = layer, n_genes = nrow(sel), elapsed_min = round(el, 2)), selected = sel, annotated = res$annotated)
        saveRDS(rec, f)
        out[[i]] <- rec
        if (verbose) {
            eta <- as.numeric(difftime(Sys.time(), t_start, units = 'mins')) / i * (n - i)
            cat(sprintf('[%2d/%2d] %-28s %d genes  %.1f min  (eta %.0f min)\n',
                        i, n, row$stem, nrow(sel), el, eta))
        }
    }

    selected <- purrr::map_dfr(out, 'selected')
    meta     <- purrr::map_dfr(out, 'meta')

    # The gene set MUST be identical across fits or per-gene support has a moving denominator.
    per_fit <- split(selected$gene_name, interaction(selected$arm, selected$intensity,selected$replicate, drop = TRUE))
    ref_set <- sort(per_fit[[1]])
    bad <- names(per_fit)[!vapply(per_fit, function(g) identical(sort(g), ref_set), logical(1))]
    if (length(bad))
        stop(sprintf('the fitted gene set differs between conditions (%s). Per-gene support would be computed over a moving denominator.',
                     paste(utils::head(bad, 3), collapse = ', ')))

    list(selected = selected, annotated = res$annotated, meta = meta, genes = tr$gene_set, tag = tag, ckpt = ckpt, fits = out)
}



#' Per gene x condition support, both flavours (see 5.5).
gene_support <- function(res, gene_map = NULL, drop_padded = TRUE) {
    sel <- res$selected
    if (drop_padded && length(res$pick$padded))
        sel <- sel %>% filter(!gene_name %in% res$pick$padded)

    ref <- sel %>% filter(arm == 'reference') %>%
        select(gene_name, ref_call = model, ref_AIC = AIC)
    if (!nrow(ref)) stop('no reference fit in the results -- support is defined relative to it.')

    supp <- sel %>% filter(arm != 'reference') %>%
        left_join(ref, by = 'gene_name') %>%
        group_by(tree, layer, gene_name, ref_call, arm, intensity) %>%
        summarise(n_rep       = n(),
                  support_ref = mean(model == ref_call),
                  modal_call  = names(which.max(table(model)))[1],
                  support_mode = max(table(model)) / n(),
                  n_distinct_calls = dplyr::n_distinct(model),
                  rf_dist     = mean(rf_dist), .groups = 'drop')

    if (!is.null(gene_map))
        supp <- supp %>% left_join(gene_map, by = c('gene_name' = 'gene'))
    supp %>% arrange(gene_name, arm, intensity)
}


#' Collapse genes to a mean +/- SE per condition -- the section-2 degradation curve's analogue.
support_by_condition <- function(supp) {
    flr <- supp %>% filter(arm == 'shuffle') %>% pull(support_ref) %>% mean()
    supp %>%
        group_by(tree, layer, arm, intensity) %>%
        summarise(n_genes = dplyr::n_distinct(gene_name), n_rep = max(n_rep),
                  support = mean(support_ref), se = stats::sd(support_ref) / sqrt(n()),
                  support_mode = mean(support_mode),
                  frac_genes_full = mean(support_ref == 1),
                  rf_dist = mean(rf_dist), .groups = 'drop') %>%
        mutate(floor = flr,
               # 1 = as stable as the reference is with itself, 0 = at the no-signal floor
               retained = if (is.finite(1 - flr) && abs(1 - flr) > 1e-9) (support - flr) / (1 - flr)
                          else NA_real_) %>%
        arrange(arm, intensity)
}


#' One row per gene: the headline annotation-confidence table.
#'
#' `support_gradient` pools the nni and collapse arms -- the honest single number for "how much do I
#' trust this call given that the tree is a reconstruction". `support_shuffle` is that gene's own
#' no-signal floor, which is NOT 1/3 and differs between genes.
gene_confidence <- function(supp) {
    grad <- supp %>% filter(arm %in% c('nni', 'collapse'))
    shuf <- supp %>% filter(arm == 'shuffle') %>%
        select(gene_name, support_shuffle = support_ref, shuffle_call = modal_call)
    worst <- grad %>% group_by(gene_name) %>%
        slice_min(support_ref, n = 1, with_ties = FALSE) %>%
        transmute(gene_name, worst_condition = sprintf('%s %.0f%%', arm, 100 * intensity),
                  support_worst = support_ref) %>% ungroup()

    grad %>%
        group_by(tree, gene_name, ref_call) %>%
        summarise(support_gradient = mean(support_ref),
                  support_mild  = mean(support_ref[intensity <= min(intensity)]),
                  support_harsh = mean(support_ref[intensity >= max(intensity)]),
                  .groups = 'drop') %>%
        left_join(worst, by = 'gene_name') %>%
        left_join(shuf,  by = 'gene_name') %>%
        mutate(above_floor = support_gradient - support_shuffle,
               confidence = cut(support_gradient, c(-Inf, 0.5, 0.7, 0.9, Inf),
                                labels = c('unstable', 'weak', 'moderate', 'strong'))) %>%
        arrange(desc(support_gradient))
}




have_obs <- function(...) {
    objs <- c(...)
    miss <- objs[!vapply(objs, exists, logical(1), where = globalenv())]
    if (length(miss)) {
        cat(sprintf('-- skipped: %s not built yet. Set CFG_OBS$expr and run 5.7 / 5.8 first.\n',
                    paste(miss, collapse = ', ')))
        return(FALSE)
    }
    TRUE
}
