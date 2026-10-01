source('/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/1_tree_robustness/scripts/scout_perturb_lib.R')
source('/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/5_real_tree_robustness/scripts/gene_support.R')

make_perturbed_trees <- function(tr, cfg, outdir) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    phy <- tr$phy
    ARM_BASE <- c(shuffle = 10000, nni = 20000, collapse = 30000)
    has_pg <- requireNamespace('phangorn', quietly = TRUE)

    grid <- rbind(
        data.frame(arm = 'reference', intensity = 0, replicate = 1L, stringsAsFactors = FALSE),
        expand.grid(arm = cfg$arms, intensity = cfg$intensities,
                    replicate = seq_len(cfg$n_rep), stringsAsFactors = FALSE))
    if (!is.null(cfg$shuffle_at)){
        grid <- rbind(grid, data.frame(arm = 'shuffle', intensity = cfg$shuffle_at,
                   replicate = seq_len(cfg$n_rep), stringsAsFactors = FALSE))
    }

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


#######################################################

perturb.config <- list(arms = c('nni', 'collapse'),
               n_rep = 30, 
               intensities = c(0.05, 0.10, 0.15, 0.2), 
               shuffle_at = NULL, 
               seed_tree = 1000, 
               cores = 46)


m5k.obj = list(
        name  = 'm5k_v2a',
        treepath  = paste0('/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/analysis/quinn_2021_science/',
                       'A549_M5K/v3/m5k_lg13_tree_hybrid_priors_pruned_convexML_resolvedMulti_cln_root.nwk'), 
        countspath = "/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/analysis/quinn_2021_science/A549_M5K/v3//m5k_lg13_tree_hybrid_priors_HVG_counts.csv",
        outpath = "/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/5_real_tree_robustness/output_support/luad_v2a",
       # genes =  "/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/5_real_tree_robustness/data/m5k_lg13_convexML_SCOUT_sample100_genes.csv"
        genes =  "/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/5_real_tree_robustness/data/260907_fullQC_m5k_lg13_convexML_SCOUT_sample90_even_genes.csv"

        
    )

#######################################################

m5k.obj$phy <- ape::read.tree(m5k.obj$treepath)
m5k_perturb = make_perturbed_trees(m5k.obj, perturb.config, m5k.obj$outpath)
cat(sprintf('%d trees written to %s\n\n', nrow(m5k_perturb$manifest),
            file.path(m5k.obj$outpath, 'perturbed_trees')))

m5k_perturb_summary <- m5k_perturb$manifest %>%
    group_by(arm, intensity) %>%
    summarise(n = n(), k = first(k), rf_mean = mean(rf_dist), rf_sd = sd(rf_dist),
              height = mean(height), depth_cv = mean(depth_cv), .groups = 'drop') %>%
    mutate(across(c(rf_mean, rf_sd, height, depth_cv), ~ signif(.x, 4)))


write.csv(m5k_perturb$manifest, paste0(m5k.obj$outpath, '/', m5k.obj$name, 'manifest.csv'))
write.csv(m5k_perturb_summary, paste0(m5k.obj$outpath, '/', m5k.obj$name, 'manifest_summary.csv'))

#######################################################

sampled_meta <- read.csv(m5k.obj$countspath)
genes_to_keep <- read.table(m5k.obj$genes)[,1]

mask_cols <- c('cellBC', 'OU4')
sampled_meta <- sampled_meta[, c(mask_cols, genes_to_keep)]
names(sampled_meta)[which(names(sampled_meta) == 'cellBC')] <- 'species'

#######################################################

obs_cfg <- list(meta = sampled_meta, gene_set = genes_to_keep, name = m5k.obj$name, regimes = c('BM1', 'OU1', 'OU4'),
    outdir = m5k.obj$outpath)

res <- run_gene_support(obs_cfg, m5k_perturb, layer = 'expression', cores = perturb.config$cores)
cat(sprintf('\ntotal compute: %.1f min over %d fits of %d genes\n', sum(res$meta$elapsed_min), nrow(res$meta), length(res$genes)))
saveRDS(res, paste0(m5k.obj$outpath, '/', m5k.obj$name, '_gene_support.rds'))


if (have_obs('res')) {
    SUPP <- gene_support(res)
    COND <- support_by_condition(SUPP)
    CONF <- gene_confidence(SUPP)
}


write.csv(SUPP, file.path(m5k.obj$outpath,  paste0(m5k.obj$name, '_observed_gene_support.csv')), row.names = FALSE)
write.csv(COND, file.path(m5k.obj$outpath, paste0(m5k.obj$name, '_observed_support_by_condition.csv')), row.names = FALSE)
write.csv(CONF, file.path(m5k.obj$outpath, paste0(m5k.obj$name, '_observed_gene_confidence.csv')), row.names = FALSE)
cat('written to', m5k.obj$outpath, '\n')



