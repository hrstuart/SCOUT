# gene_support_cele_bl_runner.R -- gene-support run on the C. elegans rpl lineage tree (scaled to
# unit height): builds nni/collapse/shuffle perturbed trees, refits SCOUT on each, and writes
# per-gene support, support-by-condition and confidence CSVs.
#
# Usage:
#   Rscript gene_support_cele_bl_runner.R    # paths, perturbation grid and cores set in cele.obj / perturb.config

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

MakeAgeTable <- function(phy, root.age=NULL){
    if(is.null(root.age)){
        node.ages <- paleotree::dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))
    }else{
        node.ages <- paleotree::dateNodes(phy, rootAge=root.age)
    }
    max.age <- max(node.ages)
    table.ages <- matrix(0, dim(phy$edge)[1], 2)
    for(row.index in 1:dim(phy$edge)[1]){
        table.ages[row.index,1] <- max.age - node.ages[phy$edge[row.index,1]]
        table.ages[row.index,2] <- max.age - node.ages[phy$edge[row.index,2]]
    }
    return(table.ages)
}



#######################################################

perturb.config <- list(arms = c('nni', 'collapse'),
               n_rep = 30, 
               intensities = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 0.75), 
               shuffle_at = 1.00, 
               seed_tree = 2000, 
               cores = 72)


cele.obj <- list(
        name  = 'cele_bl_v2',
        label = 'C. elegans lineage (rpl subset0)',
        treepath  = paste0('/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/',
                       'revision_v2/3_celegans_rpl/data/cele_cell_lineage_rpl_subset8.nwk'),
        countspath = "/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/analysis/celegans/random_precise_lineage/rpl3/0_data/celegans_random_precise_lineage_subset_0_mincells_20_metadata.csv",
        outpath = "/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/5_real_tree_robustness/output_support/celegans_rpl_v2",
        genes = '/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/5_real_tree_robustness/data/cele_bl_tree_fullQC_SCOUT_evensample120_genes.csv' # path to a list of gene names to use. 
    )


#######################################################
# manually scaling the celegans tree before running the perturbation 
cele.obj$phy <- ape::read.tree(cele.obj$treepath)
Tmax.i <- max(MakeAgeTable(cele.obj$phy, root.age=NULL))
cele.obj$phy$edge.length <- cele.obj$phy$edge.length/Tmax.i

cele_perturb = make_perturbed_trees(cele.obj, perturb.config, cele.obj$outpath)
cat(sprintf('%d trees written to %s\n\n', nrow(cele_perturb$manifest),
            file.path(cele.obj$outpath, 'perturbed_trees')))

cele_perturb_summary <- cele_perturb$manifest %>%
    group_by(arm, intensity) %>%
    summarise(n = n(), k = first(k), rf_mean = mean(rf_dist), rf_sd = sd(rf_dist),
              height = mean(height), depth_cv = mean(depth_cv), .groups = 'drop') %>%
    mutate(across(c(rf_mean, rf_sd, height, depth_cv), ~ signif(.x, 4)))


write.csv(cele_perturb$manifest, paste0(cele.obj$outpath, '/', cele.obj$name, '_manifest.csv'))
write.csv(cele_perturb_summary, paste0(cele.obj$outpath, '/', cele.obj$name, '_manifest_summary.csv'))

#######################################################

sampled_meta <- read.csv(cele.obj$countspath)
genes_to_keep <- read.table(cele.obj$genes)[,1]

mask_cols <- c('species', 'OUalt', 'OUl') # modify to match blacklist. 
sampled_meta <- sampled_meta[, c(mask_cols, genes_to_keep)]

#######################################################

obs_cfg <- list(meta = sampled_meta, gene_set = genes_to_keep, name = cele.obj$name, regimes = c('BM1', 'OU1', 'OUalt', 'OUl'),
    outdir = cele.obj$outpath)

res <- run_gene_support(obs_cfg, cele_perturb, layer = 'expression', cores = perturb.config$cores)
cat(sprintf('\ntotal compute: %.1f min over %d fits of %d genes\n', sum(res$meta$elapsed_min), nrow(res$meta), length(res$genes)))
saveRDS(res, paste0(cele.obj$outpath, '/', cele.obj$name, '_gene_support.rds'))


if (have_obs('res')) {
    SUPP <- gene_support(res)
    COND <- support_by_condition(SUPP)
    CONF <- gene_confidence(SUPP)
}


write.csv(SUPP, file.path(cele.obj$outpath,  paste0(cele.obj$name, '_observed_gene_support.csv')), row.names = FALSE)
write.csv(COND, file.path(cele.obj$outpath, paste0(cele.obj$name, '_observed_support_by_condition.csv')), row.names = FALSE)
write.csv(CONF, file.path(cele.obj$outpath, paste0(cele.obj$name, '_observed_gene_confidence.csv')), row.names = FALSE)
cat('written to', cele.obj$outpath, '\n')



