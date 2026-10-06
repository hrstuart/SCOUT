# celegans_rpl_scaled_bl_masking.R
# Held-out dropout validation of SCOUT on the C. elegans random-precise-lineage tree (branch
# lengths scaled to unit height). Masks clades and size-matched random cell sets, refits SCOUT-EM on
# the retained cells, and scores predictions of the masked cells against their observed values.
#
# Masked-cell prediction (SCOUT::runSCOUT.dropout) is based on mvMORPH::estim()
# (Clavel, Escarguel & Merceron 2015, Methods Ecol. Evol. 6:1311-1319).
#
# Usage:
#   Rscript celegans_rpl_scaled_bl_masking.R    # inputs, gene/clade counts and cores set in the config block

.libPaths(c('/dartfs/rc/lab/M/McKennaLab/projects/hannah/software/R/R-4.4.2/library', .libPaths()))

library(SCOUT)
library(progressr)
library(future.apply)
library(corpcor)
library(paleotree)
library(nloptr)
library(dplyr)
library(stringr)
library(ape)

base    <- '/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/3_celegans_rpl'
adir    <- file.path(base, 'loo_output_2')
logfile <- file.path(adir, 'celegans_rpl_loo_2.log')
dir.create(adir, showWarnings = FALSE, recursive = TRUE)

N_GENES  <- 100    # sampled from the top 500 by variance
N_CLADES <- 20     # stratified from the eligible set
CLADE_MIN <- 10
CLADE_MAX <- 50
CORES    <- 48
SEED     <- 42

REGIMES   <- c('BM1', 'OU1', 'OUalt', 'OUl')
BLACKLIST <- c('cell', 'OUcl', 'OUct')

############### Inputs -- first row of the samplesheet (the one flagged loo = TRUE) ###############
cele.obj <- list(
        name  = 'cele_bl',
        label = 'C. elegans lineage (rpl subset0)',
        treepath  = paste0('/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/', 'revision_v2/3_celegans_rpl/data/cele_cell_lineage_rpl_subset8.nwk'),
        countspath = "/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/analysis/celegans/random_precise_lineage/rpl3/0_data/celegans_random_precise_lineage_subset_0_mincells_20_metadata.csv",
        outpath = "/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/5_real_tree_robustness/output_support/celegans_rpl",
        genes = '/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/5_real_tree_robustness/data/cele_branch_lengths_SCOUT_sample100_genes.csv' # path to a list of gene names to use. 
    )


#samplesheet <- read.csv(file.path(base, '260805_rpl_resampled_min20_samplesheet.csv'), row.names = 1)
treefile   <- cele.obj$treepath #samplesheet[1, 'tree']
countsfile <- cele.obj$countspath #samplesheet[1, 'data']
testid     <- paste0(cele.obj$name, '_LOO') #paste0(samplesheet[1, 'testid'], '_LOO')

cat(sprintf('tree   : %s\ncounts : %s\ntestid : %s\n', treefile, countsfile, testid))

meta <- read.csv(countsfile, row.names = 1)
phy  <- read.tree(treefile)

# manually scaling the celegans tree before running the perturbation 

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


Tmax.i <- max(MakeAgeTable(phy, root.age=NULL))
phy$edge.length <- phy$edge.length/Tmax.i

############### Gene subset: 100 sampled from the top 500 by variance ###############

# Cost is linear in the gene count (~84 CPU-sec per gene x regime fit at this tree size), so the
# full 3,980-gene matrix is out of reach. 82.7% of the matrix is zero and the median gene variance
# is 0.071, so the low-variance tail carries almost no signal to recover anyway.
gene_all <- setdiff(colnames(meta), c('species', BLACKLIST, REGIMES))
gvar     <- apply(meta[, gene_all], 2, var)
hvg500   <- names(sort(gvar, decreasing = TRUE))[1:500]

set.seed(SEED)
genes <- sample(hvg500, N_GENES)
write.csv(data.frame(gene = genes, variance = gvar[genes]),
          file.path(adir, sprintf('%s_genes.csv', testid)), row.names = FALSE)

############### Mask set: 20 clades stratified from the eligible set ###############

# Build the tree exactly as runSCOUT.dropout will, so the node numbers below refer to the same
# object. infer_nodes = FALSE because formatSCOUT re-derives node labels per regime.
states <- dropout_states(meta, c('OUalt', 'OUl'), 'species')
phy    <- prep_dropout_tree(phy, states, infer_nodes = FALSE)

cand <- eligible_clades(phy, CLADE_MIN, CLADE_MAX, states = states, exclude_root_children = TRUE)
cat(sprintf('%d eligible clades of %d-%d tips (sizes %s)\n', nrow(cand), CLADE_MIN, CLADE_MAX,
            paste(range(cand$size), collapse = '-')))

# Keep every large clade -- they are scarce -- and spread the rest across the size range so the
# mask-size gradient survives the subsample.
big   <- cand$node[cand$size >= 25]
small <- cand[cand$size < 25, ]
small <- small[order(small$size), ]
nbins <- N_CLADES - length(big)
bins  <- cut(seq_len(nrow(small)), breaks = nbins, labels = FALSE)
set.seed(SEED)
picked_small <- sapply(seq_len(nbins), function(b) {
    idx <- which(bins == b)
    small$node[idx[sample(length(idx), 1)]]
})

mask_nodes <- sort(c(big, picked_small))
sizes <- cand$size[match(mask_nodes, cand$node)]
cat(sprintf('Masking %d clades, sizes %s (%.1f%%-%.1f%% of the tree)\n', length(mask_nodes),
            paste(sizes, collapse = ','),
            100 * min(sizes) / Ntip(phy), 100 * max(sizes) / Ntip(phy)))
write.csv(data.frame(node = mask_nodes, size = sizes),
          file.path(adir, sprintf('%s_masked_clades.csv', testid)), row.names = FALSE)

############### Run ###############

res <- runSCOUT.dropout(
    tree        = phy,
    results_dir = adir,
    data        = meta,
    species_key = 'species',
    blacklist   = BLACKLIST,
    normalize   = FALSE,
    genes       = genes,
    mask_states = c('OUalt', 'OUl'),
    mask_nodes  = mask_nodes,
    regimes     = REGIMES,
    arms        = c('clade', 'random'),
    keep_predictions = 'ref_model',
    lambda1 = 0.2,
    lambda2 = 0.2,
    fixed_root = FALSE,
    tau_prior_mean = 0.2,
    tau_prior_sd = 0.1,
    randseed = SEED,
    testid  = testid,
    cores   = CORES,
    logfile = logfile,
    verbose = TRUE
)

############### Headline ###############

# is_selected marks the regime the MASKED fit chose by AIC, so model selection never saw the
# held-out cells.
headline <- res$metrics %>%
    filter(is_selected) %>%
    group_by(arm, predictor) %>%
    summarise(n = n(),
              rmse = mean(rmse, na.rm = TRUE),
              mae = mean(mae, na.rm = TRUE),
              abs_mean_error = mean(abs(mean_error), na.rm = TRUE),
              pearson = mean(pearson, na.rm = TRUE),
              coverage95 = mean(coverage, na.rm = TRUE),
              coverage95_obs = mean(coverage_obs, na.rm = TRUE),
              .groups = 'drop') %>%
    arrange(arm, rmse)

print(as.data.frame(headline))
write.csv(headline, file.path(adir, sprintf('%s_headline.csv', testid)), row.names = FALSE)

agreement <- res$model_select %>%
    group_by(arm) %>%
    summarise(agreement_with_full = mean(best_masked == best_full, na.rm = TRUE),
              .groups = 'drop')
print(as.data.frame(agreement))

cat('Done.\n')
