# 01_make_trees.r -- Goal 1, branch-length baseline, step 01: generate every perturbed tree
# (shuffle @ 100%, nni and collapse @ 10/50/90%) and write a manifest. Seeds and perturbed edges
# match the published unit-branch-length run, so the two pair replicate-by-replicate.
#
# Usage:
#   Rscript 01_make_trees.r                 # optional env: SCOUT_ARMS, SCOUT_INTENSITIES, SCOUT_NREP

# SCOUT_LIB may be a colon-separated path LIST, so a private build (e.g. the support_clip
# SCOUT) can be prepended while its dependencies still resolve from the shared library.
.libPaths(strsplit(Sys.getenv('SCOUT_LIB',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/software/R/R-4.4.2/library'), ':')[[1]])
ROOT <- Sys.getenv('SCOUT_ROOT',
    '/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/1_tree_robustness')
source(file.path(ROOT, 'scripts', 'scout_perturb_lib.R'))
source(file.path(ROOT, 'scripts', '1_tree_perturb_bl', 'bl_paths.R'))

P <- bl_paths(ROOT)
DATADIR <- Sys.getenv('SCOUT_DATADIR', P$datadir)
# A SEPARATE tree directory. The filenames (nni_i010_rep001.nwk) are byte-identical to the published
# ones, so sharing data/1_tree_perturb/perturbed_trees/ would overwrite 105 of the 630 trees the
# published results were produced from.
TREEDIR <- file.path(DATADIR, 'perturbed_trees')
PREFIX  <- Sys.getenv('SCOUT_BL_PREFIX', '260826_simtree_256cells_truebl')
N_REP   <- as.integer(Sys.getenv('SCOUT_NREP', '15'))
INTENSITIES <- as.numeric(strsplit(Sys.getenv('SCOUT_INTENSITIES', '0.10,0.50,0.90'), ',')[[1]])
ARMS    <- strsplit(Sys.getenv('SCOUT_ARMS', 'shuffle,nni,collapse'), ',')[[1]]

guard_empty(TREEDIR, 'perturbed tree directory')

# Rank each intensity held in the FULL published grid, so a subset keeps the seeds it had there.
# shuffle is a single 100% arm and had no multiplier: its base seed is 10000 flat.
SEED_RANK <- c('0.01' = 1, '0.05' = 2, '0.1' = 3, '0.25' = 4, '0.5' = 5, '0.75' = 6, '0.9' = 7)
ARM_BASE  <- c(shuffle = 10000, nni = 20000, collapse = 30000)

base_seed_for <- function(arm, intensity) {
    if (arm == 'shuffle') return(unname(ARM_BASE['shuffle']))
    key <- as.character(intensity)
    if (!key %in% names(SEED_RANK))
        stop(sprintf('intensity %s is not in the published grid, so it has no matched seed. ',
                     key), 'Add it to SEED_RANK with a rank no other intensity uses.')
    unname(ARM_BASE[arm] + SEED_RANK[[key]] * 1000)
}

treefile <- file.path(DATADIR, paste0(PREFIX, '.nwk'))
if (!file.exists(treefile))
    stop(sprintf('baseline tree %s not found -- run 00_simulate_baseline.r first.', treefile))

# READ FROM DISK. See the header -- regenerating in memory breaks edge-matching.
phy <- read.tree(treefile)
NT  <- length(phy$tip.label)
cat(sprintf('baseline: %d tips, %d internal edges, height %.4f, ultrametric %s\n',
            NT, n_internal_edges(phy), max(node.depth.edgelength(phy)[seq_len(NT)]),
            is.ultrametric(phy)))
cat(sprintf('arms: %s | intensities: %s%% | replicates: %d\n',
            paste(ARMS, collapse = ','), paste(INTENSITIES * 100, collapse = ','), N_REP))

plan <- rbind(
    data.frame(arm = 'shuffle',  intensity = 1.00, stringsAsFactors = FALSE),
    data.frame(arm = 'nni',      intensity = INTENSITIES, stringsAsFactors = FALSE),
    data.frame(arm = 'collapse', intensity = INTENSITIES, stringsAsFactors = FALSE))
plan <- plan[plan$arm %in% ARMS, , drop = FALSE]
stopifnot(nrow(plan) > 0)
plan$base_seed <- mapply(base_seed_for, plan$arm, plan$intensity)

# A seed collision between two conditions would make them share a perturbation while looking
# independent in the manifest. Cheap to rule out, expensive to discover later.
stopifnot(!any(duplicated(plan$base_seed)))

man <- do.call(rbind, lapply(seq_len(nrow(plan)), function(i) {
    p <- plan[i, ]
    cat(sprintf('  %-9s @ %5.1f%%  seed %d  n=%d ... ', p$arm, p$intensity * 100, p$base_seed, N_REP))
    m <- perturb_tree_set(phy, type = p$arm, intensity = p$intensity, n = N_REP,
                          randseed = p$base_seed, outdir = TREEDIR)
    cat(sprintf('k=%d  mean RF=%.1f\n', m$k[1], mean(m$rf_dist, na.rm = TRUE)))
    m
}))

manfile <- file.path(DATADIR, 'perturbed_tree_manifest_bl.csv')
write.csv(man, manfile, row.names = FALSE)

stopifnot(nrow(man) == nrow(plan) * N_REP, all(file.exists(man$tree_file)),
          !any(duplicated(man$tree_file)))

# ---------------------------------------------------------------------------------------------
# The edge-matching assertion.
#
# This is the check the entire cross-run comparison rests on. RF is a topology-only distance, and
# the published trees were perturbed from the SAME topology with the SAME seeds, so if the
# perturbations really are edge-matched then every k and every per-replicate rf_dist must equal the
# published value to the integer. Anything else means the RNG stream moved -- a different phangorn
# or ape, a changed grid, a regenerated instead of round-tripped baseline -- and 04_collate.r's
# bl-vs-v2 delta would be comparing conditions that are not actually matched.
#
# Fail loudly. A quietly incomparable result is worse than no result.
# ---------------------------------------------------------------------------------------------
v2_files <- file.path(ROOT, 'data', '1_tree_perturb',
                      c('perturbed_tree_manifest.csv', 'perturbed_tree_manifest_hi.csv'))
v2_files <- v2_files[file.exists(v2_files)]

if (!length(v2_files)) {
    cat('\nWARNING: no published tree manifest found -- edge-matching to the v2 run is UNVERIFIED.\n')
} else {
    cols <- c('arm', 'intensity', 'replicate', 'k', 'rf_dist')
    v2   <- do.call(rbind, lapply(v2_files, function(f) read.csv(f, stringsAsFactors = FALSE)[, cols]))
    key  <- function(d) paste(d$arm, sprintf('%.4f', d$intensity), d$replicate, sep = '|')
    v2$key <- key(v2); man$key <- key(man)

    shared <- intersect(man$key, v2$key)
    a <- man[match(shared, man$key), ]
    b <- v2 [match(shared, v2$key ), ]

    bad_k  <- which(a$k != b$k)
    bad_rf <- which(a$rf_dist != b$rf_dist)

    cat(sprintf('\n-- edge-matching against the published v2 trees --\n'))
    cat(sprintf('%d of %d conditions comparable | k mismatches %d | RF mismatches %d\n',
                length(shared), nrow(man), length(bad_k), length(bad_rf)))

    if (length(bad_k) || length(bad_rf)) {
        show <- head(unique(c(bad_k, bad_rf)), 10)
        print(data.frame(arm = a$arm[show], intensity = a$intensity[show], replicate = a$replicate[show],
                         k_bl = a$k[show], k_v2 = b$k[show],
                         rf_bl = a$rf_dist[show], rf_v2 = b$rf_dist[show]), row.names = FALSE)
        stop(sprintf('%d k and %d RF mismatches against the published trees. The perturbations are ',
                     length(bad_k), length(bad_rf)),
             'NOT edge-matched, so the bl-vs-v2 comparison would be invalid. ',
             'Check the phangorn/ape versions and that the baseline was READ from its .nwk.')
    }
    if (length(shared) < nrow(man))
        cat(sprintf('note: %d condition(s) have no published counterpart and are unverified\n',
                    nrow(man) - length(shared)))
    cat('every shared condition reproduces the published k and RF exactly -- edge-matched.\n')
}

man$key <- NULL
write.csv(man, manfile, row.names = FALSE)

cat(sprintf('\n%d trees written to %s\n', nrow(man), TREEDIR))
cat(sprintf('manifest: %s\n\n', manfile))
print(aggregate(cbind(k, rf_dist) ~ arm + intensity, data = man,
                FUN = function(x) round(mean(x), 1)), row.names = FALSE)
