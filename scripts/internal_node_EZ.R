#!/usr/bin/env Rscript
# =============================================================================
# internal_node_EZ.R
#
# Posterior mean and sd of the latent SCOUT trait, E[Z | X], at every node of the
# tree (tips and internal nodes), using each gene's best-fit SCOUT parameters.
# Reuses SCOUT's internals; does not modify the package.
#
# Node-wise conditional estimation is reimplemented here based on mvMORPH::estim()
# (Clavel, Escarguel & Merceron 2015, Methods Ecol. Evol. 6:1311-1319).
#
# Usage:
#   Rscript internal_node_EZ.R --gene acdh.7 [--check-scaling] [--ref <results.rds>]
#   Rscript internal_node_EZ.R --all --internal-only --out all_genes_EZ.csv.gz
#   Rscript internal_node_EZ.R --help                 # all options
#
# ---------------------------------------------------------------------------
# THE MATH
# ---------------------------------------------------------------------------
# SCOUT's model (see R/SCOUT_EM.R::e_step) is
#
#     Z  ~ N(W theta, V)                      latent trait at the tips
#     X  =  Z + eps,  eps ~ N(0, tau^2 I)     observed tip values ("tip fog")
#
# so the E-step returns the posterior at the TIPS:
#
#     E[Z_tip | X]   = W theta + V (V + tau^2 I)^-1 (X - W theta)
#     Var[Z_tip | X] = V - V (V + tau^2 I)^-1 V
#
# Z is a Gaussian process defined everywhere on the tree, not just at the tips.
# Write Z_A for the trait at ALL nodes (tips 1..n, internal n+1..n+m, in ape
# numbering). Jointly,
#
#     [ Z_A ]      ( [ W_A theta ]  [ V_AA   V_AT            ] )
#     [ X   ] ~ N  ( [ W_T theta ], [ V_TA   V_TT + tau^2 I  ] )
#
# because X is observed only at tips and the fog is independent of Z. Gaussian
# conditioning then gives, for every node at once,
#
#     E[Z_A | X]   = W_A theta + V_AT S^-1 (X - W_T theta),   S = V_TT + tau^2 I
#     Var[Z_A | X] = V_AA - V_AT S^-1 V_TA
#
# Restricted to the tip rows this is *exactly* e_step(), which is the
# consistency check the script runs below.
#
# The only new ingredient is V and W evaluated at internal nodes. SCOUT's
# compute_VCV() and compute_W_matrix() are already written against a small
# summary object (defaults$parsed_alt_tree) whose only tree inputs are
#
#     leaf_dists      root->node distance t_i            (per "leaf")
#     shared_lengths  root->MRCA(i,j) distance s_ij      (per pair)
#     regime_paths    regime segments along root->node   (per "leaf")
#     unique_regimes, root.state
#
# None of that is tip-specific. So we rebuild the same object over ALL nodes
# and hand it to the unmodified SCOUT functions -- no re-derivation of the OU
# covariance, no second implementation to keep in sync:
#
#     V_ij = (sigma/2alpha) exp(-alpha d_ij) [ * (1 - exp(-2 alpha s_ij)) ]
#     d_ij = t_i + t_j - 2 s_ij
#
# which is valid for any pair of points on the tree, internal nodes included.
# This tree is NOT ultrametric; the parameterisation handles that, being
# written in terms of t_i and s_ij rather than a shared tree height.
#
# =============================================================================

suppressWarnings(suppressMessages({
  .default_lib <- "/dartfs/rc/lab/M/McKennaLab/projects/hannah/software/R/R-4.4.2/library"
  if (dir.exists(.default_lib)) .libPaths(c(.default_lib, .libPaths()))
  library(ape)
}))

# ---------------------------------------------------------------------------
# defaults -- point at the C. elegans RPL branch-length-scaled run
# ---------------------------------------------------------------------------
BASE <- "/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/3_celegans_rpl"

DEFAULTS <- list(
  gene       = NULL,
  all        = FALSE,                       # process every QC-passing gene
  qc_column  = "QC_flag",                   # filter column for --all
  qc_value   = "pass",                      # value that counts as passing
  limit      = NULL,                        # cap the gene list (testing)
  internal_only = FALSE,                    # drop tip rows from the output
  verify     = 2L,                          # genes per regime to check vs e_step()
  best_fit   = file.path(BASE, "analysis/pseudoembryo/260910_scout_EM_celegans_bl_scaled_best_fit_filtered.csv"),
  params_dir = file.path(BASE, "output_bl_scaled"),
  tid        = NULL,                        # default: taken from best_fit$dataset
  tree       = file.path(BASE, "data/cele_cell_lineage_rpl_subset8.nwk"),
  counts     = "/dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/analysis/celegans/random_precise_lineage/rpl3/0_data/celegans_random_precise_lineage_subset_0_mincells_20_metadata.csv",
  ref        = NULL,                        # optional cross-check object
  out        = NULL,
  scout_src  = NULL,                        # source R/*.R from here instead of library(SCOUT)

  # --- run settings; these must match the settings the fit was produced under -
  # The producing runner (3_celegans_rpl/scripts/260827_cele_branch_length_runner.R)
  # called SCOUT(..., scale = TRUE, normalize = FALSE,
  #               regimes = c('BM1','OU1','OUalt','OUl'),
  #               blacklist = c('cell','OUcl','OUct'), ...).
  # `scale` is not a formal argument of SCOUT(); R partial-matches it to
  # `scale_tree`, so that run had scale_tree = TRUE -> scaleHeight = TRUE (tree
  # height rescaled to 1). Verified against the stored results object:
  # defaults$scaleHeight == TRUE and leaf_dists ~ 0.54-0.68, not the raw tree
  # height of 834. --check-scaling re-confirms it for any other dataset.
  scale_tree  = TRUE,
  normalize   = FALSE,
  fixed_root  = FALSE,
  infer_anc   = "ape",
  species_key = "species",
  regimes     = c("BM1", "OU1", "OUalt", "OUl"),
  blacklist   = c("cell", "OUcl", "OUct")
)

USAGE <- "
Usage: Rscript internal_node_EZ.R (--gene <GENE> | --all) [options]

  --gene GENE          process a single gene
  --all                process every QC-passing gene in the best-fit table,
                       streamed into ONE output file
  --qc-column COL      QC filter column for --all     (default QC_flag)
  --qc-value VAL       value that counts as passing   (default pass)
  --limit N            process at most N genes (testing)
  --internal-only      omit tip rows from the output (halves the file)
  --verify N           genes per regime cross-checked against e_step() in
                       --all mode (default 2; 0 disables, -1 checks all)

  --best-fit FILE      best_fit_filtered.csv  [celegans RPL default]
  --params-dir DIR     dir holding <tid>_all_genes_<regime>_parameters.csv
  --tid ID             dataset id (default: taken from best_fit 'dataset')
  --tree FILE          newick tree used for the fit
  --counts FILE        counts/metadata csv used for the fit
  --ref FILE           results .rds (or one extracted element) to cross-check
  --out FILE           output csv; .gz extension writes gzipped
                       (default single: ./<gene>_<regime>_node_EZ.csv
                        default --all : ./<tid>_best_fit_node_EZ.csv.gz)
  --scout-src DIR      source R/*.R from this dir instead of library(SCOUT)
  --scale-tree BOOL    scaleHeight used in the fit (default TRUE)
  --normalize BOOL     log1p normalisation used in the fit (default FALSE)
  --fixed-root BOOL    fixed_root used in the fit (default FALSE)
  --species-key COL    tip-id column in the counts csv (default species)
  --regimes LIST       comma-separated regime columns, in the order the fit
                       used them  (default BM1,OU1,OUalt,OUl)
  --blacklist LIST     comma-separated metadata columns to exclude from the
                       gene set (default cell,OUcl,OUct; '' for none)
  --check-scaling      report marginal logLik under scaleHeight TRUE and FALSE
  --quiet              suppress progress messages
  -h, --help           this message
"

# ---------------------------------------------------------------------------
# argument parsing
# ---------------------------------------------------------------------------
parse_args <- function(argv, defaults) {
  opts <- c(defaults, list(check_scaling = FALSE, quiet = FALSE))
  as_logical <- function(x) {
    v <- toupper(x)
    if (v %in% c("TRUE", "T", "1", "YES")) TRUE
    else if (v %in% c("FALSE", "F", "0", "NO")) FALSE
    else stop("Expected TRUE/FALSE, got: ", x, call. = FALSE)
  }
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[[i]]
    val <- function() {
      if (i + 1L > length(argv)) stop("Missing value for ", a, call. = FALSE)
      argv[[i + 1L]]
    }
    if (a %in% c("-h", "--help")) { cat(USAGE); quit(save = "no", status = 0) }
    else if (a == "--check-scaling") { opts$check_scaling <- TRUE; i <- i + 1L }
    else if (a == "--quiet")         { opts$quiet <- TRUE;         i <- i + 1L }
    else if (a == "--all")           { opts$all <- TRUE;           i <- i + 1L }
    else if (a == "--internal-only") { opts$internal_only <- TRUE; i <- i + 1L }
    else {
      key <- switch(a,
        "--gene" = "gene", "--best-fit" = "best_fit", "--params-dir" = "params_dir",
        "--tid" = "tid", "--tree" = "tree", "--counts" = "counts",
        "--ref" = "ref", "--out" = "out", "--scout-src" = "scout_src",
        "--scale-tree" = "scale_tree", "--normalize" = "normalize",
        "--fixed-root" = "fixed_root", "--qc-column" = "qc_column",
        "--qc-value" = "qc_value", "--limit" = "limit", "--verify" = "verify",
        "--species-key" = "species_key", "--regimes" = "regimes",
        "--blacklist" = "blacklist",
        stop("Unknown argument: ", a, call. = FALSE))
      v <- val()
      opts[[key]] <-
        if (key %in% c("scale_tree", "normalize", "fixed_root")) as_logical(v)
        else if (key %in% c("limit", "verify")) as.integer(v)
        # Dataset wiring: comma-separated, so one copy of this script serves
        # every dataset instead of being forked per run (4_lg13 previously kept
        # a near-identical copy differing only in these three lines).
        else if (key %in% c("regimes", "blacklist"))
          if (!nzchar(v)) character(0) else trimws(strsplit(v, ",")[[1]])
        else v
      i <- i + 2L
    }
  }
  opts
}

opt <- parse_args(commandArgs(trailingOnly = TRUE), DEFAULTS)
say <- function(...) if (!opt$quiet) cat(...)

if (is.null(opt$gene) && !opt$all) {
  cat(USAGE); stop("give either --gene <GENE> or --all", call. = FALSE)
}
if (!is.null(opt$gene) && opt$all) {
  stop("--gene and --all are mutually exclusive.", call. = FALSE)
}

# ---------------------------------------------------------------------------
# reach SCOUT's internals (none of the ones needed here are exported)
# ---------------------------------------------------------------------------
if (!is.null(opt$scout_src)) {
  say("Sourcing SCOUT from ", opt$scout_src, "\n")
  suppressWarnings(suppressMessages({
    library(phylolm);   library(corpcor); library(nloptr)
    library(dplyr);     library(stringr); library(paleotree)
  }))
  for (f in c("SCOUT_EM_utils.R", "SCOUT_EM.R")) {
    p <- file.path(opt$scout_src, f)
    if (!file.exists(p)) stop("Not found: ", p, call. = FALSE)
    source(p)
  }
  scout <- function(nm) get(nm, envir = globalenv())
} else {
  suppressWarnings(suppressMessages(library(SCOUT)))
  scout <- function(nm) getFromNamespace(nm, "SCOUT")
}

formatSCOUT      <- scout("formatSCOUT")
preprocessTree   <- scout("preprocessTree")
preprocessGene   <- scout("preprocessGene")
compute_VCV      <- scout("compute_VCV")
compute_W_matrix <- scout("compute_W_matrix")
e_step           <- scout("e_step")

# ---------------------------------------------------------------------------
# 1. the best-fit table -> the gene list, and per gene model/regime/alpha/sigma/tau
# ---------------------------------------------------------------------------
say("Reading best-fit table: ", opt$best_fit, "\n")
bf <- read.csv(opt$best_fit, row.names = 1, stringsAsFactors = FALSE)

if (opt$all) {
  n0 <- nrow(bf)
  if (!is.null(opt$qc_column) && nzchar(opt$qc_column)) {
    if (opt$qc_column %in% colnames(bf)) {
      bf <- bf[!is.na(bf[[opt$qc_column]]) & bf[[opt$qc_column]] == opt$qc_value, , drop = FALSE]
      say(sprintf("QC filter %s == '%s': kept %d of %d rows\n",
                  opt$qc_column, opt$qc_value, nrow(bf), n0))
    } else {
      say(sprintf("NOTE: column '%s' absent from the best-fit table; no QC filter applied.\n",
                  opt$qc_column))
    }
  }
  if (anyDuplicated(bf$gene_name))
    stop("Duplicate gene_name rows in the best-fit table after filtering.", call. = FALSE)
  if (!is.null(opt$limit) && !is.na(opt$limit) && opt$limit < nrow(bf)) {
    bf <- bf[seq_len(opt$limit), , drop = FALSE]
    say(sprintf("--limit %d: restricted to the first %d genes\n", opt$limit, nrow(bf)))
  }
  if (!nrow(bf)) stop("No genes left after filtering.", call. = FALSE)
} else {
  bf <- bf[bf$gene_name == opt$gene, , drop = FALSE]
  if (nrow(bf) == 0L)
    stop(sprintf("Gene '%s' not found in %s.", opt$gene, opt$best_fit), call. = FALSE)
  if (nrow(bf) > 1L)
    stop(sprintf("Gene '%s' has %d rows in the best-fit table; expected exactly one.",
                 opt$gene, nrow(bf)), call. = FALSE)
}

tid <- if (is.null(opt$tid)) as.character(bf$dataset[1]) else opt$tid
if (length(unique(bf$dataset)) > 1L)
  stop("The best-fit table mixes datasets: ", paste(unique(bf$dataset), collapse = ", "),
       ". Filter it, or pass --tid.", call. = FALSE)

regimes_needed <- sort(unique(as.character(bf$regime)))
say(sprintf("\nGenes      : %d\nDataset    : %s\nRegimes    : %s\n",
            nrow(bf), tid,
            paste(sprintf("%s (%d)", regimes_needed,
                          as.integer(table(bf$regime)[regimes_needed])), collapse = ", ")))

# ---------------------------------------------------------------------------
# 2. theta tables, one per regime (read once, not once per gene)
# ---------------------------------------------------------------------------
theta_tables <- list()
for (rg in regimes_needed) {
  pfile <- file.path(opt$params_dir, sprintf("%s_all_genes_%s_parameters.csv", tid, rg))
  if (!file.exists(pfile)) stop("Parameter file not found: ", pfile, call. = FALSE)
  # SCOUT writes these with write.table(): whitespace separated, quoted, row names
  tb <- read.table(pfile, header = TRUE, stringsAsFactors = FALSE)
  rownames(tb) <- tb$gene_name
  theta_tables[[rg]] <- tb
}

get_gene_params <- function(row) {
  rg <- as.character(row$regime)
  tb <- theta_tables[[rg]]
  g  <- as.character(row$gene_name)
  if (!g %in% rownames(tb))
    stop(sprintf("no row for '%s' in the %s parameter table", g, rg), call. = FALSE)
  prow <- tb[g, , drop = FALSE]

  theta_cols <- grep("^theta_", colnames(prow), value = TRUE)
  th <- as.numeric(prow[1, theta_cols])
  names(th) <- sub("^theta_", "", theta_cols)

  p <- list(tau = as.numeric(row$tau), alpha = as.numeric(row$alpha),
            sigma = as.numeric(row$sigma), theta = th)

  # Guard the parameters before they reach the covariance. e_step() stops on a
  # negative tau; nothing downstream of here would, because only tau^2 is ever
  # used -- a negative tau would silently be squared into a different tip-fog
  # variance and produce plausible-looking output. Same for a non-positive
  # sigma/alpha or a non-finite theta.
  bad <- character(0)
  if (!is.finite(p$tau)   || p$tau   <  0) bad <- c(bad, sprintf("tau=%s",   format(p$tau)))
  if (!is.finite(p$sigma) || p$sigma <= 0) bad <- c(bad, sprintf("sigma=%s", format(p$sigma)))
  if (!is.finite(p$alpha) || p$alpha <= 0) bad <- c(bad, sprintf("alpha=%s", format(p$alpha)))
  if (any(!is.finite(p$theta)))            bad <- c(bad, "theta has non-finite entries")
  if (length(bad))
    stop("unusable fitted parameters: ", paste(bad, collapse = ", "), call. = FALSE)

  # alpha/sigma/tau must agree between the two tables -- it is the same fit
  for (nm in c("alpha", "sigma", "tau")) {
    a <- p[[nm]]; b <- as.numeric(prow[[nm]])
    if (!is.finite(b)) {
      warning(sprintf("[%s] %s is not finite in the parameter table", g, nm))
    } else if (abs(a - b) / max(abs(a), 1e-12) > 1e-6) {
      warning(sprintf("[%s] %s differs: best-fit %.10g vs parameter table %.10g", g, nm, a, b))
    }
  }

  p
}

# ---------------------------------------------------------------------------
# 3. rebuild the exact tree/data objects the fits ran on (once)
# ---------------------------------------------------------------------------
say("\nRebuilding SCOUT inputs (formatSCOUT)...\n")
counts <- read.csv(opt$counts, row.names = 1)
tree   <- ape::read.tree(opt$tree)

idata <- formatSCOUT(
  tree_path     = tree,
  metadata_path = counts,
  species_key   = opt$species_key,
  anc_infer     = opt$infer_anc,
  outpath       = tempdir(),
  regimes       = opt$regimes,
  normalize     = opt$normalize,
  smoothing_k   = NULL,
  blacklist     = opt$blacklist,
  logfile       = NULL
)

missing_genes <- setdiff(bf$gene_name, idata$gene_cols)
if (length(missing_genes))
  stop(sprintf("%d gene(s) in the best-fit table are absent from the counts data, e.g. %s",
               length(missing_genes), paste(head(missing_genes, 5), collapse = ", ")),
       call. = FALSE)

# preprocessGene() pulls X as meta_data[, gene] named by species; do the same
# directly so a batch run does not pay for its pic()/tapply() work per gene.
species_col <- idata$meta_data[["species"]]
gene_X <- function(gene) { x <- idata$meta_data[[gene]]; names(x) <- species_col; x }

# ---------------------------------------------------------------------------
# 4. per-regime context: everything that does not depend on the gene
# ---------------------------------------------------------------------------
# Regime segments along root->node, generalised from SCOUT's build_regime_paths()
# (which hard-codes the tips). Identical segment logic, plus a guard for the
# root, whose path is a single node.
regime_segments_all <- function(phy, node_regimes, root_dists) {
  n_tips <- length(phy$tip.label)
  n_all  <- n_tips + phy$Nnode
  root   <- n_tips + 1L

  parent_of <- integer(n_all)
  parent_of[phy$edge[, 2]] <- phy$edge[, 1]

  path_to_root <- function(node) {
    path <- node
    cur  <- node
    while (cur != root) { cur <- parent_of[cur]; path <- c(cur, path) }
    path
  }

  lapply(seq_len(n_all), function(i) {
    path_nodes <- path_to_root(i)
    segments <- list()
    current_regime     <- node_regimes[path_nodes[1]]
    segment_start_dist <- root_dists[path_nodes[1]]

    if (length(path_nodes) >= 2L) {
      for (k in 2:length(path_nodes)) {
        node <- path_nodes[k]
        if (node_regimes[node] != current_regime) {
          segments[[length(segments) + 1L]] <- list(
            regime     = current_regime,
            start_dist = segment_start_dist,
            end_dist   = root_dists[path_nodes[k - 1L]]
          )
          current_regime     <- node_regimes[node]
          segment_start_dist <- root_dists[path_nodes[k - 1L]]
        }
      }
    }
    segments[[length(segments) + 1L]] <- list(
      regime     = current_regime,
      start_dist = segment_start_dist,
      end_dist   = root_dists[path_nodes[length(path_nodes)]]
    )
    segments
  })
}

seg_key <- function(s) paste(vapply(s, function(z)
  sprintf("%s:%.12g:%.12g", as.character(z$regime), z$start_dist, z$end_dist),
  character(1)), collapse = "|")

lcp <- function(x) {
  if (!length(x)) return(NA_character_)
  p <- x[1]
  for (s in x[-1]) {
    n <- min(nchar(p), nchar(s)); k <- 0L
    while (k < n && substr(p, k + 1L, k + 1L) == substr(s, k + 1L, k + 1L)) k <- k + 1L
    p <- substr(p, 1L, k)
    if (!nchar(p)) break
  }
  if (nchar(p)) p else NA_character_
}

build_regime_context <- function(rg, example_gene) {
  tree_info  <- preprocessTree(inputs = idata, reg = rg,
                               scaleHeight = opt$scale_tree,
                               root.fixed  = opt$fixed_root,
                               root.age    = NULL,
                               skipTau     = FALSE)
  defaults   <- tree_info$defaults
  tree_const <- tree_info$tree_const
  phy        <- tree_const$tree
  edges      <- tree_const$edges

  # runSCOUT() rewrites unique_regimes from the initial theta names before
  # fitting; replicate that or W's columns will not line up with theta. Those
  # names are levels(phy$states), so they are the same for every gene in this
  # regime -- taking them from one example gene is enough.
  gi <- preprocessGene(idata, example_gene, tree_const,
                       defaults$root.state, defaults$add.root,
                       defaults$model, defaults$skipTau)
  if (any(names(gi$X) != phy$tip.label))
    stop("Tip labels and data names are out of sync.", call. = FALSE)
  stopifnot(identical(unname(gi$X), unname(gene_X(example_gene))))   # direct X matches

  if (any(names(gi$param_init$theta) != defaults$parsed_alt_tree$unique_regimes))
    defaults$parsed_alt_tree$unique_regimes <- names(gi$param_init$theta)
  if (!"root" %in% defaults$parsed_alt_tree$unique_regimes && defaults$add.root)
    defaults$parsed_alt_tree$unique_regimes <- c("root", defaults$parsed_alt_tree$unique_regimes)
  unique_regimes <- defaults$parsed_alt_tree$unique_regimes

  # run_em() stamps these onto defaults before any E-step; e_step() reads
  # defaults$verbose directly and errors on NULL.
  defaults$verbose <- FALSE

  n_tips  <- length(phy$tip.label)
  n_all   <- n_tips + phy$Nnode
  tip_idx <- seq_len(n_tips)

  root_dists <- ape::node.depth.edgelength(phy)          # t_i for every node

  # preprocessTree builds node_regimes as c(tip.states, node.label); the phylo
  # it returns carries both, so this reproduces the same vector exactly
  # (including the factor-vs-numeric coercion, which differs between OUM and
  # BM1/OU1).
  node_regimes <- c(phy$states, phy$node.label)
  names(node_regimes) <- seq_len(n_all)
  stopifnot(length(node_regimes) == n_all)

  # s_ij = root->MRCA(i,j) distance for all node pairs. dist.nodes() gives the
  # patristic distance d_ij among all nodes, so s_ij = (t_i + t_j - d_ij)/2.
  d_all      <- ape::dist.nodes(phy)
  shared_all <- (outer(root_dists, root_dists, "+") - d_all) / 2
  dimnames(shared_all) <- NULL

  alt_all <- list(
    n_leaves       = n_all,
    regime_paths   = regime_segments_all(phy, node_regimes, root_dists),
    unique_regimes = unique_regimes,
    leaf_dists     = root_dists,
    shared_lengths = shared_all,
    root.state     = defaults$parsed_alt_tree$root.state
  )

  # --- the tip block of the all-node object must reproduce SCOUT's own object -
  pa <- defaults$parsed_alt_tree
  stopifnot(isTRUE(all.equal(unname(alt_all$leaf_dists[tip_idx]),
                             unname(pa$leaf_dists), tolerance = 1e-10)))
  stopifnot(isTRUE(all.equal(unname(alt_all$shared_lengths[tip_idx, tip_idx]),
                             unname(pa$shared_lengths), tolerance = 1e-10)))
  stopifnot(identical(vapply(alt_all$regime_paths[tip_idx], seg_key, character(1)),
                      vapply(pa$regime_paths,               seg_key, character(1))))

  # --- node metadata (tree-level, so computed once per regime) ---------------
  desc <- vector("list", n_all)
  for (i in tip_idx) desc[[i]] <- i
  po <- ape::reorder.phylo(phy, "postorder")     # children always before parents
  for (e in seq_len(nrow(po$edge))) {
    p <- po$edge[e, 1]; ch <- po$edge[e, 2]
    desc[[p]] <- c(desc[[p]], desc[[ch]])
  }

  node_label <- character(n_all); n_desc <- integer(n_all)
  for (j in seq_len(n_all)) {
    n_desc[j] <- length(desc[[j]])
    # Internal nodes carry no label in this newick. On a C. elegans lineage tree
    # the longest common prefix of a clade's tip names *is* the ancestral cell
    # name (ABalaaaalp, ABalaaaarr, ... -> ABalaaaa). Clades spanning several
    # founder lineages share no prefix (the root above all); fall back to the
    # ape node number so every row stays identifiable.
    node_label[j] <- if (j <= n_tips) phy$tip.label[j] else lcp(phy$tip.label[desc[[j]]])
    if (is.na(node_label[j])) node_label[j] <- paste0("node", j)
  }

  parent <- integer(n_all)
  parent[phy$edge[, 2]] <- phy$edge[, 1]
  parent[n_tips + 1L]   <- NA_integer_
  brlen <- rep(NA_real_, n_all)
  brlen[phy$edge[, 2]]  <- phy$edge.length

  list(regime = rg, defaults = defaults, phy = phy, edges = edges,
       alt_all = alt_all, unique_regimes = unique_regimes,
       n_tips = n_tips, n_all = n_all, tip_idx = tip_idx,
       node_regimes = node_regimes, root_dists = root_dists,
       node_label = node_label, n_desc = n_desc, parent = parent, brlen = brlen)
}

# ---------------------------------------------------------------------------
# 5. E[Z] and Var[Z] at every node, for one gene
# ---------------------------------------------------------------------------
# e_step()/m_step duplicate theta in the BM1 fixed-root case; same rule here.
expand_theta <- function(paras, defaults) {
  if (isTRUE(defaults$add.root) && defaults$model == "BM1") c(paras$theta, paras$theta)
  else paras$theta
}

# extract_parameters() named the theta columns paste0('theta_', unique_regimes),
# dropping the leading 'root' entry for BM1. Undo that mapping.
order_theta <- function(theta_csv, ctx) {
  expected <- if (ctx$defaults$model == "BM1" && isTRUE(ctx$defaults$add.root))
    ctx$unique_regimes[-1] else ctx$unique_regimes
  if (!setequal(names(theta_csv), expected))
    stop(sprintf("theta columns (%s) do not match the model's regimes (%s)",
                 paste(names(theta_csv), collapse = ","),
                 paste(expected, collapse = ",")), call. = FALSE)
  theta_csv[expected]
}

node_posterior <- function(ctx, paras, X) {
  defaults <- ctx$defaults; tip_idx <- ctx$tip_idx
  V_all <- compute_VCV(ctx$alt_all, paras$alpha, paras$sigma, add.root = defaults$add.root)
  W_all <- compute_W_matrix(ctx$alt_all, paras$alpha, add.root = defaults$add.root)
  th    <- expand_theta(paras, defaults)
  if (ncol(W_all) != length(th))
    stop(sprintf("W has %d columns but theta has %d entries.",
                 ncol(W_all), length(th)), call. = FALSE)

  mu_all <- as.vector(W_all %*% th)

  V_TT <- V_all[tip_idx, tip_idx, drop = FALSE]
  V_AT <- V_all[, tip_idx, drop = FALSE]

  # same regularisation as e_step()
  S <- V_TT + diag(paras$tau^2, length(tip_idx)) + 1e-8 * diag(length(tip_idx))
  L <- tryCatch(chol(S), error = function(e) NULL)
  S_inv <- if (is.null(L)) corpcor::pseudoinverse(S) else chol2inv(L)

  resid <- as.vector(X) - mu_all[tip_idx]
  K_A   <- V_AT %*% S_inv                      # Kalman gain, all nodes x tips

  list(E_Z        = mu_all + as.vector(K_A %*% resid),
       Var_Z      = V_all - K_A %*% t(V_AT),
       prior_mean = mu_all)
}

gene_rows <- function(ctx, row, paras, X, post) {
  n_all <- ctx$n_all; n_tips <- ctx$n_tips
  sd_Z  <- sqrt(pmax(diag(post$Var_Z), 0))
  X_obs <- rep(NA_real_, n_all); X_obs[ctx$tip_idx] <- as.vector(X)

  df <- data.frame(
    gene          = as.character(row$gene_name),
    dataset       = tid,
    model         = as.character(row$model),
    regime        = ctx$regime,
    node          = seq_len(n_all),
    node_type     = ifelse(seq_len(n_all) <= n_tips, "tip",
                    ifelse(seq_len(n_all) == n_tips + 1L, "root", "internal")),
    node_label    = ctx$node_label,
    node_regime   = as.character(ctx$node_regimes),
    parent        = ctx$parent,
    branch_length = ctx$brlen,
    depth         = ctx$root_dists,
    n_desc_tips   = ctx$n_desc,
    X_obs         = X_obs,
    prior_mean    = post$prior_mean,
    E_Z           = post$E_Z,
    sd_Z          = sd_Z,
    lo95          = post$E_Z - 1.959964 * sd_Z,
    hi95          = post$E_Z + 1.959964 * sd_Z,
    stringsAsFactors = FALSE
  )
  if (opt$internal_only) df <- df[df$node_type != "tip", , drop = FALSE]
  df
}

# ---------------------------------------------------------------------------
# 6. main loop -- one file, streamed
# ---------------------------------------------------------------------------
outfile <- {
  if (!is.null(opt$out)) opt$out
  else if (opt$all) sprintf("%s_best_fit_node_EZ.csv.gz", tid)
  else sprintf("%s_%s_node_EZ.csv", bf$gene_name[1], bf$regime[1])
}

con <- if (grepl("\\.gz$", outfile)) gzfile(outfile, "wt") else file(outfile, "wt")
on.exit(try(close(con), silent = TRUE), add = TRUE)

bf <- bf[order(bf$regime, bf$gene_name), , drop = FALSE]
contexts <- list()
wrote_header <- FALSE
n_rows <- 0L
failures <- list()
checked <- setNames(integer(length(regimes_needed)), regimes_needed)
check_report <- list()
t0 <- Sys.time()

for (k in seq_len(nrow(bf))) {
  row <- bf[k, , drop = FALSE]
  g   <- as.character(row$gene_name)
  rg  <- as.character(row$regime)

  if (is.null(contexts[[rg]])) {
    say(sprintf("\n[%s] building tree context...\n", rg))
    contexts[[rg]] <- build_regime_context(rg, g)
    ctx <- contexts[[rg]]
    say(sprintf("[%s] model=%s add.root=%s assume.station=%s scaleHeight=%s | regimes: %s | root.state: %s | height %.6g\n",
                rg, ctx$defaults$model, ctx$defaults$add.root, ctx$defaults$assume.station,
                ctx$defaults$scaleHeight, paste(ctx$unique_regimes, collapse = ","),
                as.character(ctx$alt_all$root.state), max(ctx$root_dists)))
    say(sprintf("[%s] tip block reproduces SCOUT's parsed_alt_tree\n", rg))
  }
  ctx <- contexts[[rg]]

  res <- tryCatch({
    paras <- get_gene_params(row)
    paras$theta <- order_theta(paras$theta, ctx)
    if (ctx$defaults$model == "BM1") paras$alpha <- 1e-10   # run_em() pins this
    X <- gene_X(g)
    post <- node_posterior(ctx, paras, X)

    # cross-check the tip block against SCOUT's own e_step()
    do_check <- opt$verify < 0L || checked[[rg]] < opt$verify
    if (do_check) {
      e_ref <- e_step(X, ctx$phy, ctx$edges, paras, ctx$defaults, diagnose = FALSE)
      if (!is.null(e_ref)) {
        dmax <- max(abs(post$E_Z[ctx$tip_idx] - as.vector(e_ref$E_Z)))
        rel  <- dmax / max(abs(as.vector(e_ref$E_Z)), 1e-12)
        checked[[rg]] <- checked[[rg]] + 1L
        check_report[[length(check_report) + 1L]] <-
          data.frame(gene = g, regime = rg, max_abs = dmax, max_rel = rel)
        if (rel > 1e-8)
          warning(sprintf("[%s] tip E[Z] differs from e_step() by rel %.3g", g, rel))
      }
    }
    gene_rows(ctx, row, paras, X, post)
  }, error = function(e) {
    failures[[length(failures) + 1L]] <<- data.frame(gene = g, regime = rg,
                                                     message = conditionMessage(e))
    NULL
  })

  if (!is.null(res)) {
    write.table(res, con, sep = ",", row.names = FALSE,
                col.names = !wrote_header, qmethod = "double")
    wrote_header <- TRUE
    n_rows <- n_rows + nrow(res)
  }

  if (opt$all && (k %% 100L == 0L || k == nrow(bf))) {
    el <- as.numeric(Sys.time() - t0, units = "secs")
    say(sprintf("  %d/%d genes | %d rows | %.1f s elapsed | %.2f s/gene | eta %.1f min\n",
                k, nrow(bf), n_rows, el, el / k, (el / k) * (nrow(bf) - k) / 60))
  }
}

close(con)
on.exit()

# ---------------------------------------------------------------------------
# 7. optional cross-check against the stored results object (single gene)
# ---------------------------------------------------------------------------
if (!is.null(opt$ref)) {
  if (opt$all) {
    say("\nNOTE: --ref is a single-gene cross-check; skipped in --all mode.\n")
  } else {
    g <- as.character(bf$gene_name[1]); rg <- as.character(bf$regime[1])
    ctx <- contexts[[rg]]
    say("\nLoading reference object ", opt$ref, "\n")
    say("  (the full results .rds needs ~32 GB RAM and several minutes)\n")
    ref <- readRDS(opt$ref)
    hit <- if (!is.null(ref$paras) && !is.null(ref$defaults)) ref else {
      h <- NULL
      for (e in ref) {
        s <- as.character(unlist(e$settings))
        if (length(s) >= 2L && s[1] == g && s[2] == rg) { h <- e; break }
      }
      h
    }
    rm(ref); invisible(gc())

    if (is.null(hit)) {
      warning(sprintf("No element for %s / %s in %s", g, rg, opt$ref))
    } else {
      paras <- get_gene_params(bf[1, , drop = FALSE])
      theta <- order_theta(paras$theta, ctx)
      spa <- hit$defaults$parsed_alt_tree
      st  <- as.character(unlist(hit$settings))
      say(sprintf("  element     : %s / %s\n", st[1], st[2]))
      say(sprintf("  scaleHeight : stored %s | here %s\n",
                  hit$defaults$scaleHeight, ctx$defaults$scaleHeight))
      say(sprintf("  |theta_csv - theta_rds|  max = %.3g\n",
                  max(abs(unlist(hit$paras$theta) - as.numeric(theta)))))
      say(sprintf("  |leaf_dists - stored|    max = %.3g   (settles scaleHeight)\n",
                  max(abs(ctx$alt_all$leaf_dists[ctx$tip_idx] - spa$leaf_dists))))
      say(sprintf("  |shared_lengths - stored| max = %.3g\n",
                  max(abs(ctx$alt_all$shared_lengths[ctx$tip_idx, ctx$tip_idx] -
                          spa$shared_lengths))))
      # The regimes painted on internal nodes come from ape::ace() ancestral-state
      # reconstruction, re-run here. Confirm it reproduced the original painting,
      # otherwise internal-node W rows attribute to the wrong theta.
      say(sprintf("  root.state  : stored %s | here %s\n",
                  as.character(spa$root.state), as.character(ctx$alt_all$root.state)))
      say(sprintf("  tip regime_paths identical : %s   (ancestral states reproduced)\n",
                  identical(vapply(ctx$alt_all$regime_paths[ctx$tip_idx], seg_key, character(1)),
                            vapply(spa$regime_paths, seg_key, character(1)))))
    }
  }
}

# ---------------------------------------------------------------------------
# 8. optional: which scaleHeight setting do these parameters belong to?
# ---------------------------------------------------------------------------
if (opt$check_scaling) {
  g   <- as.character(bf$gene_name[1]); rg <- as.character(bf$regime[1])
  ctx <- contexts[[rg]]
  paras <- get_gene_params(bf[1, , drop = FALSE])
  paras$theta <- order_theta(paras$theta, ctx)
  if (ctx$defaults$model == "BM1") paras$alpha <- 1e-10
  X <- gene_X(g)
  say(sprintf("\n--check-scaling (%s): marginal logLik of X under each setting\n", g))
  marg_ll <- function(scale_tree) {
    ti <- preprocessTree(inputs = idata, reg = rg, scaleHeight = scale_tree,
                         root.fixed = opt$fixed_root, root.age = NULL, skipTau = FALSE)
    dd <- ti$defaults
    dd$parsed_alt_tree$unique_regimes <- ctx$unique_regimes
    V  <- compute_VCV(dd$parsed_alt_tree, paras$alpha, paras$sigma, add.root = dd$add.root)
    W  <- compute_W_matrix(dd$parsed_alt_tree, paras$alpha, add.root = dd$add.root)
    mu <- as.vector(W %*% expand_theta(paras, dd))
    S  <- V + diag(paras$tau^2, nrow(V)) + 1e-8 * diag(nrow(V))
    ch <- chol(S)
    r  <- as.vector(X) - mu
    q  <- backsolve(ch, r, transpose = TRUE)
    -0.5 * (length(r) * log(2 * pi) + 2 * sum(log(diag(ch))) + sum(q^2))
  }
  for (s in c(TRUE, FALSE)) {
    ll <- tryCatch(marg_ll(s), error = function(e) NA_real_)
    say(sprintf("  scaleHeight = %-5s  logLik(X) = %.4f\n", s, ll))
  }
  say("  The setting the parameters were fitted under should score far higher.\n")
}

# ---------------------------------------------------------------------------
# 9. summary
# ---------------------------------------------------------------------------
n_ok <- nrow(bf) - length(failures)
say(sprintf("\nWrote %d rows for %d gene(s) to %s\n",
            n_rows, n_ok, normalizePath(outfile, mustWork = FALSE)))
say(sprintf("Nodes per gene: %d (%s)\n",
            if (n_ok > 0L) n_rows / n_ok else 0L,
            if (opt$internal_only) "internal + root only" else "all tips + internal"))

if (length(check_report)) {
  cr <- do.call(rbind, check_report)
  say(sprintf("Cross-check vs e_step() on %d gene(s): max |dE[Z]| = %.3g, max rel = %.3g\n",
              nrow(cr), max(cr$max_abs), max(cr$max_rel)))
}

if (length(failures)) {
  fr <- do.call(rbind, failures)
  say(sprintf("\n%d gene(s) FAILED and were skipped:\n", nrow(fr)))
  for (i in seq_len(min(nrow(fr), 20L)))
    say(sprintf("  %-20s [%s] %s\n", fr$gene[i], fr$regime[i], fr$message[i]))
  if (nrow(fr) > 20L) say(sprintf("  ... and %d more\n", nrow(fr) - 20L))
  ffile <- paste0(sub("\\.gz$", "", sub("\\.csv$", "", outfile)), "_failures.csv")
  write.csv(fr, ffile, row.names = FALSE)
  say("  failure list: ", ffile, "\n")
}
