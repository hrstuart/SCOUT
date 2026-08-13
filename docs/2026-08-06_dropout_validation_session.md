# Session report — held-out dropout validation for SCOUT-EM

**Date:** 2026-08-06
**Branch:** `scout_em_redo`
**Starting commit:** `a1c9213` (updating simulation strategy to include outside generated tree)

---

## 1. Goal

Build a script that takes a newick tree, randomly drops a clade, runs SCOUT-EM on data
simulated *before* the dropout, and then runs the E-step to compare the recovered latent
expression against the true simulated expression — to demonstrate that the latent process can
recover masked values, i.e. that the model is appropriate for these data.

---

## 2. Design decisions

### 2.1 The E-step must be block-conditional, not full-X

**Decision:** implement a block-conditional estimator rather than calling `e_step` on the full
tree with the complete `X`.

**Why:** `e_step` (`R/SCOUT_EM.R:73-74`) computes `E_Z = Wθ + K(X − Wθ)` over *all* tips with
`K = V(V + τ²I)⁻¹`. Handing it the full `X` puts the masked tips' own observations inside the
conditioning set, and since τ is bounded at ≥ 0.05 while V is typically much larger, `K ≈ I` and
`E_Z` at the masked tips essentially returns their own `X`. The plot would look perfect
regardless of whether the model fits. The block form conditions only on retained tips:

```
E[Z_h | X_o]   = W_hθ + V_ho (V_oo + τ²I)⁻¹ (X_o − W_oθ)
Var[Z_h | X_o] = V_hh   − V_ho (V_oo + τ²I)⁻¹ V_oh
```

`e_step` is the special case `h = o = all tips`, which became the primary correctness test.
The variance term is free and gives calibrated prediction intervals, so coverage is reportable
alongside RMSE.

### 2.2 Two masking arms, because clade dropout is rank-one

**Decision:** every replicate runs a `clade` arm and a size-matched `random` arm.

**Why:** for masked tip `i` and retained tip `o`, `MRCA(i,o)` is the same node for all `i` in
the clade, so

```
V[i,o] = (σ²/2α)·exp(−α(t_i + t_o − 2s(o)))·(1 − exp(−2α·s(o))) = exp(−α·t_i) · g(o)
```

— **`V_ho` is rank one** (holds for `add.root` TRUE and FALSE, `compute_VCV:361,363`). The
correction collapses to one scalar: `E[Z_h|X_o]_i = (W_hθ)_i + c·exp(−α·t_i)`. On an ultrametric
tree (which `generate_tree_hbd_reverse` produces) predictions within the clade are *constant* up
to the regime-weight term. Per-tip RMSE inside a masked clade is therefore floored by
irreducible within-clade variance, and a flat scatter is correct model behaviour, not evidence
against the model.

Consequence for the writeup: score the clade arm on **clade-mean error and coverage**, and use
the random arm (full-rank `V_ho`) for **per-tip RMSE and correlation**. Verified numerically —
`qr(V[h,o])$rank == 1` for clade, `== k` for random.

### 2.3 Ground truth is the EVF matrix, not the counts matrix

**Decision:** score against `sim$simulated_OU$evfs`, fit with `normalize = FALSE`.

**Why:** the EVF matrix (`R/simulate.R:163`) is the raw `OUwie.sim.edited` output — `Z` exactly,
no observation noise — so RMSE against `E[Z_h|X_o]` is well-posed and on the same scale. The
counts matrix runs latent factors through `GeneEffects` mixing and beta-Poisson sampling
(`R/simulate.R:130-137`), and only the first of three kinetic slices (`eevfs[[i]][[1]]`) is
returned as "truth", so there is no 1:1 gene↔latent map and only rank correlation would be
interpretable.

### 2.4 Baselines are mandatory

Absolute RMSE is uninterpretable. Four predictors are reported per gene per arm:

| Predictor | Definition | Role |
|---|---|---|
| `masked` | `E[Z_h\|X_o]`, masked-fit params | the estimand |
| `prior` | `W_hθ` only | **the null to beat** — no phylogenetic borrowing |
| `oracle` | same estimator, full-tree params | isolates "masking hurt the params" from "this is hard" |
| `global` | `mean(X_o)` | crude floor |

### 2.5 Other decisions

- **`scaleHeight` hardcoded FALSE.** `preprocessTree:789` divides edge lengths by `Tmax`; if the
  masked and full trees differ in height, α and σ² are on different time scales and cannot
  transfer. Note this diverges from `260804_sim_test_2.ipynb`, which uses `scale_tree = TRUE`.
- **Root's children excluded from clade sampling.** Dropping one collapses the root and shifts
  every root-to-tip distance.
- **Regime-eliminating masks rejected.** An unrepresented regime has no estimable θ.
- **AIC computed inline**, not via `annotate_history` — the experiment varies tip count per
  replicate by construction. (You independently fixed the hardcoded `ntips = 256` during the
  session; the new code is unaffected either way since it never calls `annotate_history`.)
- **Entry point takes `sim_res = NULL`.** Supplying a simulation reuses it across replicates so
  variation reflects clade choice, not simulation noise.

---

## 3. Code changes

### New files

| File | Contents |
|---|---|
| `R/dropout_validation.R` | `runSCOUT.dropout` (exported) + helpers: `predict_masked_latent`, `sample_clade`, `sample_random_tips`, `descendant_tips`, `align_regimes`, `check_regime_coverage`, `fit_paras`, `fit_aic`, `score_recovery`, `prep_dropout_tree`, `index_fits` |
| `tests/testthat.R`, `tests/testthat/test-dropout-validation.R` | 89 assertions — the package previously had no tests |
| `examples/CladeDropoutVignette.ipynb` | worked example |
| `man/*.Rd` (12 files) | roxygen-generated |

### Modified files

| File | Change | Why |
|---|---|---|
| `R/simulate.R:323` | `reorder.phylo(phy, "cladewise")` in `OUwie.sim.edited` | **Bug fix — see §4.1** |
| `R/SCOUT_EM_utils.R:1165,1272` | `plan(...)` → `future::plan(future::multisession, ...)` | **Bug fix — see §4.2** |
| `R/SCOUT_EM_utils.R:322-330` | error on unmatched regime in `compute_W_matrix` | `which(regimes == seg$regime)` returning `integer(0)` made `W[i, integer(0)] <- numeric(0)` a **silent no-op**, yielding a wrong W with no error |
| `R/SCOUT_EM_utils.R:354-372` | vectorised `compute_VCV` | bit-identical output, 8× faster; it runs every EM iteration × gene × regime, and this experiment multiplies that by arms × replicates |
| `R/SCOUT_EM_utils.R:831` | vectorised `shared_lengths` | same, `mrca()` already computed |
| `DESCRIPTION` | added `future`, `future.apply`, `progressr`, `stringr` to Imports; `testthat` to Suggests | those four were `@import`ed in NAMESPACE but absent from DESCRIPTION |
| `NAMESPACE` | `export(runSCOUT.dropout)` | roxygen |

### Return value of `runSCOUT.dropout`

```r
list(clades       = replicate, arm, node, n_masked, masked_tips, status
     predictions  = replicate, arm, gene, tip, z_true, z_masked, z_prior, z_oracle, pred_sd
     metrics      = replicate, arm, gene, regime, predictor, rmse, mae, bias,
                    pearson, spearman, mean_error, coverage, converge
     params       = masked-fit and full-fit alpha/sigma/tau side by side, ll, AIC, converge
     model_select = replicate, arm, gene, true_model, best_full, best_masked
     settings     = call arguments + seed)
```

All five tables also written to `results_dir` as `<runid>_<date>_<table>.csv` plus one `.rds`.

---

## 4. Bugs found

### 4.1 The hbd_reverse tree silently reduced the simulated traits to i.i.d. noise

**The precondition.** `OUwie.sim.edited` integrates the OU process by walking the edge matrix in
row order, reading the ancestor's already-simulated value:

```r
x <- matrix(0, n, 1)                                   # R/simulate.R:399
x[ROOT,] <- theta0                                     # R/simulate.R:402
for (i in 1:length(edges[,1])) {                       # R/simulate.R:404
    x[edges[i,3],] <- x[edges[i,2],]*exp(-alpha*dt) + theta*(1-exp(-alpha*dt)) + noise
}                                                      # R/simulate.R:419
```

`x` is initialised to zeros, so this is correct **only if every parent edge appears before its
children**. When a parent has not yet been visited, `x[edges[i,2],]` reads `0` instead of the
parent's trait, and the child's trajectory restarts from zero.

**Why `generate_tree_hbd_reverse` specifically.** It constructs the tree backwards in time and
returns the edge matrix in that construction order — descending node number, with the root
processed *last*:

```
row:   1     2     3     4     5     6     7     8
anc:  199   199   198   198   197   197   196   196      root is node 101 -> visited LAST
```

**196 of 198 edges (99%) have an unvisited parent.** This is not a mild perturbation; it is
close to the worst possible ordering. `ape::read.tree` returns cladewise edges (0 bad edges),
which is why every file-based workflow in the package was unaffected and this went unnoticed.

**What the data actually was.** With nearly every parent read as zero, each tip's value is just
its own terminal branch's increment. 600 BM replicates, θ₀ = 10, σ² = 1, tree height 3.80:

| | before fix | after fix | theory |
|---|---|---|---|
| mean tip value | **−0.004** | 10.015 | 10 |
| mean tip variance | **0.459** | 3.950 | 3.803 (σ²·height) |
| mean terminal branch length | — | — | **0.458** |
| mean off-diagonal covariance | **−0.000** | 0.829 | 0.752 |
| cor(empirical *V_ij*, theoretical *V_ij*) | **+0.017** | +0.987 | 1 |

Tip variance matched the mean *terminal branch length* to three decimals, and off-diagonal
covariance was zero. **The simulated traits were i.i.d. Gaussian noise, one independent draw per
tip.** θ₀ never propagated and the tree contributed nothing beyond a per-tip variance scale.

**Where it bit in the notebook.** `260804_sim_test_2.ipynb` passes the *in-memory* castor tree
to the simulator, but writes it to newick and reads it back for fitting:

```r
simtree <- simtree.obj$trees[[1]]                             # cell 5  castor order, NOT preorder
simtest <- simulate_test_data(30, ncells, simtree, ...)       # cell 7  <- BROKEN
ape::write.tree(simtree, '...260805_..._fixed.nwk')           # cell 8  round-trip reorders it
scout.res <- SCOUT(tree.file = '...260805_..._fixed.nwk', ...)# cell 9  <- correctly preordered
```

So the *fitting* side received a correctly ordered tree and computed the right phylogenetic
covariance for that topology, while the *data* had none. SCOUT was being asked to discriminate
BM1 / OU1 / OUM on white noise.

**Consequence for model selection.** Same tree, same seed, same fitting code; only the simulator
differs (EVF matrix, 30 genes, n = 100, α = 3, σ² = 1, t₀ = 8):

| | broken simulation | fixed simulation |
|---|---|---|
| overall accuracy | **0.63** | **0.90** |
| BM1 recovered | **2 / 10** (8 called OU1) | 10 / 10 |
| OU1 recovered | 7 / 10 | 7 / 10 |
| OUM recovered | 10 / 10 | 10 / 10 |
| median fitted τ | 0.60 (all classes) | 0.30 |

The bias is systematic and predictable: BM covariance is σ²·*s_ij* and cannot be diagonal at any
parameter value, whereas OU approaches a diagonal as α grows, and the tip-fog term τ²**I** is
exactly diagonal. On white noise, OU can fit and BM cannot — so **BM1 genes are misclassified as
OU1**, and the inflated τ shows the model absorbing the noise. OUM still scored 10/10 because a
residual regime-mean offset survives — each tip retains θ_r·(1 − e^{−α·dt}) ≈ 11% of its regime
optimum — so OUM was being detected through a mean shift, not through phylogeny.

**Affected outputs.** Both prefixes in that notebook — `260804_simtree_100cells_hbd_rev_*` and
`260805_simtree_100cells_hbd_rev_fixed_*` — were produced this way and should be regenerated.
Anything built from a newick file on disk is unaffected.

**The fix.** `reorder.phylo(phy, "cladewise")` at the top of `OUwie.sim.edited`
(`R/simulate.R:323`), enforcing the precondition where it is required rather than assuming it of
the caller. Covered by a regression test.

### 4.2 `plan()` was never imported (pre-existing)

`runSCOUT` and `runSCOUT.batches` call `plan(multisession, ...)` from **future**, which is
imported nowhere — `future.apply` does not re-export it (`"plan" %in%
getNamespaceExports("future.apply")` is `FALSE`). It only ever worked because interactive
sessions happened to have `library(future)` loaded. Calling `runSCOUT` from inside the package
failed immediately with `could not find function "plan"`.

### 4.3 Ordering bug in my own estimator (found by testing, now fixed)

`predict_masked_latent` originally sorted `held_idx`, while the caller compared against
`X[held_idx]` in the *original* order. Clade masks come out near-sorted so the damage was small;
random masks were fully misaligned. This produced the false result that the `random` arm was
*worse* than baseline while the `clade` arm looked fine — an inverted signal that would have
been easy to misread as "SCOUT can't do per-tip prediction."

Caught by a known-parameter test: drawing `Z ~ N(μ, V)` directly from SCOUT's own model and
predicting with the *true* parameters still gave masked 1.964 vs global 1.667 and coverage 0.65
— impossible for a correct conditional Gaussian, which is Bayes-optimal by construction. Fixed
by preserving caller order, plus a name-alignment assertion at the call site and a dedicated
regression test.

**Lesson worth keeping:** the oracle test (known parameters, data from the model's own
generative process) is the one that isolates estimator bugs from fitting behaviour. It should be
the first thing run on any future change to the prediction code.

---

## 4b. Is `generate_tree_hbd_reverse` an appropriate simulation design?

Separating the software failure from the modelling question, because they are independent.

**The tree model itself is fine.** `generate_tree_hbd_reverse(n, rho = 1, lambda = 1, mu = 0)`
is a Yule process — pure birth, no extinction, complete sampling — which is the canonical model
for a dividing cell population and produces a valid ultrametric tree. Nothing about §4.1 was a
defect of the tree; it was a violated software precondition, now enforced. Keep using it.

**But α = 3 is inappropriate against these trees.** Tree height grows with *n*:

| n | height | mean branch | α·H at α = 3 | half-lives across the tree |
|---|---|---|---|---|
| 50 | 4.76 | 0.577 | 14.3 | 20.6 |
| 100 | 4.99 | 0.516 | 15.0 | 21.6 |
| 256 | 6.21 | 0.507 | 18.6 | 26.9 |
| 500 | 6.41 | 0.510 | 19.2 | 27.8 |
| 1000 | 7.55 | 0.502 | 22.7 | 32.7 |

The phylogenetic half-life is log(2)/α. At α = 3 that is 0.23 time units — about **6% of the
tree height**, so an OU trait forgets its ancestry roughly 22 half-lives before reaching the
tips. The OU1 and OUM arms of the benchmark are therefore near-white-noise by construction, and
the held-out validation measures this directly: for OU1 genes, masked RMSE 0.385 against a
constant-mean baseline of 0.392 — there is essentially nothing to recover. At α = 0.3 the same
genes give 0.812 against 1.121, a real effect.

| α | half-life | as % of a height-3.8 tree |
|---|---|---|
| 0.1 | 6.93 | 182% |
| 0.3 | 2.31 | 61% |
| 1.0 | 0.69 | 18% |
| 3.0 | 0.23 | **6%** |
| 10.0 | 0.07 | 2% |

**Recommendation:** choose α relative to tree height, targeting α·H ≈ 0.5–3 (half-life between
roughly a quarter and twice the tree height) so that the OU classes are distinguishable from
both BM and white noise. For these trees that means **α in the 0.1–0.5 range**, not 3.

> ### CORRECTION (added later the same day)
>
> **The recommendation immediately above is right for held-out prediction and wrong for model
> classification.** The two objectives have opposite difficulty gradients and were conflated here:
>
> | | α·H → 0 | α·H → ∞ |
> |---|---|---|
> | **classification** (BM vs OU) | **impossible** — OU converges to BM | easy — OU is white noise, BM is not |
> | **held-out prediction** of OU genes | easy — strong autocorrelation to borrow from | **impossible** — nothing left to predict |
>
> Acting on the α = 0.1–0.5 advice dropped classification accuracy to 0.52. Use α·H ≈ 0.5–3 for
> the dropout/prediction analysis, and a **higher** α·H for the classification benchmark — or,
> better, sweep α·H and report accuracy as a curve, which is what `runSCOUT.identifiability`
> in `R/identifiability.R` now does. See §4c.

**Second, subtler issue: the benchmark's difficulty drifts with `ncells`.** Height rises from
4.76 to 7.55 between n = 50 and n = 1000, so at fixed α a larger simulation is a *harder*
problem. Any comparison across cell numbers is confounded unless α is scaled with height or the
height is fixed. Note that setting `scale_tree = TRUE` at fit time does **not** fix this — it
rescales the fitted α by H but does not change the information content of data simulated on the
unscaled tree. The notebook is also inconsistent here: cell 9 has `scale_tree = TRUE` commented
out while cell 11 uses it.

**Scope caveat worth stating in the paper.** These are clean, ultrametric, correctly-resolved
trees. Real scLT trees reconstructed from character matrices (e.g. the `casNJ` trees used
elsewhere in the package) are typically non-ultrametric, contain polytomies, and carry
reconstruction error. A birth-death tree is the right *idealised* arm of a tree-robustness
analysis, but conclusions from it are an upper bound on real-data performance, not a
substitute for it. This also matters for the clade-dropout arm specifically: the rank-one
collapse makes within-clade predictions exactly constant only on an ultrametric tree, so on a
reconstructed tree that arm will behave slightly differently.

## 4c. Why classification accuracy regressed vs the TedSim era

Follow-up investigation of `examples/260804_sim_test_2.ipynb`, which reported 0.52 and 0.36.

**The 0.36 run was misconfigured, not a model result.** `scout.res.evf` read
`..._OUEVF_counts.csv` (the counts matrix) rather than `..._OUEVF_evf_counts.csv`, with
`normalize = FALSE`. Raw beta-Poisson counts span 0–41,342 and **12% of values exceed 1000, the
hard upper bound on `theta`** (`R/SCOUT_EM.R:214`), so no model can place its optimum near the
data mean, every fit is poor, and BM1 wins on parsimony — 86/90 genes called BM1. It also passed
the `260804_...nwk` tree with `260805_...` data. Both fixed; the notebook now derives every path
from one `prefix` variable.

**The 0.52 run is the real finding: every previously-working tree had unit branch lengths.**

| tree | branch lengths | height | ultrametric | α·H at the α used |
|---|---|---|---|---|
| TedSim vignette | all = 1 | 7.0 | yes | **21** (a=3) |
| casNJ n256 | **absent** → `formatSCOUT:1056` fills 1s | 9.0 | **no** (depths 6–9) | **27** (a=3) |
| hbd_rev (current) | real Yule, 0.006–2.34 | 3.80 | yes | **1.9** (a=0.5) |

A 14× drop in α·H, from the tree getting shorter *and* α being lowered on the §4b advice. Two
signals were lost at once:

1. **Correlation shape.** At α·H = 1.9 the OU and BM correlation structures still agree to
   `cor_shape = 0.886` on this tree, so AIC's parsimony preference tips OU1 genes into BM1 —
   exactly the observed pattern (16/30 OU1 and 17/30 OUM → BM1, BM1 itself fine at 26/30).
2. **Tip variance.** Under BM the variance of a tip is σ²·depth; under OU it is stationary. On a
   tree with varying tip depths that is a direct discriminator. The casNJ tree had it
   (`cv_var_bm = 0.106`); an **ultrametric tree has none of it** (`cv_var_bm = 0`), leaving only
   signal 1.

Both are measurable without any model fitting via `identifiability_reference` in
`R/identifiability.R`. Evaluated on the **real** trees at the α each was actually run with, the
whole regression is visible in one table:

| tree | α·H | `cor_shape` (BM vs OU) | `cv_var_bm` (2nd signal) |
|---|---|---|---|
| casNJ n256 | 27.0 | **0.282** | **0.068** |
| TedSim vignette | 21.0 | **0.356** | 0.000 |
| hbd_rev n100 (current) | 1.9 | **0.886** | 0.000 |

The BM and OU correlation structures now agree at 0.886 where they previously agreed at 0.28–0.36
— roughly a threefold increase in similarity — and the casNJ tree's tip-variance discriminator is
gone entirely. Neither number requires fitting a single model.

`runSCOUT.identifiability` sweeps α·H across branch-length variants of a fixed topology
(`yule_ultra`, `unit_ultra`, `unit_vary`, `yule_vary`), all rescaled to a common height so tree
shape is never confounded with tree height.

### Reduced sweep result (2 arms × 3 α·H × 2 matrices, 30 genes, 1 replicate, 46 min on 1 core)

| matrix | arm | α·H = 1.9 | α·H = 5 | α·H = 27 |
|---|---|---|---|---|
| **evf** | yule_ultra | **0.50** | **0.87** | 0.80 |
| **evf** | unit_vary | 0.73 | 0.90 | 0.90 |
| **counts** | yule_ultra | 0.60 | 0.57 | 0.60 |
| **counts** | unit_vary | 0.63 | 0.67 | 0.63 |

30 genes per cell gives ±0.15–0.19 binomial CIs, so only the largest contrasts are resolved.

**Established.** On the EVF matrix α·H is a strong driver: `yule_ultra` goes 0.50 → 0.87 between
α·H = 1.9 and 5, and those CIs do not overlap ([0.31, 0.69] vs [0.69, 0.96]). The empirical curve
tracks the analytic `cor_shape` (0.886 → 0.565) as predicted, which is the check that the harness
is measuring what it claims to.

**Suggestive, not resolved.** At α·H = 1.9 the depth-varying arm beats the ultrametric arm on EVF
(0.73 vs 0.50), the direction predicted by `cv_var_bm` = 0.113 vs 0, but the CIs overlap. Needs
the full grid with replicates.

**Revision to the diagnosis — this matters.** On the **counts** matrix α·H has essentially no
effect: 0.57–0.67 across a 14-fold range of α·H, on both arms. The α·H story is therefore *not*
the binding constraint for the counts benchmark, and **raising α will not repair it**. The
likely constraint is the counts pipeline itself: `Get_params` (`R/sim_utilities.R:100-105`)
replaces the EVF values with order statistics drawn from a fixed reference density — a global
rank transform that discards the variance signature separating BM from OU — and the beta-Poisson
draw then adds heavy shot noise on top. The gene-effect mixing also uses random-signed weights,
so a gene's regime offset scales with `sum(w_i)` and can cancel to near zero, which blunts OUM
specifically.

**Still open.** Why the historical counts runs scored well. The leading untested candidate is
**tip count**: those runs used n = 256, this sweep used n = 100, and discrimination power scales
with the number of tips. That should be the next thing varied.

**Practical recommendation.** Use the EVF matrix for the model-selection benchmark, where the
signal is real and α·H behaves as theory predicts. Report counts separately as the realistic-noise
case, and do not tune α expecting it to help there.

**Full production grid** (≈30 min on 32 cores; the dev box used here has 1 core, so only a
reduced grid was run):

```r
res <- runSCOUT.identifiability(
    tree     = simtree, results_dir = './identifiability',
    alphaH   = c(0.5, 1, 2, 5, 10, 20, 30),
    arms     = c('yule_ultra','unit_ultra','unit_vary','yule_vary'),
    matrices = c('evf','counts'),
    height = 7, ngenes = 30, nevfs = 30, n_replicates = 2,
    randseed = 42, cores = 32)
```

## 5. Verification performed

**Unit tests** — `tests/testthat/`, 89 assertions, 1.6 s:

- reduction to `e_step` when conditioning on all tips (BM1/OU1/OUM): agreement **5 × 10⁻¹⁵**
- single held-out tip vs brute-force `solve()`: **4 × 10⁻¹⁵**
- predictions invariant to perturbing the masked observations by +1000: exactly 0
- output order follows `held_idx` order; duplicates rejected
- true-parameter conditional beats the mean baseline and gives **coverage 0.95**
- `V_ho` rank 1 (clade) vs full rank (random)
- `sample_clade` respects size bounds, keeps all regimes, never picks a root child
- `descendant_tips` agrees with `ape::extract.clade`
- `compute_W_matrix` errors on an unmatched regime
- `OUwie.sim.edited` propagates `theta0` on a non-preordered tree

**`compute_VCV` rewrite** — `max |old − new| = 0.000e+00` across n ∈ {20, 120, 400},
α ∈ {0.05, 1, 7}, both root parameterisations; `shared_lengths` bit-identical. 8× faster at
n = 400.

**End-to-end** — 60-tip tree, 2 replicates, 6 genes × 3 regimes, 1.8 min on 4 cores:

| arm | masked | prior | global | pearson | coverage |
|---|---|---|---|---|---|
| random | **0.558** | 1.230 | 1.136 | 0.66 | 0.90 |
| clade | **0.883** | 1.211 | 1.266 | ~0 (rank-1) | 0.88 |

`masked ≈ oracle` throughout → masking did not degrade the parameter estimates. Model selection
agreed with the full-tree fit on 23/24 genes.

**Breakdown by generating model** (3 replicates, two α regimes, same simulation):

| α | arm | model | masked | prior | global | r |
|---|---|---|---|---|---|---|
| 0.3 | random | BM1 | 0.916 | 1.836 | 1.712 | 0.83 |
| 0.3 | random | OU1 | 0.812 | 1.121 | 1.088 | 0.47 |
| 0.3 | random | OUM | 0.736 | 0.986 | 1.001 | 0.70 |
| 3.0 | random | BM1 | 0.969 | 1.458 | 1.470 | 0.58 |
| 3.0 | random | OU1 | 0.385 | 0.392 | 0.391 | 0.27 |
| 3.0 | random | OUM | 0.691 | 0.873 | 0.926 | 0.69 |
| 0.3 | clade | BM1 | 1.601 | 2.296 | 2.193 | −0.12 |
| 3.0 | clade | BM1 | 1.007 | 1.633 | 1.491 | 0.01 |

---

## 6. Things to watch

- **Choose `a` deliberately.** At α = 3 on a height-4.4 tree, OU1 is nearly white noise
  (correlation ~ exp(−3d)) and there is genuinely almost nothing to recover — masked 0.385 vs
  baseline 0.392. That is a correct result, not a failure, but it makes a weak figure. Lower α
  gives a much clearer demonstration.
- **`t0` needs headroom.** `m_step_tipfog` bounds θ ≥ 10⁻¹⁰, so genes with non-positive values
  cannot be fit. BM1 drifts with variance σ²·height, so at height 4.4 even `t0 = 3` leaves 7.6%
  of BM1 values negative. `runSCOUT.dropout` warns; the vignette uses `t0 = 8`.
- **`runSCOUT` re-plans the future cluster on every call**, so N replicates × 2 arms means
  2N + 1 cluster spin-ups. Works fine, but for large N it is wasted time. Would need a change
  inside `runSCOUT` to fix.
- **`formatSCOUT` fails with a single gene column** — `apply(meta[, gene_cols], 2, as.numeric)`
  errors when `gene_cols` has length 1. Hit this while writing tests; not fixed.
- **`W[, root.state]` indexes by factor level code when `root.state` is a factor**
  (`compute_W_matrix:332`, OUM with `fixed_root = FALSE`). Self-consistent within a run, so it
  does not affect these results, but it is fragile. Not touched.
- **Not fixed, flagged in the plan:** `simulate_test_data` returns only the kon slice of the
  latent (`eevfs[[i]][[1]]`); `generate_EvoCounts:119` hardcodes `randseed = 123`;
  `run_em`'s `history` schema declares `ll_ou`/`ll_noise` columns that are never produced.

---

## 7. How to re-run

```r
library(SCOUT)
res <- runSCOUT.dropout(
    tree         = simtree,        # newick path or phylo with $states
    results_dir  = './dropout_output',
    n_replicates = 10, clade_min = 5, clade_max = 20,
    arms         = c('clade','random'),
    ngenes = 10, nevfs = 10, a = 3, s = 1, t0 = 8, theta_step = 2,
    randseed = 42, cores = 8)
```

Pass `sim_res = <simulate_test_data output>` to reuse one simulation across replicates, and
`sim_key` to select an (α, σ²) combination. Tests: `testthat::test_dir("tests/testthat",
package = "SCOUT")`. Full walkthrough in `examples/CladeDropoutVignette.ipynb`.
