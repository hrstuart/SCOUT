# Leave-one-out validation of SCOUT-EM on simulated trees, with a mask-size sweep

*Results of the simulation dropout run `SCOUT_dropout_20260811_v2`, launched 2026-08-11 and
completed 2026-08-12 (6 h 46 min on 24 cores). Companion to `dropout_validation_methods.md`
(the estimator) and `2026-08-10_celegans_rpl_dropout_results.md` (the real-data application).
Unlike the *C. elegans* run, the latent truth is known here, so absolute error is
interpretable and the generating model of each trait is available for stratification.*

---

## Methods

### Tree and traits

A 256-tip ultrametric tree was simulated with `castor::generate_tree_hbd_reverse`
(ρ = 1, λ = 1, μ = 0, seed 1). Three discrete regimes were painted on it with a Markov model
whose transition matrix

```
P = [0.90 0.08 0.02; 0.10 0.80 0.10; 0.02 0.08 0.90],   Q = logm(P)
```

was chosen to give persistent, contiguous regime blocks rather than isolated tips. The
resulting tip-state counts were 106 / 61 / 89 for states 1 / 2 / 3; internal states were
reconstructed with `infer_anc`.

Latent expression was simulated once, up front, with `simulate_test_data`
(α = 3, σ = 1, t₀ = 8, θ-step = 2, 25 EVFs per model class, seed 42) and passed to every fit,
so mask choice is the only factor that varies across the run. This yields **75 latent traits**
in three equal blocks of 25 generated under BM1, OU1 and OUM respectively, giving a known
`true_model` for each. Predictions are scored against the noiseless EVFs, not against
beta-Poisson counts; `ngenes = 50` therefore affects only the counts matrix written for
provenance and has no effect on the fits.

### Masking design

Twenty replicates were run. Each replicate drew one monophyletic clade of 5–20 tips, and
**six masking arms** were scored against it:

| arm | masked tips | description |
|---|---|---|
| `clade` | 5–19 (mean 9.75) | the drawn monophyletic clade |
| `random` | 5–19 (mean 9.75) | scattered tips, size-matched to that replicate's clade |
| `random_05pct` | 13 | scattered, 5% of the tree |
| `random_10pct` | 26 | scattered, 10% |
| `random_15pct` | 38 | scattered, 15% |
| `random_20pct` | 51 | scattered, 20% |

The size-matched `random` arm is the like-for-like control for `clade` — the two differ only
in *where* the masked tips sit, not how many there are — while the four fixed-fraction arms
give a dose-response over mask size. Clade sizes cannot be dialled freely to serve that
purpose: enumerating `eligible_clades()` on this tree gives 251 eligible clades whose sizes
jump 43 → 54 → 70 → 99, so no clade exists near 51 tips. This is why the fraction arms are
random-only, and it required a new `random_frac` argument to `runSCOUT.dropout()` (see
*Implementation note* below).

Every draw is constrained to keep all three regimes represented among the retained tips, since
a regime with no retained tips has no estimable optimum.

### Fitting

For each mask, SCOUT-EM was refitted on the retained tips alone under BM1, OU1 and OUM with
λ₁ = λ₂ = 0.2, a free root, and a tip-fog prior of mean 0.2 and SD 0.1. Masked latent
expression was predicted by the conditional Gaussian expectation

  **E**[**Z**_h | **X**_o] = **W**_h**θ** + **V**_ho(**V**_oo + τ²**I**)⁻¹(**X**_o − **W**_o**θ**),
  **Var**[**Z**_h | **X**_o] = **V**_hh − **V**_ho(**V**_oo + τ²**I**)⁻¹**V**_oh.

Four predictors were compared: `masked` (parameters from the retained tips), `prior`
(**W**_h**θ** alone, so `masked` − `prior` isolates the phylogenetic covariance term),
`oracle` (full-tree parameters, still conditioned only on retained tips) and `global` (mean of
the retained tips). Reported results are filtered to the regime minimising AIC on the *masked*
tree, so model selection never saw the held-out tips. Masks were drawn `upfront` rather than
inline, to keep the mask sequence independent of the RNG consumed by `runSCOUT` in proportion
to trait count.

The run comprised 20 replicates × 6 arms × 75 traits × 3 regimes ≈ 27,000 fits (120 masked
fits plus one full-tree reference), producing 108,000 metric rows. **No fit failed, no regime
was lost, no prediction failed and no gene was dropped** (all four tallies 0 across all 120
masks), and realised mask sizes matched the design exactly.

---

## Results

### Refitting on the retained tips costs almost nothing

The central validation question is whether having to re-estimate parameters from a depleted
tree degrades prediction. It does not: `masked` tracks `oracle` to within 0.1–2.8% of RMSE in
every arm × generating-model cell.

| generating model | arm | `masked` | `oracle` | gap |
|---|---|---|---|---|
| BM1 | random (matched) | 0.781 | 0.779 | +0.3% |
| BM1 | random_20pct | 0.859 | 0.858 | +0.1% |
| OU1 | random (matched) | 0.371 | 0.369 | +0.5% |
| OU1 | random_20pct | 0.383 | 0.381 | +0.5% |
| OUM | random (matched) | 0.407 | 0.396 | +2.8% |
| OUM | random_20pct | 0.424 | 0.413 | +2.7% |
| BM1 | clade | 1.317 | 1.316 | +0.1% |

The gap does not widen with mask fraction. Removing a fifth of the tree leaves enough
information to recover parameters that predict as well as the full-tree parameters do.

### Whole clades are much harder than scattered tips at the same mask size

Comparing `clade` against the size-matched `random` arm — identical mask sizes,
replicate by replicate, mean 9.75 tips — isolates the effect of mask geometry:

| arm | predictor | RMSE | MAE | \|mean error\| | Pearson *r* | coverage | coverage (obs.) |
|---|---|---|---|---|---|---|---|
| clade | oracle | 0.709 | 0.607 | 0.395 | 0.366 † | 0.699 | 0.857 |
| clade | **masked** | **0.713** | **0.611** | **0.402** | **0.318 †** | **0.704** | **0.872** |
| clade | prior | 0.904 | 0.803 | 0.626 | 0.722 † | — | — |
| clade | global | 1.303 | 1.188 | 1.000 | — | — | — |
| random | oracle | 0.515 | 0.409 | 0.146 | 0.733 | 0.699 | 0.869 |
| random | **masked** | **0.520** | **0.413** | **0.148** | **0.731** | **0.705** | **0.886** |
| random | prior | 0.956 | 0.797 | 0.326 | 0.848 | — | — |
| random | global | 1.342 | 1.146 | 0.417 | — | — | — |

† *See "Pearson is not interpretable in the clade arm" below.*

Masking a clade costs **37% more RMSE than masking the same number of scattered tips**
(0.713 vs 0.520). The mechanism is visible in the `prior` row: for scattered tips the
phylogenetic term nearly halves the error relative to the regime optimum alone
(0.956 → 0.520), whereas for a clade it recovers far less (0.904 → 0.713). A monophyletic
group joins the rest of the tree through a single branch, so the observed tips constrain it
only through that one connection.

### Mask size matters much less than mask geometry

Across the fraction sweep, error rises only slightly as the masked fraction quadruples:

| arm | masked tips | `masked` | `prior` | `global` | `masked` − `prior` |
|---|---|---|---|---|---|
| random_05pct | 13 | 0.525 | 0.976 | 1.359 | 0.451 |
| random_10pct | 26 | 0.548 | 0.978 | 1.377 | 0.430 |
| random_15pct | 38 | 0.543 | 0.971 | 1.375 | 0.428 |
| random_20pct | 51 | 0.555 | 0.985 | 1.379 | 0.430 |

RMSE rises 5.7% from 5% to 20% masking, and the advantage over the regime optimum is
essentially flat. Removing a fifth of the tips at random leaves the retained tips informative
about the removed ones. Note that the 37% clade-vs-random penalty at ~4% masking is more than
six times the penalty for quadrupling the amount of scattered data removed.

(`global`'s |mean error| falls steeply across the sweep, from 0.42 to 0.16, but this is an
artifact: a larger random sample of held-out tips has a mean closer to the tree-wide mean.
Its RMSE, the honest number, is flat.)

### Stratified by generating model, the recovery is almost entirely a BM1 phenomenon

This is the result that matters for interpretation, and it is invisible in the pooled tables.

| generating model | arm | `masked` RMSE | `prior` RMSE | gain from covariance | Pearson *r* |
|---|---|---|---|---|---|
| **BM1** | random (matched) | 0.781 | 1.997 | **61%** | 0.895 |
| BM1 | random_20pct | 0.859 | 2.067 | 58% | 0.897 |
| BM1 | clade | 1.317 | 1.876 | 30% | −0.037 † |
| **OU1** | random (matched) | 0.371 | 0.392 | **5%** | 0.329 |
| OU1 | random_20pct | 0.383 | 0.402 | 5% | 0.335 |
| OU1 | clade | 0.398 | 0.396 | −1% | −0.044 † |
| **OUM** | random (matched) | 0.407 | 0.481 | **15%** | 0.968 |
| OUM | random_20pct | 0.424 | 0.485 | 13% | 0.970 |
| OUM | clade | 0.424 | 0.440 | 4% | 0.818 † |

Three distinct regimes of behaviour:

- **BM1 traits are where held-out prediction works.** The phylogenetic term cuts error by
  ~60% and Pearson *r* ≈ 0.90. Brownian motion has no stationary distribution, so the regime
  optimum is nearly uninformative (`prior` RMSE ≈ 2.0, worse than the tree-wide mean) and
  essentially all the signal is carried by the covariance with neighbouring tips.
- **OU1 traits carry almost no cross-tip signal.** `masked` beats `prior` by 5%, and the clade
  arm shows no gain at all. At α = 3 the process reaches stationarity long before the tree
  depth, so a tip is nearly independent of its neighbours and the best available prediction is
  the optimum itself. **This reproduces, on simulated data with known truth, the effect
  inferred from the fitted α ≈ 1 in the *C. elegans* RPL analysis** — there it was a
  conjecture from the fitted half-life; here it is demonstrated against the generating model.
- **OUM traits look excellent by correlation and mediocre by RMSE.** Pearson *r* = 0.97 for
  `masked` — but also 0.97 for `prior`, which uses no tip information at all. The correlation
  is driven almost entirely by between-regime spread (θ differs by 2 per regime step), so
  reporting *r* alone for a multi-optimum model would overstate what the estimator is doing.
  RMSE, which is not inflated by between-regime variance, shows a real but modest 15% gain.

The practical implication for the revision: **any pooled dropout statistic is a mixture whose
value depends on the BM/OU composition of the trait set**, and should be reported stratified.

### Model selection on the masked tree is as accurate as on the full tree

Accuracy of AIC model recovery, masked fit vs full-tree reference (500 gene × mask
observations per cell):

| generating model | clade | random | 05pct | 10pct | 15pct | 20pct | full tree |
|---|---|---|---|---|---|---|---|
| BM1 | 1.000 | 0.998 | 0.998 | 1.000 | 0.998 | 1.000 | 1.000 |
| OU1 | 0.824 | 0.822 | 0.812 | 0.828 | 0.802 | 0.790 | 0.840 |
| OUM | 0.972 | 0.944 | 0.942 | 0.946 | 0.964 | 0.950 | 0.880 |

BM1 is recovered essentially perfectly. OU1 recovery sits at ~0.80–0.84 and degrades mildly
with mask fraction (0.822 → 0.790 from matched to 20%); its misclassifications go to OUM, an
over-specification that costs little in prediction. OUM recovery on the masked tree
(0.94–0.97) is *above* the full-tree figure (0.88), which is a small-sample AIC effect rather
than an improvement — the penalty term bites harder on the smaller retained sample and
disfavours OU1 for genuinely multi-optimum traits.

The key point is that model selection does not collapse under masking, so the AIC-filtered
results above are not an artifact of degraded selection.

### Prediction intervals are well calibrated under BM1 and OU1, and badly miscalibrated under OUM

Empirical coverage of nominal 95% intervals, against the latent truth:

| generating model | `masked` coverage | `oracle` coverage | mean `pred_sd` | actual error SD | ratio |
|---|---|---|---|---|---|
| BM1 | 0.885 | 0.884 | 0.665 | 0.908 | 1.37 |
| OU1 | 0.906 | 0.910 | 0.331 | 0.381 | 1.15 |
| OUM | **0.311** | **0.297** | 0.077 | 0.399 | **5.18** |

BM1 and OU1 intervals are mildly narrow (89–91% against a nominal 95%). OUM intervals are
wrong by a factor of five: the model reports a conditional SD of 0.077 while errors actually
scatter with SD 0.399.

Two observations constrain the explanation. First, the **`oracle` under-covers just as badly
(0.297)**, so this is not error in the masked refit — it is present with full-tree parameters.
Second, the point predictions are **unbiased**: mean error is 0.002, and the spread of
per-gene bias is 0.010 against a within-gene error SD of 0.399. The estimator gets the mean
right and the variance wrong. When an OUM-generated trait happens to be fitted under BM1
(141 of 3,000 cases), coverage returns to 0.887, confirming that the miscalibration is
specific to the OUM variance calculation rather than to the traits themselves.

#### The mechanism: OUM fits move variance out of the process and into tip fog

Inspecting the fitted parameters identifies the cause. Median estimates for traits fitted
under their own generating model (the reported `sigma` is the variance rate σ², which the OU1
row confirms):

| generating model | α̂ | σ̂² | τ̂ | implied process SD √(σ̂²/2α̂) | implied total SD |
|---|---|---|---|---|---|
| OU1 | 1.24 | 0.405 | 0.077 | **0.407** | 0.414 |
| OUM | 1.46 | **0.010** | **0.286** | **0.059** | 0.294 |
| *truth* | 3 | 1 | 0 | **0.408** | 0.408 |

The OU1 fit recovers the process SD essentially exactly (0.407 against a true 0.408) and
correctly puts τ̂ near zero — the EVF targets are noiseless, so there is no tip fog to find.

The OUM fit does not. It drives **σ̂² to 0.010 against a true 1**, a hundredfold
under-estimate, and compensates by inflating **τ̂ to 0.286 against a true 0**. The model
concludes that the evolutionary process barely moves and that nearly all the observed scatter
is measurement noise. Total variance is roughly preserved (0.294 against 0.408), which is why
point predictions remain unbiased and RMSE stays competitive — under either story, the best
guess for a masked tip is its regime optimum. Only the *split* between process and noise is
wrong.

That split is exactly what the latent interval depends on. `Var[Z_h | X_o]` counts the process
variance alone, so it comes out at 0.077 instead of ~0.40. Adding τ̂² back gives
`coverage_obs` = 0.807, which is why this stayed invisible in the *C. elegans* run: real-data
scoring uses `coverage_obs` because the target carries tip fog. **In simulation the target is
the noiseless EVF, so `coverage` is the correct column and `coverage_obs` is the misleading
one** — the reverse of the real-data convention.

Note that α̂ ≈ 1.24–1.46 against a true α = 3 in *both* OU fits, while the OU1 process SD is
nonetheless right. α and σ² are only weakly identified separately on an ultrametric tree; what
the data constrain is their ratio. OU1 lands on the correct point of that ridge, OUM slides
down it to a near-degenerate corner (σ̂² = 0.010 is close to a boundary, suggesting the
optimiser is hitting a floor rather than finding an interior optimum). Why the extra freedom
of three optima induces this is not yet established.

#### Scope: confined to one cell, and absent from real data

The collapse is not a general property of the OUM implementation. Fitted σ̂² across every
generating-model × fitted-regime cell, with the *C. elegans* real-data run for comparison:

| data | fitted regime | σ̂² | τ̂ | fraction with σ̂² < 0.05 |
|---|---|---|---|---|
| OUM (sim) | **OUM** | **0.010** | **0.286** | **100% of 3000** |
| OU1 (sim) | OUM | 0.403 | 0.077 | 0% |
| BM1 (sim) | OUM | 1.214 | 0.374 | 0% |
| OUM (sim) | OU1 | 0.861 | 0.326 | 0% |
| *C. elegans* counts | OU1 / OUalt / OUl | 0.392 | 0.106 | **0 of 4000** |

OU1-generated traits fit perfectly well *under the OUM model*, so the OUM code path is not
itself broken. And on the real *C. elegans* data — 391 cells of log-normalised counts — not a
single fit of 4,000 collapses, with per-regime coverage of 0.933 (OUalt) and 0.927 (OUl)
against a nominal 95%. **The real-data results are unaffected.**

The failure therefore requires something specific to this simulation's OUM configuration. Two
candidates, not yet separated:

- **Regime signal is far stronger than in real data.** θ-step = 2 against a process SD of
  0.408 is roughly 5:1. Once the fit has the correct three optima, the residual is pure
  stationary scatter — and at a half-life of 5% of tree depth that scatter is very nearly iid
  across tips, which is exactly what τ models. The two components become near-degenerate.
  When the fitted optima are *wrong* (OUM data under OU1) the residual retains the contiguous
  between-regime block structure, which is not iid, and σ̂² is forced up to 0.861 instead.
- **The tree is exactly ultrametric** — all 256 root-to-tip depths are 4.739. Tips at differing
  depths are what let a fit separate accumulated process variance from constant observation
  noise; the *C. elegans* tree has depths of 5–11, which breaks the degeneracy.

**Consequences, in scope order:**

1. Do not quote `coverage` or the fitted variance components from the OUM arm of *this
   simulation*. Point predictions, RMSE, the clade-vs-random contrast, the BM/OU
   stratification and model selection are unaffected — none depend on the process/noise split.
2. The manuscript's real-data results need no re-checking on this account.
3. The residual risk is methodological: **an ultrametric tree with well-separated optima is the
   standard way a simulation-based calibration check gets built**, so any future "SCOUT's
   intervals are calibrated" claim resting on simulation would fail here — possibly for reasons
   that are an artifact of the simulation design rather than a property of the estimator.
   Worth resolving before such a claim is made, and cheap to probe: refit with a
   non-ultrametric tree, or with θ-step reduced to ~0.5, and see whether σ̂² recovers.
   Planned in `2_dropout_output/plans/260812_sigma_collapse_probe.md`.

#### Follow-up (2026-08-12): it is the τ path, not the bound

σ̂² = 0.0100 sits suspiciously close to the optimiser's floor (`sigma_min <- 0.01`,
`R/SCOUT_EM.R:220`), but a direct test rules the bound out:

| fit | σ̂² (OUM traits) | implied process SD | truth |
|---|---|---|---|
| `sigma_min = 0.01` (shipped) | 0.010 | 0.059 | 0.408 |
| `sigma_min = 1e-6` (patched) | 0.010 — **bit-identical** | 0.059 | 0.408 |
| **`skipTau = TRUE`** (no tip fog) | **0.53–0.64** | **0.42–0.46** ✓ | 0.408 |

Lowering the floor by four orders of magnitude changed nothing to any decimal place, so
σ̂² ≈ 0.010 is a genuine *interior* optimum rather than a clamp. (The patch was verified live by
deparsing the loaded `m_step_tipfog`, not by reading the file.) Removing tip fog from the model
entirely, however, recovers the process SD almost exactly.

**The collapse is therefore a τ-path phenomenon.** With τ free, the OUM likelihood genuinely
prefers σ² ≈ 0.01 with τ ≈ 0.29 over σ² ≈ 0.55 with τ ≈ 0. At a half-life of 5% of tree depth
the stationary OU scatter is very nearly iid across tips, which is exactly what τ models, and
the fit chooses the tip-fog explanation.

Why OUM and not OU1, on the same tree with the same α and τ equally free: under OUM the three
regime optima already account for the between-regime spread, so the model can afford to call
all remaining within-regime scatter noise. Under OU1 there is a single optimum, so that spread
*must* be carried by the process and σ² is held up. The two cross-fits corroborate it — OUM
data under OU1 inflates σ̂² to 0.861 (the process absorbing regime differences), OU1 data under
OUM sits correctly at 0.403 (no regime structure to exploit).

Script `scripts/260812_sigma_min_test.R`; results
`output/260812_sigma_min_test_{baseline,patched}.csv`. This makes the factorial probe in
`plans/260812_sigma_collapse_probe.md` more relevant rather than less, with the question
restated as *under what conditions does the τ path steal the process variance*.

### Pearson is not interpretable in the clade arm

For a monophyletic mask, **V**_ho is rank one — every masked tip connects to the observed
tips through the same branch — so under a single-optimum model the predictions within the
clade are constant and correlation is undefined. The run bears this out: Pearson is `NA` for
**50.7% of clade masks and 0% of random masks**, rising to **99.5% of clade masks whose
selected regime is OU1**. Where it is not `NA` under a single-optimum model, `sd(pred)` is
non-zero only through floating-point residue, and the resulting values are noise around zero —
hence the −0.037 and −0.044 in the tables above.

The exception is OUM, where regime variation inside the clade makes **W**_h**θ** vary across
masked tips and *r* = 0.82 is meaningful — though, as noted, it then measures regime
assignment rather than tip-level prediction.

**Report `mean_error` or RMSE for the clade arm, never Pearson.**

---

## Implementation note

The random arm of `runSCOUT.dropout()` was previously hard-coded to be size-matched to each
replicate's clade, which made the fraction sweep impossible. A backwards-compatible
`random_frac` argument was added:

```r
runSCOUT.dropout(..., random_frac = c(NA, 0.05, 0.10, 0.15, 0.20))
```

`NULL` (the default) reproduces the previous behaviour exactly; an `NA` entry means the
size-matched arm, so the matched control and the fraction sweep can run together; numeric
entries in (0, 1) produce arms named `random_05pct`, `random_10pct` and so on. Validation
rejects out-of-range values, duplicates and fractions that would leave too few tips to fit on.

Two minor issues were noticed and not acted on:

- `simulate_test_data` emits a stray column literally named `OUM` alongside `OUM_1…OUM_25` in
  the EVF matrix. The fitter correctly ignores it (75 traits scored), but the runner script's
  manifest counts it, so `n_traits` reads 76 where 75 were scored.
- `runSCOUT.dropout` still does not forward `param_mode` / `signal_gain` / `burst_scale` to
  its internal `simulate_test_data` call. Irrelevant here, since dropout scores EVFs rather
  than counts, but it will matter for counts-based work.

## Reproducing

```bash
cd /dartfs/rc/lab/M/McKennaLab/projects/hannah/OU/revisions_analysis/revision_v2/2_dropout_output/scripts
Rscript 260811_sim_v2_loo_analysis.R          # SCOUT_SMOKE=1 for a ~6 min 64-tip version
```

Every parameter, both seeds, package versions and `sessionInfo()` are recorded in
`output/SCOUT_dropout_20260811_v2_manifest.{rds,csv}`. Tree seed 1, simulation seed 42,
SCOUT 0.1.0, R 4.5.3.
