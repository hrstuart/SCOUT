# Held-out validation of SCOUT-EM on the *C. elegans* random-precise-lineage dataset

*Results of the production dropout run of 2026-08-07. Companion to
`dropout_validation_methods.md`, which describes the estimator and the simulation-based
verification; the methods subsection below covers only what is specific to the real-data
application.*

---

## Methods (real-data application)

Held-out validation was applied to the *C. elegans* random-precise-lineage dataset
(`0_RPL_MVk13qbQ`; 391 cells, lineage tree non-ultrametric with root-to-tip depths of 5–11
branch units, median 9). The expression matrix was already log-normalised, so the internal
normalisation step was disabled (`normalize = FALSE`); applying it would have log-transformed
the matrix a second time.

Because the cost of the procedure is linear in the number of genes (≈67 CPU-seconds per
gene × regime fit at this tree size), the full 3,980-gene matrix was not tractable. Genes were
therefore restricted to the 500 most variable, from which 100 were sampled at random
(seed 42); these have per-gene variances of 0.203–2.897 against a matrix-wide median of 0.071,
and 62.2% zeros against 82.7% matrix-wide.

Two masking arms were run. In the **clade** arm a whole monophyletic group was removed; in the
**random** arm the same number of cells was removed at random from across the tree. Candidate
clades were enumerated with `eligible_clades()` subject to three constraints: 10–50 tips, not a
child of the root, and every regime still represented among the retained cells (a regime with
no retained tips has no estimable optimum). Forty-two clades satisfied these constraints, with
sizes 10–36 — the tree's topology is such that no clade falls between 37 and 56 tips. Twenty
were selected by keeping all clades of ≥25 tips, which are scarce, and drawing one from each of
ten equal-count size bins below that threshold, so that the mask-size gradient survived the
subsample. The selected masks span 10–36 cells, that is 2.6%–9.2% of the tree (mean 5.8%).

For each mask, SCOUT-EM was refitted on the retained cells alone under four regimes (BM1, OU1,
OUalt, OUl) with λ₁ = λ₂ = 0.2, a stationary root, and a tip-fog prior of mean 0.2 and standard
deviation 0.1 — the settings established in the replicate analyses. The masked cells' latent
expression was then predicted from the retained cells by the conditional Gaussian expectation

  **E**[**Z**_h | **X**_o] = **W**_h**θ** + **V**_ho(**V**_oo + τ²**I**)⁻¹(**X**_o − **W**_o**θ**),
  **Var**[**Z**_h | **X**_o] = **V**_hh − **V**_ho(**V**_oo + τ²**I**)⁻¹**V**_oh.

Four predictors were compared: `masked`, the estimator above using parameters fitted to the
retained cells only; `prior`, the regime optimum **W**_h**θ** alone, so that `masked` − `prior`
isolates the contribution of the phylogenetic covariance; `oracle`, the same estimator using
full-tree parameters but still conditioned only on retained cells; and `global`, the mean of the
retained cells with no model. Per gene and mask, the reported model is the one minimising AIC on
the **masked** tree, so model selection never saw the held-out cells.

Since no latent truth exists for real data, predictions were scored against the observed values
at the masked cells, which carry tip fog; root-mean-square error therefore has an irreducible
floor and the interpretation rests on the comparison between predictors rather than on absolute
error. Interval coverage is reported both against the latent posterior variance (`coverage`) and
against that variance widened by τ² (`coverage_obs`), the latter being the appropriate nominal
95% target for an observed-scale comparison.

The run comprised 20 masks × 2 arms × 100 genes × 4 regimes = 16,000 fits and completed in
16 h 12 min on 16 cores. No fit failed, no regime was lost, and no gene was dropped.

## Results

### Recovery is substantial for scattered cells and marginal for whole clades

| arm | predictor | RMSE | MAE | \|mean error\| | Pearson *r* | coverage | coverage (obs.) |
|---|---|---|---|---|---|---|---|
| clade | oracle | 0.558 | 0.467 | 0.162 | 0.001 | 0.929 | 0.936 |
| clade | **masked** | **0.563** | **0.471** | **0.176** | **0.008** | **0.924** | **0.930** |
| clade | prior | 0.578 | 0.485 | 0.199 | 0.044 | — | — |
| clade | global | 0.591 | 0.505 | 0.220 | — | — | — |
| random | oracle | 0.521 | 0.390 | 0.095 | 0.356 | 0.920 | 0.927 |
| random | **masked** | **0.524** | **0.392** | **0.098** | **0.350** | **0.917** | **0.925** |
| random | prior | 0.583 | 0.462 | 0.111 | 0.270 | — | — |
| random | global | 0.606 | 0.491 | 0.115 | — | — | — |

Expressed as a pooled skill score, 1 − MSE_masked/MSE_baseline:

| arm | vs. `global` | vs. `prior` | vs. `oracle` |
|---|---|---|---|
| clade | 0.045 | 0.022 | −0.016 |
| random | **0.282** | **0.237** | −0.011 |

In the random arm the estimator explains 28% of the error variance of the model-free baseline
and 24% beyond the regime optimum, and beats `global` on 73.4% of gene × mask combinations. In
the clade arm it explains 4.5% and 2.2% respectively, and the median per-gene gain over `prior`
is 0.0002% — that is, for a typical gene the phylogenetic covariance contributes nothing at all
once a whole clade is removed.

### Masking barely degrades parameter estimation

`masked` trails `oracle` by 1.6% of error variance in the clade arm and 1.1% in the random arm.
Refitting on the retained cells therefore recovers very nearly the parameters obtained from the
complete tree, and the gap between the two arms is a property of what the retained cells can say
about the held-out ones, not of estimation quality. Consistently, the AIC-selected model on the
masked tree matched the full-tree selection for 94.9% (clade) and 94.9% (random) of
gene × mask combinations.

Posterior intervals are close to calibrated but mildly overconfident: nominal 95% intervals
achieved 92.5%–93.0% observed-scale coverage.

### The two arms differ because of mean reversion, not clade size

The clade-arm result is not a size effect. Across masks of 10 to 36 cells the skill over
`global` fluctuates between −0.9% and 15% with no trend, and in the random arm the gain is flat
at 10%–17% across the same range. Splitting instead by the selected model class is decisive:

| arm | selected model | skill vs. `prior` | skill vs. `global` |
|---|---|---|---|
| clade | BM1 | 0.149 | 0.185 |
| clade | OU (OU1/OUalt/OUl) | **0.00003** | 0.047 |
| random | BM1 | **0.456** | 0.463 |
| random | OU (OU1/OUalt/OUl) | 0.020 | 0.128 |

For genes assigned an OU model, removing a clade leaves *literally* no recoverable information
beyond the regime optimum: the median absolute difference between the `masked` and `prior`
predictions is 2.6 × 10⁻⁵, and the standard deviation of the predictions across the masked cells
of a given gene is 4.3 × 10⁻⁶ against an observed standard deviation of 0.507. The conditional
expectation collapses onto a constant equal to the regime optimum.

This follows directly from the fitted selection strengths. The median α of the selected OU
models is 0.995, giving a phylogenetic half-life of log 2/α = 0.70 branch units against a median
root-to-tip depth of 9 — a memory of 7.7% of the tree's depth. Covariance between a masked cell
and a retained cell decays as exp(−α·*d*), so across a clade's stem branch and the path to its
nearest retained relative it is numerically zero. In the random arm each masked cell typically
retains a sister at short distance, which is why some signal survives there (skill 0.020 over
`prior`) but not much. Under BM1, where α is fixed at 10⁻¹⁰ and covariance is shared path length
that never decays, information transfer is large in both arms.

Two consequences follow. First, the near-zero Pearson correlation in the clade arm (*r* = 0.008)
is a structural property, not a diagnostic of poor fit: when the prediction vector is constant,
its correlation with anything is noise about zero. Mean error, not correlation, is the
informative statistic for that arm. Second, for an OU-selected gene the model's own claim is
that expression is pinned near its regime optimum with fast mean reversion, and the best
prediction for a held-out cell therefore *is* that optimum. The held-out test confirms the
estimator behaves as the fitted model says it should, with approximately calibrated intervals;
it does not indicate that phylogenetic borrowing has failed.

The practical implication for imputation is that SCOUT can meaningfully reconstruct expression
for cells scattered among retained relatives, particularly for BM-like genes, but should not be
expected to reconstruct an entire missing clade unless the gene is BM-like.

## Reproduction

```
Rscript revision_v2/3_celegans_rpl/scripts/260806_celegans_rpl_loo_runner.r
```

Outputs in `revision_v2/3_celegans_rpl/loo_output/`, prefixed `0_RPL_MVk13qbQ_LOO`:
`_genes.csv` and `_masked_clades.csv` (the sampled genes and selected clades),
`_20260807_{clades,predictions,metrics,params,model_select}.csv`, `_20260807_dropout.rds`,
`_headline.csv`, and `_full_fit.rds` (the cached full-tree reference pass, reusable via the
`full_fit =` argument to skip the 28-minute reference stage). Seeds: 42 throughout.
