# Methods: held-out validation of SCOUT-EM by tree masking

*Draft methods section, Genome Biology style. Simulation and masking parameters below are
those used in the verification runs of 2026-08-06; update the italicised values to match the
final production run before submission.*

---

## Phylogenetic Ornstein–Uhlenbeck model

SCOUT models the expression of a single gene across the tips of a single-cell lineage tree as
a latent Ornstein–Uhlenbeck (OU) process observed with independent measurement error. Writing
**Z** for the latent expression at the *n* tips and **X** for the observed values,

  **X** = **Z** + **ε**,  **Z** ~ N(**Wθ**, **V**),  **ε** ~ N(**0**, τ²**I**),

where **θ** collects the regime-specific optima, **W** is the *n* × *k* OU weight matrix giving
each tip's exposure to each regime along its root-to-tip path, and τ is the tip-fog standard
deviation. The covariance **V** follows the closed form of Ho and Ané [1]. Letting *t_i* denote
the root-to-tip distance of tip *i*, *s_ij* the root-to-most-recent-common-ancestor distance of
tips *i* and *j*, and *d_ij* = *t_i* + *t_j* − 2*s_ij*,

  *V_ij* = (σ²/2α)·exp(−α*d_ij*)·(1 − exp(−2α*s_ij*))  (fixed root)
  *V_ij* = (σ²/2α)·exp(−α*d_ij*)          (stationary root),

with α the strength of selection and σ² the diffusion variance. Three model classes were
fitted: BM1 (Brownian motion, α fixed at 10⁻¹⁰, single optimum), OU1 (single optimum) and OUM
(one optimum per regime).

Parameters were estimated by expectation–maximisation. The E-step computes the posterior mean
of the latent process given all tips, **E**[**Z**|**X**] = **Wθ** + **K**(**X** − **Wθ**) with
Kalman gain **K** = **V**(**V** + τ²**I**)⁻¹. The M-step maximises the OU log-likelihood of the
current posterior mean, evaluated by the three-point algorithm of Ho and Ané as implemented in
phylolm [2], over (α, σ², **θ**, τ) on the log scale using the subplex algorithm
(NLOPT_LN_SBPLX) in nloptr, with a maximum of 200 function evaluations and relative and
absolute tolerances of 10⁻⁸ and 10⁻¹⁰. Parameters were bounded to α ∈ [10⁻³, 10²],
σ² ∈ [10⁻², 10²], θ ∈ [10⁻¹⁰, 10³] and τ ∈ [0.05, 2]. Weakly informative penalties were applied
to α (Gamma, shape 2, rate 1) and to τ (Normal, mean 0.2·sd(**X**), standard deviation 0.1),
each scaled by *n* and by weights λ₂ = λ₁ = *0.2* so as to be commensurate with the
log-likelihood. EM was run for at most 100 iterations and stopped when the change in total
log-likelihood fell below both an absolute tolerance of 0.01 and a relative tolerance of 10⁻⁴,
with an additional early-stopping rule when the likelihood trace oscillated.

## Simulation of expression data on lineage trees

Lineage trees were simulated as reverse-time homogeneous birth–death trees using
`generate_tree_hbd_reverse` in castor [3] with *n = 100* tips, sampling fraction ρ = 1,
speciation rate λ = 1 and extinction rate μ = 0; the resulting trees are ultrametric. Discrete
cell-state regimes were painted onto the tips by simulating a three-state continuous-time
Markov model along the tree (`simulate_mk_model`, castor), with generator **Q** = log(**P**)
for a transition matrix **P** with diagonal 0.90, 0.80, 0.90 and off-diagonal mass distributed
as described in the accompanying vignette. Ancestral regimes at internal nodes were
reconstructed by maximum likelihood under an equal-rates model (`ace`, ape [4]).

For each of *nevfs = 10* independent latent factors and each of the three generative model
classes, tip values were simulated by numerical integration of the OU process down the tree.
BM1 factors used α = 10⁻¹⁰ and a single optimum θ = *t₀*; OU1 factors used α = *a* and a single
optimum θ = *t₀ + Δ*; OUM factors used α = *a* and regime optima spaced arithmetically,
θ_r = *t₀* + (r − 1)·*Δ* over the three regimes. All factors used σ² = *s* and root state
θ₀ = *t₀*. Simulations used *a = 3*, *s = 1*, *t₀ = 8* and *Δ = 2*.

Two constraints govern the choice of simulation parameters. First, because the estimation step
constrains regime optima to be positive (θ ≥ 10⁻¹⁰), *t₀* must be large enough that simulated
trajectories remain positive; for Brownian factors the tip variance is σ²·*h*, where *h* is the
tree height, so *t₀* was set to provide approximately three standard deviations of headroom, and
simulations were checked for non-positive values. Second, the selection strength must be chosen
relative to tree height rather than in absolute units. The phylogenetic half-life of an OU
process is log(2)/α, and trees generated under the above model have height increasing with tip
number (approximately 5.0 at *n* = 100 and 7.6 at *n* = 1000 for λ = 1); values of α for which
the half-life is a small fraction of the tree height yield tip values that are effectively
independent, rendering the OU classes indistinguishable from uncorrelated noise. Conversely, as
α·*h* approaches zero the OU covariance converges to the Brownian covariance and the two model
classes become unidentifiable. The two analyses reported here therefore occupy opposite ends of
this axis and require different settings: **held-out prediction** of an OU trait requires
residual phylogenetic autocorrelation and so needs α·*h* in the range 0.5–3, whereas **model-class
classification** requires visible mean reversion and degrades as α·*h* falls, being impossible in
the limit. Rather than fix a single value, classification accuracy is reported as a function of
α·*h* (see below), which exposes the identifiability floor directly.

> **Parameter values to finalise.** The verification runs used *a = 3* (α·*h* ≈ 15) for the
> classification benchmark and *a = 0.5* for the prediction analysis. Confirm the production
> values against the identifiability sweep before submission and update these numbers.

## Model-class identifiability

Because the discriminability of BM and OU depends jointly on α·*h* and on tree shape, and because
the trees used across the benchmark differ in both, classification accuracy was characterised as
a function of α·*h* rather than reported at a single operating point. Four branch-length variants
were constructed from a single fixed topology, so that tree shape is never confounded with
topology: the tree as generated (heterogeneous branch lengths, ultrametric); uniform internal
branch lengths with tips equalised to a common depth; uniform branch lengths with tips at varying
depths; and heterogeneous branch lengths with tips at varying depths. All four were rescaled to a
common height, making α·*h* comparable across arms.

Two mechanisms separate the model classes and were quantified analytically, without model
fitting. The first is the shape of the covariance decay, summarised as the correlation between the
off-diagonal entries of the Brownian and OU correlation matrices; this tends to unity as α·*h*
tends to zero, marking the region in which no method can classify. The second is the marginal tip
variance, which under Brownian motion is proportional to a tip's root-to-tip depth and under an
OU process is independent of it; on an ultrametric tree this signal is absent entirely, whereas on
a tree with varying tip depths it provides an additional discriminator. Both quantities were
computed from the covariance matrices directly and are reported alongside the empirical accuracy
curves.

These latent factors, rather than the derived count matrix, were used as ground truth. The
count matrix produced by the simulator applies a gene-effect mixing step followed by
beta-Poisson sampling, so individual count columns are mixtures of several latent factors and
carry no one-to-one correspondence with any single trajectory; the latent factor matrix by
contrast is the OU process itself, observed without noise, and is therefore directly
comparable in scale to the model's posterior estimate of **Z**.

## Held-out prediction of masked tips

Validation proceeds by masking a subset of tips, refitting the model on the reduced tree, and
predicting the masked tips from the retained tips alone. Partitioning the tips into a retained
set *o* and a held-out set *h*, and noting that the measurement error is independent of the
latent process so that Cov(**Z**_h, **X**_o) = **V**_ho and Var(**X**_o) = **V**_oo + τ²**I**,
the conditional distribution of the held-out latent expression is Gaussian with

  **E**[**Z**_h | **X**_o] = **W**_h**θ** + **V**_ho(**V**_oo + τ²**I**)⁻¹(**X**_o − **W**_o**θ**)
  Var[**Z**_h | **X**_o] = **V**_hh − **V**_ho(**V**_oo + τ²**I**)⁻¹**V**_oh.

This is the standard E-step restricted to the cross-covariance block; setting *h* = *o* = all
tips recovers **E**[**Z**|**X**] exactly, and this identity was used as a numerical check of
the implementation (agreement to within 5 × 10⁻¹⁵). Crucially, **X**_h enters neither the
conditioning set nor the parameter estimates, so the prediction is free of leakage.

Covariance and weight matrices were evaluated on the geometry of the *complete* tree but at the
parameters estimated from the *masked* tree. Systems were solved by Cholesky factorisation,
falling back to a Moore–Penrose pseudoinverse (corpcor) when the factorisation failed, with a
ridge of 10⁻⁸ added to the conditioning block for numerical stability. Because the latent
factors carry no observation noise, interval coverage was assessed against Var[**Z**_h|**X**_o]
rather than the predictive variance of a new observation.

## Masking schemes

Two masking schemes were applied to each replicate. In the **clade** scheme, all tips
descending from a single internal node were masked. In the **random** scheme, the same number
of tips was drawn uniformly at random from across the tree.

The two schemes are not interchangeable, and the distinction is a property of the model rather
than of the implementation. For a masked monophyletic clade, the most recent common ancestor of
any masked tip *i* and any retained tip *o* is the same node irrespective of *i*, so *s_io*
depends only on *o*. Substituting into the covariance above,

  *V_io* = (σ²/2α)·exp(−α(*t_i* + *t_o* − 2*s(o)*))·(1 − exp(−2α*s(o)*)) = exp(−α*t_i*)·*g(o)*,

for both root parameterisations; that is, **V**_ho is of rank one. The correction term
therefore collapses to a single scalar *c*, and **E**[**Z**_h|**X**_o]_i = (**W**_h**θ**)_i +
*c*·exp(−α*t_i*). On an ultrametric tree all *t_i* are equal and predictions within the masked
clade are constant up to the regime-weight term. This is the correct behaviour of a
well-specified model — a monophyletic clade communicates with the remainder of the tree through
a single branch — but it means that within-clade variation is not recoverable in principle, and
that per-tip error inside a masked clade is bounded below by the irreducible within-clade
variance. Clade masking was therefore scored on recovery of the clade's overall level and on
interval calibration, while per-tip recovery was assessed using the random scheme, for which
**V**_ho is of full rank. Both properties were verified numerically on each simulated tree.

Candidate clades were enumerated over all internal nodes and restricted to those with between
*5* and *20* descendant tips. Two classes of clade were excluded. Clades descending from either
child of the root were excluded because removing one collapses the root and shifts every
root-to-tip distance, which would place the masked-tree parameters on a different time scale
from the complete tree. Clades whose removal would eliminate a regime from the retained tips
were excluded because the optimum of an unrepresented regime is not estimable, leaving part of
the weight matrix undefined; the same constraint was applied to the random scheme. A clade was
then drawn uniformly from the remaining candidates. Masked trees were produced with `drop.tip`
(ape), which collapses resulting singleton nodes and so preserves the root-to-tip distances of
the retained tips.

Ancestral regimes were reconstructed afresh on each masked tree, and all three model classes
were refitted independently on it. Tree-height rescaling was disabled throughout, since
rescaling edge lengths by tree height would place α and σ² on tree-specific time scales and
invalidate the transfer of parameters from the masked tree to the complete tree.

## Baselines and evaluation

Each set of held-out predictions was compared against three reference predictors. The **prior**
predictor uses **W**_h**θ** alone, that is the fitted regime optima with no phylogenetic
borrowing from the retained tips; it is the null against which any claim of phylogenetic
information must be judged. The **oracle** predictor applies the same conditional estimator but
with parameters estimated on the complete tree, and so separates degradation of the parameter
estimates caused by masking from the intrinsic difficulty of the prediction. The **global**
predictor uses the mean of the retained observations. Performance was summarised by root mean
squared error, mean absolute error, bias, Pearson and Spearman correlation, the error in the
mean of the masked set, and empirical coverage of 95% posterior intervals, computed per
replicate, per gene and per fitted model class.

Model selection was assessed separately by comparing the AIC-preferred model class on the
complete and masked trees. AIC and AICc were computed from the log-likelihood at convergence,
the number of free parameters, and the number of tips actually retained in each masked tree.

*N = 10* masking replicates were performed, each drawing an independent clade and an
independent size-matched random set. All replicates were scored against a single simulated
dataset so that variation across replicates reflects the choice of masked tips rather than
simulation noise; random number generation was seeded once per run.

## Implementation and availability

All analyses were performed in R 4.4.2. The validation procedure is implemented as
`runSCOUT.dropout` in the SCOUT package and depends on ape 5.8.1, castor 1.8.3, phylolm 2.6.5,
nloptr 2.2.1, corpcor 1.6.10, paleotree 3.4.7 and future.apply 1.20.2. Fits across genes and
model classes were parallelised with future.apply. A regression test suite covering the
conditional estimator, the masking schemes and the simulator accompanies the implementation.

## References

1. Ho LST, Ané C. A linear-time algorithm for Gaussian and non-Gaussian trait evolution models.
   Syst Biol. 2014;63:397–408.
2. Ho LST, Ané C. phylolm: phylogenetic linear regression. R package.
3. Louca S, Doebeli M. Efficient comparative phylogenetics on large trees. Bioinformatics.
   2018;34:1053–5.
4. Paradis E, Schliep K. ape 5.0: an environment for modern phylogenetics and evolutionary
   analyses in R. Bioinformatics. 2019;35:526–8.
