# ALSFRS-R - handover

## 1. Where things are

| What | Where |
|--------------------------------|----------------------------------------|
| Code | branch `reorganise-r-directory`, pushed, 4 commits ahead of `master`, uncommitted work on top (§2) |
| Prior × backend results | [Google document](https://docs.google.com/document/d/1u-B9efjk9vD1P1Zp4iXq2D2AjKFOMyyJkPp3RLQEklM/edit?usp=sharing), partial results, needs to be rerun for true time estimates |
| s2z and nutpie exploration | [Google document](https://docs.google.com/document/d/1u-B9efjk9vD1P1Zp4iXq2D2AjKFOMyyJkPp3RLQEklM/edit?usp=sharing) |
| Registry data | `private_data/F_PROACT_ALSFRS.csv` The new version, you have it in your e-mail. |

## 2. Repository state relative to `master`

### 2.1 What is new

| path | what changed/is new |
|------------------------------------|------------------------------------|
| `R/brm_parallel.R` | nutpie backend compile branch, `backend` default fix |
| `R/sampling_parallel.R` | nutpie backend, cache fix, cache keys (§2.2) |
| `R/nutpie_to_stanfit.R` | nutpie to stanfit converter, does not depend on rstan internals |
| `analyses/Misa_sims/Sim1_wide.Rmd` | the last server run file with wide logit models |
| `analyses/Misa_sims/Sim1_long.Rmd` | the long-format model on the Sim1 replicates; contains `effect_hier_from_sim_study` (marginal estimand code) |
| `analyses/prior_backend/nutpie_stan_benchmark.R` | long-format model with cmdstanr and nutpier |
| `analyses/prior_backend/prior_backend_small.R` | prior × backend on Q5, Q12 and the twelve-item model; the script behind §4.5 |
| `analyses/prior_backend/prior_backend_small_report.R` | Claude's file that creates a report form prior_backend_small.R |
| `analyses/s2z` | s2z exploration on both wide and long models |

### 2.2 What changed in the main pipeline

**`R/sampling_parallel.R`**

- A `nutpier` backend. New arguments `nutpie_model`, `nutpie_args` (for example `list(adaptation = "low_rank")` or `list(num_warmup = 1000)`), `nutpie_converter` (defaults to `nutpie_draws_to_stanfit()`), `extra_cache_key`. brms/cmdstanr arguments are translated: `chains → num_chains`, `iter − warmup → num_draws`, `adapt_delta → target_accept`, `max_treedepth`, `init` (same meaning in both backends), `seed`; `store_divergences = TRUE` is forced; cmdstanr-only arguments are dropped by a whitelist against `formals(nutpie_sample)`. **`warmup` is deliberately not mapped to `num_warmup`**: nuts-rs uses its own default, 400 draws for diagonal adaptation and 800 for low-rank. Override with `nutpie_args$num_warmup`. A message says which is in force. Every per-time comparison of the two backends carries this asymmetry.
- **The fit cache never worked from commit `99ed820` (2026-05-18, "Update renv, use mori") until now.** `rlang::hash()` of a `mori::share()`d object is not reproducible, so the data half of every cache key was random and every run refitted everything. Hashes are now computed in the master process before sharing. Verified: two identical calls produce one cache file. The summary cache key also includes `deparse(body(summarise_fun))`, so editing the summary function invalidates summaries without discarding fits. Sampler arguments are deliberately **not** in the fit key, so two runs differing only in `iter`, `warmup`, `adapt_delta` or `seed` share a key and the first wins: use a separate `cache_dir` per sampler configuration.
- nutpie fits get their own cache namespace, `fit_nutpie-diag_<code>_<data>.rds` and `fit_nutpie-lowrank_…`, because they are also `stanfit` objects.
- Two pre-existing issues untouched: `chains` is translated to cmdstanr's deprecated `num_chains`, and `control$max_treedepth` to `max_depth`, which is not cmdstanr's argument name. No current call site passes `control`.

**`R/brm_parallel.R`**: the compile branch for `backend = "nutpier"`; the `backend` default changed from `options("brms.backend")`, which returns a list, to `getOption("brms.backend", "rstan")`. `rename_pars()` still runs in the worker, so a worker needs the same brms as the master.

**`R/nutpie_to_stanfit.R`**: `nutpie_draws_to_stanfit(draws, model_name, seed, adapt_delta, max_treedepth, metric)`. Modelled on `brms::read_csv_as_stanfit` and assembled from plain R objects, no un-exported rstan internals. Carries the NUTS sampler parameters (`n_leapfrog__`, `treedepth__`, `divergent__`, …) so `rstan::get_sampler_params()` and `get_divergent_iterations()` work, and `lp__` from `nutpie_diagnostics()$logp` (may differ from CmdStan's by a constant). It does **not** set the per-chain `elapsed_time` attribute, so `rstan::get_elapsed_time()` returns nothing for nutpie fits; the fix is in §6.

**`R/effect_functions.R`** : `effect_from_sim_study()` is the estimand used everywhere, the conditional 12-month difference in differences of the median subject's expected total (`posterior_epred` at `alsfrs_dly_mnths ∈ {0, 12}` × arm, `re_formula = NA`); `effects_lmer()`, `effects_lmer_splines()`, `effects_logit()`; `collect_diags()`, `collect_diags_all()`, `summarise_diags()` with `worst_*` and `pct_*` column names.

**`R/simulate_data.R`**: `visit_id` added; `is_first` and `is_last` exist and the simulator asserts one of each per subject.

## 3. Results to date

### 3.1 Sim1_wide: 100 replicates, p ∈ {0.5, 0.8}, n = 50 per arm, run on server

Models: random-intercept-only `logit`, `logit_TPR`, `logit_TPR_disc`, `logit_disc`; random-slope `logit_rslope`, `logit_rslope_visit`, `logit_rslope_TPR` and their `_disc` twins; LMM baselines `lmer4_splines` and `lmer4_splines_rslope`. Null rejection at p = 0.5 is 1 − coverage:

| tier | rejection at p = 0.5 | mean interval width |
|------------------------|------------------------|------------------------|
| random intercept only (4 models) | 0.29 to 0.31 | 3.6 to 4.4 |
| random slope (3 models and `_disc` twins) | 0.07 to 0.09 | 6.5 to 8.9 |
| LMM control, RI → random slope | 0.26 → 0.06 | 3.8 → 8.5 |

- **The random slope needs to be in the linear predictor.** Omitting slope heterogeneity halves the credible (confidence) interval; the LMM control reproduces the pattern without MCMC, so it is not an ordinal or sampler behaviour.
- **It seems that random-slope models are consistent with nominal at n = 100 but underpowered to prove it:** Settling it needs about 400 to 500 replicates.
- **`logit_rslope_disc` and `logit_rslope_TPR_disc` are computationally invalid** at both p: rhat 2.6 to 3.8, minimum ESS about 3. `disc ~ 0 + alsfrs_dly_mnths` and a random slope on time both make the spread grow with time; jointly unidentified. The non-random-slope `_disc` models converge and show no benefit.
- **Sampling quality was marginal even where valid** under default priors and cmdstan: random-slope rhat 1.13 to 1.17, minimum ESS 14 to 16, 30 to 60 % of fits with divergences, up to about 2 000 in one fit. This is what the nutpie switch and the prior work address.
- **TPR splines** gave narrower intervals (6.8 vs 8.8), better null coverage (0.07 vs 0.09) and lower power at p = 0.8 (0.86 vs 0.94) than the linear-time model. Bias or regularisation is undecided.

### 3.2 s2z (brms PR #1919, sum-to-zero group-level effects): not useful right now

Tested in three configurations against the base parameterisation; ratios are effect ESS per minute relative to base on the same backend: It is correct (equivalence passes everywhere) but useless. Under cmdstan it helps only where the sampler already struggles, under nutpie it is redundant or harmful. What it buys is stability on Q12 under nutpie (rhat 1.022, ESS 137, against base's 1.077 and 28), which the threshold prior now buys more cheaply. It cannot be applied to the multivariate model at all: `mvbind(...) ~ ... (1 + t | p | subject_id)` is refused with "ordinal sum-to-zero group-level ID cannot span multiple linear predictors". More detailed results can be found here: [Google document](https://docs.google.com/document/d/1u-B9efjk9vD1P1Zp4iXq2D2AjKFOMyyJkPp3RLQEklM/edit?usp=sharing).

### 3.3 Prior × backend on the small models (`analyses/prior_backend`)

One replicate (seed 468652233, null, 50 per arm), models × {default, WI} × {cmdstanr, nutpier}, 3 chains, 2 000 sampling draws.

- **Q12 needs the WI block on either backend**: default + cmdstan 12 min with 64 divergences, default + nutpie 1 min but rhat 1.075 and ESS 28 on that threshold, WI + cmdstan clean in 4 min, WI + nutpie clean in 1.2 min. Q5 is fine in all four cells.
- **Twelve items**: the WI block removes cmdstan's divergences (77 → 0) and halves its time (112 → 60 min); nutpie takes 36 min under either prior with 1.5 times the ESS per gradient under WI, and its Q12 threshold failure disappears (rhat 1.034 → 1.002). With WI priors the bottleneck moves from tresholds to the 24 × 24 correlation matrix (ESS 400 to 700 of 6 000).
- **The estimand on the twelve-item total moves 0.31 to 0.34 SD and narrows 20 % under the WI block**, more than on any single item, in the same direction on both backends. The mechanism, as far as the tables show: the twelve group main effects shrink 25 to 38 % under N(0, 3), because they share the treated arm's baseline severity and the twelve priors act jointly on it, and the DID on the points scale depends on the group effects through the curvature of the expected-score function. The interaction SDs are unchanged while the group-effect SDs narrow 11 %..
- Results can be found here: [Google document](https://docs.google.com/document/d/1u-B9efjk9vD1P1Zp4iXq2D2AjKFOMyyJkPp3RLQEklM/edit?usp=sharing).

## 4. Decisions and their evidence

| decision | evidence |
|------------------------------------|------------------------------------|
| Random slope on time is required in the ordinal models | Sim1_wide tiers, LMM control (§4.1) |
| `logit_rslope_disc` and `logit_rslope_TPR_disc` are discarded | rhat 2.6 to 3.8, ESS ≈ 3 (§4.1) |
| nutpie replaces cmdstan as the sampler in `brm_parallel` | measured 3.66 times the effect ESS per minute on long model `nutpie_stan_benchmark.R`, 2 to 10× on wide models `prior_backend_small.R` |
| Diagonal adaptation, nutpie's default 400 warm-up | low-rank worse; results in `nutpie_stan_benchmark.R` |
| s2z is closed | no benefit under nutpie, incompatible with `\| p \|` (§3.2) |
| The WI prior block replaces default priors | `prior_backend_small.R:`no divergences, half the cmdstan cost, threshold non-convergence removed |

## 5. Caveats and bugs

- **The nutpie sampler clock.** `rstan::get_elapsed_time()` reads an `elapsed_time` attribute the converter never sets. Fix: `system.time()` around `nutpie_sample()` in `sampling_parallel.R`, pass `elapsed_sec` to the converter (guard with `"elapsed_sec" %in% names(formals(nutpie_converter))`), and set `attr(samples[[i]], "elapsed_time") <- c(warmup = 0, sample = elapsed_sec)` per chain. The number then includes nutpie's 400 warm-up draws; ESS per gradient is the clean measure.

- **cmdstanr translation**: `chains → num_chains` (deprecated) and `control$max_treedepth → max_depth` (wrong name). No call site passes `control`.

- **Cache keys** exclude sampler arguments (§2.2).

- **The TPR construction** in Sim1 (`s(t, by = group)`, unordered factor) has no treatment parameter and centres each arm's curve on its own rows; with a real effect the treated arm's baseline offset leaks into the random slopes, a candidate explanation for TPR's narrower intervals and lower power at p = 0.8.

- **`min_ess_bulk` in Sim1's table** is dominated by Q12's near-unidentified thresholds and says nothing about the effect.

- **Conditional versus marginal estimand.** `effect_from_sim_study` uses `re_formula = NA`, the median subject (u = 0). The LMM's 12·γ₃ is both conditional and marginal because random effects cancel in a linear model; in the ordinal model they do not. Measured on the Sim1 fits: conditional ≈ 1.7 to 1.9 × marginal (6.81 vs 3.98 points; −7.55 vs −3.96). Under p = 0.5 both are zero, so null calibration is unaffected, but power at p = 0.8 is not comparable across model families and any bias against the truth is on the wrong scale. Marginal version: `posterior_epred(re_formula = NULL, allow_new_levels = TRUE, sample_new_levels = "gaussian")` with about 200 new subject ids in `newdata`, averaging within each draw; code in `effect_hier_from_sim_study` in `Sim1_long.Rmd`. Averaging over the fitted subjects instead is invariant to where the thresholds sit relative to the random-intercept mean, which the threshold prior otherwise shifts by about a logit.

- **Dependence of the estimand on the group main effect.** The DID on the expected-score scale reads the treated arm's 12-month change from its own baseline position on the S-shaped response curve, so a chance baseline imbalance between arms enters it through curvature even when the interactions are zero, and any prior on the group effect moves it. On the latent scale the DID is exactly 12 Σ β_j3 and does not see the group effect.

## 6. Recommendations

- **Models to keep.** The wide model with a random slope, in two variants: linear time (`logit_rslope`) and spline time. The second formulation worth keeping for its logical structure is the long-format hierarchical model in `Sim1_long.Rmd`, whose single `time_treat` parameter is the right geometry for a systemic treatment and the natural place for a prior calibrated on the total. Drop the `_disc` random-slope variants and, unless a specific question needs them, the random-intercept-only models.
- **Sampler.** nutpie, diagonal adaptation, through `brm_parallel(backend = "nutpier")`. It is equivalent and 1.7 to 10 times faster.
- **Priors.** The part that matters mostly for sampling is the N(0, 10) on thresholds; it is what removes the divergences and the heavy tail on unobserved-category thresholds. The line that matters for the estimand is the group main effect in wide models, which I would remove from the model rather than keep with any prior.
- **The group main effect.** The simulator matches pairs on baseline total before assigning arms, so the arm difference at month 0 is zero by design. Fitting a parameter whose truth is known to be zero costs efficiency and, in the ordinal model, contaminates the change contrast. Exclude it from the linear predictor, especially when exploring priors. The LMM baselines can keep theirs, its estimand does not see it, or be constrained too, which makes them more efficient.
- **Marginal versus conditional.** For comparability with the LMM, compute the marginal estimand, preferably by averaging over the fitted subjects, which also removes the dependence on the threshold location. Keep the conditional one alongside until the two have been compared on a replicate.
- **Unobserved categories.** Write a function that checks each generated dataset item by item and merges a category that is never observed with the one above it, so that no threshold is left identified by the prior alone. Make the rule fixed per item from the registry frequencies rather than decided per replicate, otherwise the model changes between replicates.
- **Replicates.** Increase `nsims` in `prior_backend_small.R` to 5 or 10 for Q5, Q12 and the twelve-item model on nutpie (about 24 minutes per replicate for the single items, 70 for the twelve-item model with both priors); the report script summarises paired contrasts over replicates. For the calibration study, 400 to 500 replicates on the server.
- **Then rerun Sim1_wide and Sim1_long** with the "to keep" models only, with nutpie, the WI block, no group main effect and the marginal estimand, at both eff \_probs, on the server, with a fresh `cache_dir`.
