# Sampler benchmark on the `hier` (long) ordinal model.
#
# Three arms on the SAME model, SAME data, SAME number of post-warmup draws:
#   1. cmdstanr, diag_e metric, 1000 warmup   (the current production sampler)
#   2. nutpieR,  diag adaptation, 400 warmup  (nuts-rs default)
#   3. nutpieR,  low_rank adaptation, 800 warmup
# All three go through brm_parallel(backend = ...), so this also exercises the
# pipeline integration at real scale, not just the sampler.
#
# Decision metric is ess_gamma_per_min: effective draws of the treatment-effect
# coefficient b_time_treat per wall minute.

# ---- 1. SETUP -------------------------------------------------------------
library(dplyr)
library(tidyr)
library(brms)
library(readr)
library(posterior)
library(nutpieR)
source(here::here("R", "data_prep.R"))
source(here::here("R", "simulate_data.R"))
source(here::here("R", "nutpie_to_stanfit.R"))
source(here::here("R", "sampling_parallel.R"))
source(here::here("R", "brm_parallel.R"))

versions <- c(nutpieR = as.character(packageVersion("nutpieR")),
              brms = as.character(packageVersion("brms")),
              rstan = as.character(packageVersion("rstan")),
              posterior = as.character(packageVersion("posterior")),
              cmdstan = as.character(cmdstanr::cmdstan_version()))
print(versions)

# one fit at a time, three chains inside it: no contention between arms
future::plan(future::multisession, workers = 1)
cores_per_fit <- 3

results_file <- here::here("local_temp_data", "nutpie_benchmark_results.rds")
results <- list(versions = versions, started = Sys.time(), arms = list())
checkpoint <- function() saveRDS(results, results_file)

# ---- 2. DATA --------------------------------------------------------------
sim_seed <- 468652233
sampler_seed <- 1
n_subj_per_group <- 50
question_col_names <- sprintf("Q%02d", 1:12)

proact <- read_csv(here::here("private_data/F_PROACT_ALSFRS.csv"),
                   show_col_types = FALSE)
proact_standardized <- data_prep(proact)

set.seed(sim_seed)
sim <- simulate_data_from_registry(proact_standardized,
                                   max_duration = 12,
                                   min_measurements_per_subject = 4,
                                   max_measurements_per_subject = 8,
                                   n_subjects_per_group = n_subj_per_group,
                                   effect_prob = 0.5)   # exact strong null

d_long <- sim |>
    select(subject_id, visit_id, group, alsfrs_dly_mnths, all_of(question_col_names)) |>
    pivot_longer(all_of(question_col_names), names_to = "item", values_to = "score") |>
    mutate(score = score + 1,
           item = factor(item, levels = question_col_names),
           treat = as.integer(group == "Treatment"),
           time = alsfrs_dly_mnths,
           time_treat = time * treat)

results$sim_seed <- sim_seed
results$sampler_seed <- sampler_seed
results$n_obs <- nrow(d_long)
results$n_subjects <- dplyr::n_distinct(d_long$subject_id)
checkpoint()

# ---- 3. MODEL -------------------------------------------------------------
hier_formula <- bf(score | thres(4, gr = item) ~ item:time + time_treat +
                       (1 + time | subject_id) + (1 + time | subject_id:item) +
                       (1 | subject_id:visit_id),
                   family = cumulative("logit"))

hier_priors <- c(
    prior(normal(0, 2), class = "Intercept"),
    prior(normal(0, 0.5), class = "b"),
    prior(normal(0, 0.25), class = "b", coef = "time_treat"),
    prior(student_t(3, 0, 2.5), class = "sd"),
    prior(normal(0, 0.5), class = "sd", coef = "time", group = "subject_id"),
    prior(normal(0, 0.5), class = "sd", coef = "time", group = "subject_id:item"),
    prior(normal(0, 1), class = "sd", group = "subject_id:visit_id"),
    prior(lkj(2), class = "cor")
)

# iter/warmup are the brms-style arguments Sim1.Rmd already uses. Under cmdstanr
# they become iter_warmup = 1000 / iter_sampling = 1000; under nutpie they become
# num_draws = 1000 with num_warmup left at the nuts-rs default (400 diag, 800
# low_rank). Same post-warmup draws, different adaptation budget - that
# difference is the thing being measured.
args_sh <- list(formula = hier_formula, prior = hier_priors,
                adapt_delta = 0.95, init = 0.1, seed = sampler_seed,
                chains = 3, iter = 2000, warmup = 1000)
args_per_fit <- list(list(data = d_long))

# ---- 4. PRE-COMPILE (outside the timed region) -----------------------------
# Both backends cache compiled artefacts by source content, so the compile
# inside the timed brm_parallel() call becomes a near-instant cache hit and the
# wall times compare sampling, not compilation.
stancode <- as.character(make_stancode(hier_formula, data = d_long, prior = hier_priors))
cat("pre-compiling cmdstan ...\n")
invisible(cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode)))
cat("pre-compiling nutpie/BridgeStan (slow on first ever use) ...\n")
invisible(nutpieR::nutpie_compile_model(code = stancode, verbose = 0))

# ---- 5. WHAT EACH ARM REPORTS ---------------------------------------------
# Conditional (re_formula = NA) DID of the expected total score at 12 months --
# the median-patient estimand, matching effect_from_sim_study.
effect_hier_cond <- function(bfit, s_star = 12) {
    items <- levels(bfit$data$item)
    nd <- expand.grid(item = factor(items, levels = items),
                      time = c(0, s_star), treat = c(0, 1), KEEP.OUT.ATTRS = FALSE)
    nd$time_treat <- nd$time * nd$treat
    pred <- brms::posterior_epred(bfit, newdata = nd, re_formula = NA)
    cats <- as.integer(dimnames(pred)[[3]])
    esc <- 0
    for (k in seq_along(cats)) esc <- esc + pred[, , k] * cats[k]
    keys <- paste(nd$time, nd$treat, sep = "_")
    tot <- sapply(unique(keys), function(kk) rowSums(esc[, keys == kk, drop = FALSE]))
    (tot[, paste0(s_star, "_1")] - tot[, "0_1"]) -
        (tot[, paste0(s_star, "_0")] - tot[, "0_0"])
}

# Runs in the worker and returns a few MB, not a multi-GB brmsfit.
bench_summary <- function(brmsfit) {
    dr <- posterior::as_draws(brmsfit)
    eff <- effect_hier_cond(brmsfit)
    sp <- rstan::get_sampler_params(brmsfit$fit, inc_warmup = FALSE)
    list(vars = posterior::variables(dr),
         n_draws = posterior::ndraws(dr),
         summ = posterior::summarise_draws(dr, "mean", "sd", "ess_bulk", "rhat"),
         effect = as.numeric(eff),
         ess_effect = posterior::ess_bulk(matrix(eff, ncol = posterior::nchains(dr))),
         divergences = sum(rstan::get_divergent_iterations(brmsfit$fit)),
         max_treedepth_hits = sum(vapply(sp, function(x) sum(x[, "treedepth__"] >= 10),
                                         numeric(1))),
         metric = brmsfit$fit@stan_args[[1]]$metric)
}

run_arm <- function(label, backend, nutpie_args = list()) {
    cat("\n=== arm:", label, "=== started", format(Sys.time(), "%H:%M:%S"), "\n")
    t <- system.time(
        res <- brm_parallel(args_shared = args_sh, args_per_fit = args_per_fit,
                            backend = backend,
                            cores_per_fit = cores_per_fit,
                            cache_fits = FALSE, cache_summaries = FALSE,
                            summarise_fun = bench_summary,
                            future.globals = list(
                                "bench_summary" = bench_summary,
                                "effect_hier_cond" = effect_hier_cond),
                            nutpie_args = nutpie_args)
    )
    out <- res[[1]]
    out$wall_sec <- unname(t["elapsed"])
    out$label <- label
    out$backend <- backend
    out$nutpie_args <- nutpie_args
    cat(sprintf("    %.1f min | %d divergences | metric %s\n",
                out$wall_sec / 60, out$divergences, out$metric))
    out
}

# ---- 6. THE THREE ARMS ----------------------------------------------------
results$arms$cmdstan_diag <- run_arm("cmdstan diag", "cmdstanr")
checkpoint()

results$arms$nutpie_diag <- run_arm("nutpie diag", "nutpier")
checkpoint()

results$arms$nutpie_lowrank <- run_arm("nutpie low_rank", "nutpier",
                                       nutpie_args = list(adaptation = "low_rank"))
checkpoint()

# ---- 7. VALIDATION --------------------------------------------------------
# Nothing below decides anything until these pass: a faster sampler that
# disagrees with the current one is not a faster sampler, it is a bug.
arms <- results$arms
ref <- arms$cmdstan_diag

cat("\n--------------- VALIDATION vs cmdstan diag ---------------\n")
for (nm in setdiff(names(arms), "cmdstan_diag")) {
    a <- arms[[nm]]
    cat("\n##", a$label, "\n")

    cat("  variables identical:", identical(sort(ref$vars), sort(a$vars)),
        sprintf("(%d vs %d)\n", length(ref$vars), length(a$vars)))
    cat("  post-warmup draws equal:", identical(ref$n_draws, a$n_draws),
        sprintf("(%d vs %d)\n", ref$n_draws, a$n_draws))

    # lp__ is excluded: nutpie's log density can differ from CmdStan's by an
    # additive constant, so it is not comparable across samplers.
    cmp <- dplyr::inner_join(a$summ, ref$summ, by = "variable",
                             suffix = c("_a", "_ref")) |>
        dplyr::filter(sd_a > 0 | sd_ref > 0, variable != "lp__") |>
        dplyr::mutate(z = (mean_a - mean_ref) /
                          sqrt(sd_a^2 / ess_bulk_a + sd_ref^2 / ess_bulk_ref)) |>
        dplyr::filter(!is.na(z)) |>
        dplyr::arrange(dplyr::desc(abs(z)))
    # with this many variables, max |z| ~ sqrt(2 log N) ~ 4 is expected
    cat(sprintf("  z agreement over %d variables: max |z| = %.2f, frac |z| > 3 = %.4f\n",
                nrow(cmp), max(abs(cmp$z)), mean(abs(cmp$z) > 3)))
    cat("  top offenders (healthy = scattered random effects, both signs):\n")
    print(head(cmp[, c("variable", "mean_a", "mean_ref", "z")], 5))

    cat(sprintf("  effect (cond. DID @12m): %s %.2f [%.2f, %.2f] vs cmdstan %.2f [%.2f, %.2f]\n",
                a$label, mean(a$effect), quantile(a$effect, .025), quantile(a$effect, .975),
                mean(ref$effect), quantile(ref$effect, .025), quantile(ref$effect, .975)))
    cat(sprintf("  max rhat: %.3f (cmdstan %.3f)\n",
                max(a$summ$rhat, na.rm = TRUE), max(ref$summ$rhat, na.rm = TRUE)))
}

# ---- 8. DECISION TABLE ----------------------------------------------------
gamma_row <- function(a) {
    g <- a$summ[a$summ$variable == "b_time_treat", ]
    if (nrow(g) != 1) {
        stop("b_time_treat not found in ", a$label,
             " -- the decision metric needs it; check the formula")
    }
    g
}

table_bench <- do.call(rbind, lapply(names(arms), function(nm) {
    a <- arms[[nm]]
    g <- gamma_row(a)
    mins <- a$wall_sec / 60
    data.frame(arm = a$label,
               wall_min = round(mins, 1),
               ess_gamma = round(g$ess_bulk),
               ess_gamma_per_min = round(g$ess_bulk / mins, 1),
               ess_effect = round(a$ess_effect),
               ess_effect_per_min = round(a$ess_effect / mins, 1),
               min_ess_bulk = round(min(a$summ$ess_bulk, na.rm = TRUE)),
               max_rhat = round(max(a$summ$rhat, na.rm = TRUE), 3),
               divergences = a$divergences,
               treedepth_hits = a$max_treedepth_hits,
               stringsAsFactors = FALSE)
}))

ref_rate <- table_bench$ess_gamma_per_min[table_bench$arm == "cmdstan diag"]
table_bench$ratio_vs_cmdstan <- round(table_bench$ess_gamma_per_min / ref_rate, 2)

cat("\n--------------- BENCHMARK ---------------\n")
print(table_bench, row.names = FALSE)

cat("\npre-registered rule: >= 2.00 switch | < 1.50 stay | between = judgement\n")
for (i in seq_len(nrow(table_bench))) {
    r <- table_bench$ratio_vs_cmdstan[i]
    if (table_bench$arm[i] == "cmdstan diag") next
    cat(sprintf("  %-16s ratio %.2f -> %s\n", table_bench$arm[i], r,
                if (r >= 2) "SWITCH" else if (r < 1.5) "STAY" else "judgement call"))
}

results$table <- table_bench
results$finished <- Sys.time()
checkpoint()
future::plan(future::sequential)
cat("\nresults saved to", results_file, "\n")
