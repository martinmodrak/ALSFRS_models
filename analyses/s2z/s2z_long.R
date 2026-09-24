# This file answers three questions, in this order. Question 3 does not matter if
# question 2 fails: a faster sampler that answers a different question is worthless.
#
#   1. CONVERTS?   nutpie -> stanfit -> rename_pars -> posterior_epred, now with
#                  177 lines of generated quantities reconstructing b, the
#                  thresholds and every r_*. This is the only layer of nutpie/s2z
#                  compatibility that could not be checked statically.
#   2. SAME POSTERIOR?  s2z claims to be an EXACT reparameterisation for Gaussian
#                  group effects. On the same backend, b_time_treat and the 12-month
#                  DID effect must agree with the baseline within Monte Carlo error.
#   3. FASTER?     ess per wall minute, bulk AND tail. Tail is reported because the
#                  simulation study lives on interval coverage, not posterior means.
#
# Four fits at reduced settings: three on cmdstanr (equivalence + speed, sampler
# held constant) and one on nutpie (conversion + cross-backend agreement).
#

S2Z_LIB <- path.expand("~/R/brms-s2z")
if ("brms" %in% loadedNamespaces()) {
    stop("brms is ALREADY LOADED (from ", dirname(system.file(package = "brms")), ").\n",
         "  Switching .libPaths() now has no effect. Run this from a terminal:\n",
         "     Rscript analyses/s2z_stage_c.R")
}
.libPaths(c(S2Z_LIB, .libPaths()))

suppressMessages({
    library(dplyr); library(tidyr); library(brms); library(readr); library(posterior)
})
library(nutpieR)
source(here::here("R", "data_prep.R"))
source(here::here("R", "simulate_data.R"))
source(here::here("R", "nutpie_to_stanfit.R"))
source(here::here("R", "sampling_parallel.R"))
source(here::here("R", "brm_parallel.R"))

stopifnot("s2z" %in% names(formals(brms::gr)))
cat("brms:", as.character(packageVersion("brms")), "from", dirname(system.file(package = "brms")),
    "| commit", substr(packageDescription("brms")$RemoteSha, 1, 8), "\n")

# The worker must see the SAME brms. parallelly does not always propagate a
# .libPaths() set at runtime, and a worker running brms 2.23.1 would call
# rename_pars() on an s2z fit -- wrong results, no error.
cl <- parallelly::makeClusterPSOCK(1L, rscript_libs = .libPaths())
future::plan(future::cluster, workers = cl)
worker_lib <- future::value(future::future(dirname(system.file(package = "brms"))))
cat("worker brms lib:", worker_lib, "\n")
if (!identical(normalizePath(worker_lib), normalizePath(S2Z_LIB))) {
    stop("The future worker does not see the s2z brms (", worker_lib, "). Aborting: ",
         "rename_pars() would run under the wrong brms.")
}
cores_per_fit <- 3

# ---- data: real simulator, reduced n so Stage C is minutes not hours --------
sim_seed <- 468652233
sampler_seed <- 1
n_subj_per_group <- 20
question_col_names <- sprintf("Q%02d", 1:12)

proact <- read_csv(here::here("private_data/F_PROACT_ALSFRS.csv"), show_col_types = FALSE)
set.seed(sim_seed)
sim <- simulate_data_from_registry(data_prep(proact),
                                   max_duration = 12,
                                   min_measurements_per_subject = 4,
                                   max_measurements_per_subject = 8,
                                   n_subjects_per_group = n_subj_per_group,
                                   effect_prob = 0.5)

d_long <- sim |>
    select(subject_id, visit_id, group, alsfrs_dly_mnths, all_of(question_col_names)) |>
    pivot_longer(all_of(question_col_names), names_to = "item", values_to = "score") |>
    mutate(score = score + 1,
           item = factor(item, levels = question_col_names),
           treat = as.integer(group == "Treatment"),
           time = alsfrs_dly_mnths,
           time_treat = time * treat)

# ---- the three parameterisations -------------------------------------------
mk <- function(rhs) bf(as.formula(paste("score | thres(4, gr = item) ~", rhs)),
                       family = cumulative("logit"))

f_base <- mk("item:time + time_treat + (1 + time | subject_id) +
              (1 + time | subject_id:item) + (1 | subject_id:visit_id)")

# center = FALSE: non-centered s2z coordinates. NOTE the inverted default --
# omitting `center` on an s2z block gives CENTERED coordinates, which is the
# configuration reported as a large slowdown on the ordinal verbagg example.
f_s2z_nc <- mk("item:time + time_treat +
              (1 + time | gr(subject_id, s2z = TRUE, center = FALSE)) +
              (1 + time | gr(subject_id:item, s2z = TRUE, center = FALSE)) +
              (1 | gr(subject_id:visit_id, s2z = TRUE, center = FALSE))")

f_s2z_c <- mk("item:time + time_treat +
              (1 + time | gr(subject_id, s2z = TRUE)) +
              (1 + time | gr(subject_id:item, s2z = TRUE)) +
              (1 | gr(subject_id:visit_id, s2z = TRUE))")

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

sampler_args <- list(adapt_delta = 0.95, init = 0.1, seed = sampler_seed,
                     chains = 2, iter = 1000, warmup = 500)
args_per_fit <- list(list(data = d_long))

results_file <- here::here("local_temp_data", "s2z_stage_c_results.rds")
results <- list(started = Sys.time(), n_subj_per_group = n_subj_per_group,
                brms_sha = packageDescription("brms")$RemoteSha,
                sampler_args = sampler_args, arms = list())
checkpoint <- function() saveRDS(results, results_file)

# ---- pre-compile outside the timed region -----------------------------------
cat("\npre-compiling (cmdstan + nutpie) ...\n")
for (nm in c("base", "s2z_nc", "s2z_c")) {
    f <- switch(nm, base = f_base, s2z_nc = f_s2z_nc, s2z_c = f_s2z_c)
    sc <- as.character(make_stancode(f, data = d_long, prior = hier_priors))
    invisible(cmdstanr::cmdstan_model(cmdstanr::write_stan_file(sc)))
    if (nm == "s2z_nc") invisible(nutpieR::nutpie_compile_model(code = sc, verbose = 0))
    cat("  ", nm, "compiled\n")
}

# ---- what each fit reports --------------------------------------------------
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

stage_c_summary <- function(brmsfit) {
    dr <- posterior::as_draws(brmsfit)
    eff <- effect_hier_cond(brmsfit)
    em <- matrix(eff, ncol = posterior::nchains(dr))
    list(vars = posterior::variables(dr),
         n_draws = posterior::ndraws(dr),
         summ = posterior::summarise_draws(dr, "mean", "sd", "ess_bulk", "ess_tail", "rhat"),
         effect = as.numeric(eff),
         ess_effect_bulk = posterior::ess_bulk(em),
         ess_effect_tail = posterior::ess_tail(em),
         divergences = sum(rstan::get_divergent_iterations(brmsfit$fit)),
         metric = brmsfit$fit@stan_args[[1]]$metric)
}

run_arm <- function(label, formula, backend) {
    cat("\n=== ", label, " === started", format(Sys.time(), "%H:%M:%S"), "\n")
    t <- system.time(
        res <- brm_parallel(
            args_shared = c(list(formula = formula, prior = hier_priors), sampler_args),
            args_per_fit = args_per_fit,
            backend = backend,
            cores_per_fit = cores_per_fit,
            cache_fits = FALSE, cache_summaries = FALSE,
            summarise_fun = stage_c_summary,
            future.globals = list("stage_c_summary" = stage_c_summary,
                                  "effect_hier_cond" = effect_hier_cond))
    )
    out <- res[[1]]
    out$wall_sec <- unname(t["elapsed"]); out$label <- label; out$backend <- backend
    cat(sprintf("    %.1f min | %d vars | %d divergences | metric %s\n",
                out$wall_sec / 60, length(out$vars), out$divergences, out$metric))
    out
}

# ---- the four fits ----------------------------------------------------------
results$arms$base_cs   <- run_arm("base / cmdstanr",            f_base,   "cmdstanr"); checkpoint()
results$arms$s2z_nc_cs <- run_arm("s2z center=FALSE / cmdstanr", f_s2z_nc, "cmdstanr"); checkpoint()
results$arms$s2z_c_cs  <- run_arm("s2z centered / cmdstanr",     f_s2z_c,  "cmdstanr"); checkpoint()
results$arms$s2z_nc_np <- run_arm("s2z center=FALSE / nutpie",   f_s2z_nc, "nutpier");  checkpoint()

a <- results$arms

# ---- QUESTION 1: does the nutpie + s2z fit convert? -------------------------
cat("\n\n############ 1. CONVERSION (nutpie + s2z) ############\n")
np <- a$s2z_nc_np
cat("  brmsfit built and summarised     :", !is.null(np$summ), "\n")
cat("  posterior_epred ran (effect draws):", length(np$effect), "\n")
cat("  get_divergent_iterations worked   :", is.numeric(np$divergences),
    "| divergences:", np$divergences, "\n")
cat("  variables                         :", length(np$vars),
    "(cmdstanr s2z arm has", length(a$s2z_nc_cs$vars), ")\n")
cat("  reconstructed GQ present          :",
    all(c("b_time_treat", "sd_subject_id__Intercept") %in% np$vars), "\n")
cat("  metric label                      :", np$metric, "\n")

# ---- QUESTION 2: is it the same posterior? ---------------------------------
cat("\n############ 2. EQUIVALENCE vs base (same backend) ############\n")
ref <- a$base_cs
# Three quantities are computed in the model's OWN coordinates and are therefore
# meaningless to compare across parameterisations:
#   lp__             - log posterior, differs by the Jacobian/normalising constant
#   lprior           - log prior, evaluated on theta_*/z_s2z_* rather than b/z
#   merged_Intercept - brms's INTERNAL threshold vector, which under s2z holds the
#                      shifted coordinates, not the reconstructed thresholds. The
#                      public, reconstructed thresholds are b_Intercept[...], and
#                      those DO agree (max |z| = 2.6 in the first Stage C run).
PARAM_DEPENDENT <- c("lp__", "lprior")
zcmp <- function(x, y) {
    dplyr::inner_join(x$summ, y$summ, by = "variable", suffix = c("_a", "_ref")) |>
        dplyr::filter(sd_a > 0 | sd_ref > 0,
                      !variable %in% PARAM_DEPENDENT,
                      !grepl("^merged_Intercept", variable)) |>
        dplyr::mutate(z = (mean_a - mean_ref) /
                          sqrt(sd_a^2 / ess_bulk_a + sd_ref^2 / ess_bulk_ref)) |>
        dplyr::filter(!is.na(z)) |> dplyr::arrange(dplyr::desc(abs(z)))
}
eff_z <- function(x, y) {
    (mean(x$effect) - mean(y$effect)) /
        sqrt(sd(x$effect)^2 / x$ess_effect_bulk + sd(y$effect)^2 / y$ess_effect_bulk)
}
for (nm in c("s2z_nc_cs", "s2z_c_cs", "s2z_nc_np")) {
    x <- a[[nm]]
    cmp <- zcmp(x, ref)
    g <- cmp[cmp$variable == "b_time_treat", ]
    cat("\n##", x$label, "\n")
    cat(sprintf("   shared variables: %d of %d (s2z adds theta_*/fixed_s2z/q_recovered_*)\n",
                nrow(cmp), length(x$vars)))
    cat(sprintf("   max |z| = %.2f | frac |z| > 3 = %.4f   (expect max ~4 by chance at this size)\n",
                max(abs(cmp$z)), mean(abs(cmp$z) > 3)))
    cat(sprintf("   b_time_treat: %.4f (sd %.4f) vs base %.4f (sd %.4f) -> z = %.2f\n",
                g$mean_a, g$sd_a, g$mean_ref, g$sd_ref, g$z))
    cat(sprintf("   effect @12m : %.3f [%.3f, %.3f] vs base %.3f [%.3f, %.3f] -> z = %.2f\n",
                mean(x$effect), quantile(x$effect, .025), quantile(x$effect, .975),
                mean(ref$effect), quantile(ref$effect, .025), quantile(ref$effect, .975),
                eff_z(x, ref)))
    cat("   worst-agreeing variables:\n")
    print(head(cmp[, c("variable", "mean_a", "mean_ref", "z")], 3))
}

# ---- QUESTION 3: is it faster? ---------------------------------------------
cat("\n############ 3. EFFICIENCY ############\n")
gm <- function(x, col) {
    r <- x$summ[x$summ$variable == "b_time_treat", ]
    if (nrow(r) != 1) NA_real_ else r[[col]]
}
tab <- do.call(rbind, lapply(names(a), function(nm) {
    x <- a[[nm]]; mins <- x$wall_sec / 60
    data.frame(arm = x$label,
               wall_min = round(mins, 1),
               ess_gamma_bulk = round(gm(x, "ess_bulk")),
               ess_gamma_tail = round(gm(x, "ess_tail")),
               bulk_per_min = round(gm(x, "ess_bulk") / mins, 1),
               tail_per_min = round(gm(x, "ess_tail") / mins, 1),
               eff_bulk_per_min = round(x$ess_effect_bulk / mins, 1),
               eff_tail_per_min = round(x$ess_effect_tail / mins, 1),
               min_ess_bulk = round(min(x$summ$ess_bulk, na.rm = TRUE)),
               max_rhat = round(max(x$summ$rhat, na.rm = TRUE), 3),
               divergences = x$divergences,
               stringsAsFactors = FALSE)
}))
base_rate <- tab$bulk_per_min[tab$arm == "base / cmdstanr"]
tab$ratio_vs_base <- round(tab$bulk_per_min / base_rate, 2)
print(tab, row.names = FALSE)

cat("\nRead question 2 first. If any s2z arm disagrees with base on b_time_treat or\n",
    "the effect (|z| > 3 with no MC explanation), s2z is NOT reproducing your model\n",
    "and the timings below it are irrelevant. Only if equivalence holds is a ratio\n",
    ">= 2 worth taking to Stage D.\n")

results$table <- tab
results$finished <- Sys.time()
checkpoint()
parallel::stopCluster(cl); future::plan(future::sequential)
cat("\nsaved to", results_file, "\n")
