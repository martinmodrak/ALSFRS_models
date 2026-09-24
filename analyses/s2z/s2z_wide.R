# WIDE format, ONE question at a time -- Q05 and Q12.
#
#   Qxx | thres(4) ~ 1 + alsfrs_dly_mnths * group + (1 + alsfrs_dly_mnths | subject_id)
#
# Why these two items:
#   Q05  steepest decline in PROACT (-0.119 pts/month), well spread over categories
#   Q12  ceiling item, 84% in the top category (94.5% at baseline): thresholds
#        barely identified, which is where a location ridge bites hardest --
#        and, under Sim1.Rmd's flat prior on b, where sampling is slowest
#


S2Z_LIB <- path.expand("~/R/brms-s2z")
if ("brms" %in% loadedNamespaces()) {
    stop("brms is ALREADY LOADED (from ", dirname(system.file(package = "brms")), ").\n",
         "  Run this from a terminal:  Rscript analyses/s2z_stage_c_wide.R")
}
.libPaths(c(S2Z_LIB, .libPaths()))

suppressMessages({
    library(dplyr); library(brms); library(readr); library(posterior)
})
library(nutpieR)
source(here::here("R", "data_prep.R"))
source(here::here("R", "simulate_data.R"))
source(here::here("R", "effect_functions.R"))
source(here::here("R", "nutpie_to_stanfit.R"))
source(here::here("R", "sampling_parallel.R"))
source(here::here("R", "brm_parallel.R"))

stopifnot("s2z" %in% names(formals(brms::gr)))
cat("brms:", as.character(packageVersion("brms")), "| commit",
    substr(packageDescription("brms")$RemoteSha, 1, 8), "\n")

cl <- parallelly::makeClusterPSOCK(1L, rscript_libs = .libPaths())
future::plan(future::cluster, workers = cl)
if (!identical(normalizePath(future::value(future::future(dirname(system.file(package = "brms"))))),
               normalizePath(S2Z_LIB))) {
    stop("The future worker does not see the s2z brms.")
}
cat("worker brms lib: ok\n")
cores_per_fit <- 3

QUESTIONS <- c("Q05", "Q12")

# ---- data: WIDE, exactly as Sim1.Rmd builds it, at production n -------------
sim_seed <- 468652233
sampler_seed <- 1
n_subj_per_group <- 50
question_col_names <- sprintf("Q%02d", 1:12)

proact <- read_csv(here::here("private_data/F_PROACT_ALSFRS.csv"), show_col_types = FALSE)
set.seed(sim_seed)
d_wide <- simulate_data_from_registry(data_prep(proact),
                                      max_duration = 12,
                                      min_measurements_per_subject = 4,
                                      max_measurements_per_subject = 8,
                                      n_subjects_per_group = n_subj_per_group,
                                      effect_prob = 0.5) |>            # exact strong null
    mutate(across(all_of(question_col_names), ~ .x + 1))                # brms wants 1..K

cat("\nrows:", nrow(d_wide), "| subjects:", n_distinct(d_wide$subject_id), "\n")
for (q in QUESTIONS) {
    cat(sprintf("  %s category counts: %s\n", q,
                paste(sprintf("%d:%d", 1:5, tabulate(d_wide[[q]], 5)), collapse = "  ")))
}

# ---- the three parameterisations, per item ---------------------------------
GR <- list(
    base   = "subject_id",
    s2z_nc = "gr(subject_id, s2z = TRUE, center = FALSE)",
    s2z_c  = "gr(subject_id, s2z = TRUE)")     # centered: the PR's DEFAULT for s2z

mk <- function(q, gr_spec) {
    bf(as.formula(sprintf(
        "%s | thres(4) ~ 1 + alsfrs_dly_mnths * group + (1 + alsfrs_dly_mnths | %s)",
        q, gr_spec)), family = cumulative())
}

# Sim1.Rmd passes no `prior`, so class b gets brms's flat default. The PR states
# s2z supports normal and student-t priors, so test rather than assume.
prior_sets <- list(
    brms_default = NULL,
    weakly_informative = c(prior(normal(0, 2), class = "Intercept"),
                           prior(normal(0, 1), class = "b"),
                           prior(student_t(3, 0, 2.5), class = "sd"),
                           prior(lkj(2), class = "cor")))

cat("\n=========== STAGE A: which prior set works for every item x parameterisation? ===========\n")
grid <- expand.grid(item = QUESTIONS, par = names(GR), pri = names(prior_sets),
                    stringsAsFactors = FALSE)
grid$ok <- NA; grid$err <- ""
for (i in seq_len(nrow(grid))) {
    sc <- tryCatch(make_stancode(mk(grid$item[i], GR[[grid$par[i]]]), data = d_wide,
                                 prior = prior_sets[[grid$pri[i]]]),
                   error = function(e) e)
    grid$ok[i] <- !inherits(sc, "error")
    if (!grid$ok[i]) grid$err[i] <- substr(conditionMessage(sc), 1, 200)
}
print(grid[, c("item", "par", "pri", "ok")], row.names = FALSE)
for (i in which(!grid$ok)) cat("  !", grid$item[i], grid$par[i], grid$pri[i], ":", grid$err[i], "\n")

usable <- names(prior_sets)[vapply(names(prior_sets),
                                   function(p) all(grid$ok[grid$pri == p]), logical(1))]
if (!length(usable)) stop("No prior set works for all parameterisations; see errors above.")
PRIOR_SET <- usable[1]
cat("\nusing prior set:", PRIOR_SET,
    if (PRIOR_SET != "brms_default") " (NOTE: NOT what Sim1.Rmd uses)" else " (as in Sim1.Rmd)", "\n")
the_prior <- prior_sets[[PRIOR_SET]]

sd0 <- make_standata(mk(QUESTIONS[1], GR$base), data = d_wide)
cat("sampled dimensions (base, one item):",
    sd0$nthres + sd0$Kc + sd0$M_1 * sd0$N_1 + sd0$M_1 + sd0$M_1 * (sd0$M_1 - 1) / 2, "\n")

# ---- what each fit reports --------------------------------------------------
# effect_from_sim_study() is YOUR function, unchanged: on a univariate ordinal
# model pred[,1,] is already [draws x categories], so it returns the DID of this
# one item's expected score at 12 months, conditional (re_formula = NA).
wide_summary <- function(brmsfit) {
    b <- effect_from_sim_study(brmsfit)
    dr <- posterior::as_draws(brmsfit)
    em <- matrix(b$effect_samples, ncol = posterior::nchains(dr))
    list(vars = posterior::variables(dr),
         n_draws = posterior::ndraws(dr),
         summ = posterior::summarise_draws(dr, "mean", "sd", "ess_bulk", "ess_tail", "rhat"),
         effect = as.numeric(b$effect_samples),
         ess_effect_bulk = posterior::ess_bulk(em),
         ess_effect_tail = posterior::ess_tail(em),
         divergences = b$diags$n_divergent,
         metric = brmsfit$fit@stan_args[[1]]$metric)
}

# Sim1.Rmd's production budget: 3 chains x (1000 warmup + 2000 draws).
sampler_args <- list(adapt_delta = 0.95, init = 0.1, seed = sampler_seed,
                     chains = 3, iter = 3000, warmup = 1000)
args_per_fit <- list(list(data = d_wide))

results_file <- here::here("local_temp_data", "s2z_stage_c_single_results.rds")
results <- list(started = Sys.time(), questions = QUESTIONS, prior_set = PRIOR_SET,
                n_subj_per_group = n_subj_per_group,
                brms_sha = packageDescription("brms")$RemoteSha,
                sampler_args = sampler_args, arms = list())
checkpoint <- function() saveRDS(results, results_file)

cat("\npre-compiling (one compile per item x parameterisation) ...\n")
for (q in QUESTIONS) for (p in names(GR)) {
    sc <- as.character(make_stancode(mk(q, GR[[p]]), data = d_wide, prior = the_prior))
    invisible(cmdstanr::cmdstan_model(cmdstanr::write_stan_file(sc)))
    if (p == "s2z_nc" && q == QUESTIONS[1]) invisible(nutpieR::nutpie_compile_model(code = sc, verbose = 0))
}
cat("  done\n")

run_arm <- function(item, par, backend) {
    label <- sprintf("%s %s / %s", item, par, backend)
    cat("\n=== ", label, " === started", format(Sys.time(), "%H:%M:%S"), "\n")
    t <- system.time(
        res <- brm_parallel(
            args_shared = c(list(formula = mk(item, GR[[par]]), prior = the_prior), sampler_args),
            args_per_fit = args_per_fit,
            backend = backend, cores_per_fit = cores_per_fit,
            cache_fits = FALSE, cache_summaries = FALSE,
            summarise_fun = wide_summary,
            future.globals = list("wide_summary" = wide_summary,
                                  "effect_from_sim_study" = effect_from_sim_study))
    )
    out <- res[[1]]
    out$wall_sec <- unname(t["elapsed"]); out$label <- label
    out$item <- item; out$par <- par; out$backend <- backend
    cat(sprintf("    %.1f min | %d vars | %d divergences | rhat %.3f\n",
                out$wall_sec / 60, length(out$vars), out$divergences,
                max(out$summ$rhat, na.rm = TRUE)))
    out
}

# ---- the fits: 3 cmdstanr arms per item, plus one nutpie arm on Q05 ---------
for (q in QUESTIONS) {
    for (p in names(GR)) {
        results$arms[[paste(q, p, "cs", sep = "_")]] <- run_arm(q, p, "cmdstanr")
        checkpoint()
    }
}
k <- paste(QUESTIONS[1], "s2z_nc", "np", sep = "_")
results$arms[[k]] <- run_arm(QUESTIONS[1], "s2z_nc", "nutpier"); checkpoint()

a <- results$arms

# ---- 1. conversion ----------------------------------------------------------
cat("\n\n############ 1. CONVERSION (nutpie + s2z, wide) ############\n")
np <- a[[k]]; ref_np <- a[[paste(QUESTIONS[1], "s2z_nc", "cs", sep = "_")]]
cat("  brmsfit built         :", !is.null(np$summ), "\n")
cat("  effect draws          :", length(np$effect), "\n")
cat("  divergences (numeric) :", is.numeric(np$divergences), "|", np$divergences, "\n")
cat("  variables vs cmdstanr :", length(np$vars), "vs", length(ref_np$vars), "\n")

# ---- 2. equivalence ---------------------------------------------------------
# lp__ / lprior / merged_Intercept live in each model's OWN coordinates.
cat("\n############ 2. EQUIVALENCE vs base, per item ############\n")
zcmp <- function(x, y) {
    dplyr::inner_join(x$summ, y$summ, by = "variable", suffix = c("_a", "_ref")) |>
        dplyr::filter(sd_a > 0 | sd_ref > 0,
                      !variable %in% c("lp__", "lprior"),
                      !grepl("^merged_Intercept", variable)) |>
        dplyr::mutate(z = (mean_a - mean_ref) /
                          sqrt(sd_a^2 / ess_bulk_a + sd_ref^2 / ess_bulk_ref)) |>
        dplyr::filter(!is.na(z)) |> dplyr::arrange(dplyr::desc(abs(z)))
}
GAMMA <- "b_alsfrs_dly_mnths:groupTreatment"
gamma_of <- function(x) {
    g <- x$summ[x$summ$variable == GAMMA, ]
    if (nrow(g) != 1) stop("treatment coefficient not found; have: ",
                           paste(grep("groupTreatment", x$summ$variable, value = TRUE), collapse = ", "))
    g
}
eff_z <- function(x, y) (mean(x$effect) - mean(y$effect)) /
    sqrt(sd(x$effect)^2 / x$ess_effect_bulk + sd(y$effect)^2 / y$ess_effect_bulk)

for (q in QUESTIONS) {
    ref <- a[[paste(q, "base", "cs", sep = "_")]]
    for (nm in names(a)) {
        x <- a[[nm]]
        if (x$item != q || x$par == "base") next
        cmp <- zcmp(x, ref); g <- cmp[cmp$variable == GAMMA, ]
        cat("\n##", x$label, "\n")
        cat(sprintf("   %d shared vars | max |z| = %.2f | frac |z|>3 = %.4f (chance 0.0027)\n",
                    nrow(cmp), max(abs(cmp$z)), mean(abs(cmp$z) > 3)))
        cat(sprintf("   gamma : %.4f (sd %.4f) vs base %.4f (sd %.4f) -> z = %.2f\n",
                    g$mean_a, g$sd_a, g$mean_ref, g$sd_ref, g$z))
        cat(sprintf("   effect: %.4f [%.4f, %.4f] vs base %.4f [%.4f, %.4f] -> z = %.2f\n",
                    mean(x$effect), quantile(x$effect, .025), quantile(x$effect, .975),
                    mean(ref$effect), quantile(ref$effect, .025), quantile(ref$effect, .975),
                    eff_z(x, ref)))
        cat("   worst-agreeing:\n"); print(head(cmp[, c("variable", "mean_a", "mean_ref", "z")], 3))
    }
}

# ---- 3. efficiency ----------------------------------------------------------
cat("\n############ 3. EFFICIENCY ############\n")
tab <- do.call(rbind, lapply(names(a), function(nm) {
    x <- a[[nm]]; mins <- x$wall_sec / 60; g <- gamma_of(x)
    data.frame(item = x$item, arm = paste(x$par, x$backend),
               wall_min = round(mins, 1),
               ess_gamma_bulk = round(g$ess_bulk), ess_gamma_tail = round(g$ess_tail),
               bulk_per_min = round(g$ess_bulk / mins, 1),
               tail_per_min = round(g$ess_tail / mins, 1),
               eff_tail_per_min = round(x$ess_effect_tail / mins, 1),
               min_ess_bulk = round(min(x$summ$ess_bulk, na.rm = TRUE)),
               max_rhat = round(max(x$summ$rhat, na.rm = TRUE), 3),
               divergences = x$divergences, stringsAsFactors = FALSE)
}))
tab$ratio_vs_base <- NA_real_
for (q in QUESTIONS) {
    b <- tab$bulk_per_min[tab$item == q & tab$arm == "base cmdstanr"]
    tab$ratio_vs_base[tab$item == q] <- round(tab$bulk_per_min[tab$item == q] / b, 2)
}
print(tab, row.names = FALSE)

cat("\nReference: `hier` gave 0.94 (bulk) / 1.07 (tail) with s2z on all three terms.\n",
    "Here the single constraint removes the only ridge, so this is s2z's best case.\n",
    "Read section 2 first -- a ratio only counts if the posterior is reproduced.\n")

results$table <- tab; results$finished <- Sys.time(); checkpoint()
parallel::stopCluster(cl); future::plan(future::sequential)
cat("\nsaved to", results_file, "\n")
