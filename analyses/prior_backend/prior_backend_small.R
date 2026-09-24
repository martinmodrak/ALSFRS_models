#PRIOR VS BACKEND SMALL SCALE
#testing the weakly informative vs. default prior in small models
#varies between two time points and the standard (roughly 6 time points) settings, logit_rslope models
#first 1 question, Q5 - not expecting a big influence of WI priors
#second 1 questiom, Q12 - expecting bigger influence of WI priors
#third all questions

# LOAD LIBRARIES
library(dplyr)
library(tidyr)
library(brms)
library(readr)
library(posterior)
library(nutpieR)

# SOURCES
source(here::here("R", "data_prep.R"))
source(here::here("R", "simulate_data.R"))
source(here::here("R", "effect_functions.R"))
source(here::here("R", "nutpie_to_stanfit.R"))
source(here::here("R", "sampling_parallel.R"))
source(here::here("R", "brm_parallel.R"))

# PARAMETERS
sim_seed <- 468652233
nsims <-  1
sampler_seed <- 1
n_subj_per_group <- 50
question_col_names <- sprintf("Q%02d", 1:12)
questions_to_fit <- question_col_names
cache_dir <- here::here("local_temp_data", "prior_backend_small")

if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

# SAMPLING ARGUMENTS
args_sh <- list(adapt_delta = 0.95,
                init = 0.1,
                seed = sampler_seed,
                chains = 3,
                iter = 3000,
                warmup = 1000)
cores_per_fit <- 3

# DATA
proact <- read_csv(here::here("private_data/F_PROACT_ALSFRS.csv"), show_col_types = FALSE)
proact_std <- data_prep(proact)
set.seed(sim_seed)
sims <- list()
for(i in 1:nsims) {
    sims[[i]] <- simulate_data_from_registry(proact_std,
                                             max_duration = 12,
                                             min_measurements_per_subject = 4,
                                             max_measurements_per_subject = 8,
                                             n_subjects_per_group = n_subj_per_group,
                                             effect_prob = 0.5) #exploring the null
}

sims_two <- purrr::map(sims, ~ dplyr::filter(.x, is_first | is_last)) #two time points only, first and last

#### 1 QUESTION, Q5 (we know Q5 is well-conditioned so no significant improvement is expected with weakly informative priors)
questions_to_fit <- question_col_names[5]
#MODEL
logit_formula_rslope <- as.formula(paste0(questions_to_fit,
                                          " ~ 1 + alsfrs_dly_mnths*group + (1 + alsfrs_dly_mnths | subject_id)"))
model_formula <- bf(logit_formula_rslope, family = cumulative())

#DATA
args_per_fit <- list()
for(i in 1:nsims) {
    data_for_fit <- sims[[i]] %>%
        mutate(across(all_of(question_col_names), function(x) { x + 1 }))
    args_per_fit[[i]] <- list(data = data_for_fit)
}

#DESIGN
prior_WI <- c(
    set_prior("normal(0, 10)", class = "Intercept"),
    set_prior("normal(0, 1)", class = "b", coef = "alsfrs_dly_mnths"),
    set_prior("normal(0, 3)", class = "b", coef = "groupTreatment"),
    set_prior("normal(0, 1)", class = "b", coef = "alsfrs_dly_mnths:groupTreatment"),
    set_prior("normal(0, 10)", class = "sd", coef = "Intercept", group = "subject_id"),
    set_prior("normal(0, 1)", class = "sd", coef = "alsfrs_dly_mnths", group = "subject_id"),
    set_prior("lkj(2)", class = "cor"))

priors <- list(default = NULL, WI = prior_WI)
backends = c("nutpier", "cmdstanr")
grid <- expand.grid(prior = names(priors), backend = backends, stringsAsFactors = FALSE)

#FIT FUNCTIOIN WITH DIAGNOSTICS
fit_summary <- function(fit, effect_fun) {
    eff <- effect_fun(fit)$effect_samples #estimand
    dr <- posterior::as_draws(fit)
    em <- matrix(eff, ncol = posterior::nchains(dr))
    sp <- do.call(rbind, rstan::get_sampler_params(fit$fit, inc_warmup = FALSE))
    list(vars = posterior::variables(dr),
         n_draws = posterior::ndraws(dr),
         nchains = posterior::nchains(dr),
         summ = posterior::summarise_draws(dr, "mean", "sd", "ess_bulk", "ess_tail", "rhat"),
         effect = as.numeric(eff),
         ess_effect_bulk = posterior::ess_bulk(em),
         ess_effect_tail = posterior::ess_tail(em),
         rhat_effect = posterior::rhat(em),
         divergences = sum(rstan::get_divergent_iterations(fit$fit)),
         mean_n_leapfrog = mean(sp[, "n_leapfrog__"]),
         mean_treedepth = mean(sp[, "treedepth__"]),
         elapsed = tryCatch(rstan::get_elapsed_time(fit$fit), error = function(e) NULL))
}

#RESULTS FUNCTION
get_results <- function(priors_par,
                        model_par,
                        args_per_fit,
                        grid,
                        results_file
                        ) {

   results <- list(sim_seed = sim_seed,
              nsims = nsims,
              args_sh = args_sh,
              priors = priors_par,
              cells = list(),
              wall_min = list())

    for (i in seq_len(nrow(grid))) {

        p <- grid$prior[i]
        b <- grid$backend[i]
        key <- paste(p, b, sep = "_")
        cat("\n===", key, "=== started", format(Sys.time(), "%H:%M:%S"), "\n") #info while running

        t <- system.time(
            cell <- tryCatch(
                brm_parallel(
                    args_shared = c(list(formula = model_par,
                                         prior = priors_par[[p]],
                                         summarise_fun_args = list(effect_fun = effect_from_sim_study)),
                                    args_sh),
                    args_per_fit = args_per_fit,
                    backend = b,
                    cores_per_fit = cores_per_fit,
                    cache_dir = cache_dir, cache_fits = TRUE, cache_summaries = TRUE,
                    summarise_fun = fit_summary,
                    future.globals = list("fit_summary" = fit_summary, "effect_from_sim_study" = effect_from_sim_study)),
                error = function(e) e))

        if (inherits(cell, "error")) { #if brm_parallel ends with error this message pops up
            cat("FAILED:", conditionMessage(cell), "\n")
            results$cells[[key]] <- list(error = conditionMessage(cell))
            next
        }
        for (j in seq_along(cell)) { #
            cell[[j]]$prior <- p
            cell[[j]]$backend <- b
            cell[[j]]$replicate <- j}
        results$cells[[key]] <- cell
        results$wall_min[[key]] <- unname(t["elapsed"]) / 60
        cat(sprintf("    %.1f min for %d replicate(s)\n", results$wall_min[[key]], length(cell)))
        saveRDS(results, results_file)
         }
   results
}

#CREATE TABLE FUNCTION
make_table <- function(results) {
    table <- purrr::map_df(names(results$cells), function(key) {
        cell <- results$cells[[key]]
        if (!is.null(cell$error)) return(data.frame(cell = key, error = cell$error)) #show errors explicitely
        purrr::map_df(cell, function(x) {
            s <- as.data.frame(x$summ)
            st <- s[grepl("^(b_|sd_|cor_)", s$variable), ] #only structural parameters: thresholds, b, sd, cor
            el <- if (is.null(x$elapsed)) NA_real_ else max(rowSums(x$elapsed))
            # outer clock: the whole brm_parallel call for this cell (compile-cache
            # lookup, worker start, sampling, conversion, summary). Per fit only while
            # a cell holds one fit run on its own; with nsims > 1 it is the cell total.
            wall <- results$wall_min[[key]] * 60
            data.frame(cell = key,
                       replicate = x$replicate,
                       eff_mean = round(mean(x$effect), 3),
                       eff_sd = round(sd(x$effect), 3),
                       ess_eff = round(x$ess_effect_bulk),
                       rhat_eff = round(x$rhat_effect, 3),
                       max_rhat = round(max(st$rhat, na.rm = TRUE), 3),
                       min_ess = round(min(st$ess_bulk, na.rm = TRUE)),
                       med_ess = round(median(st$ess_bulk, na.rm = TRUE)),
                       div = x$divergences,
                       sampler_s = round(el),
                       wall_s = round(wall),
                       ess_eff_per_s = round(x$ess_effect_bulk / el, 2),
                       ess_eff_per_wall_s = round(x$ess_effect_bulk / wall, 2))
        })
    })
    return(table)
}

#RESULTS, TABLE
results_Q5 <- get_results(priors, model_formula, args_per_fit, grid, file.path(cache_dir, "small_Q5.rds"))
table_Q5 <- make_table(results_Q5)
print(table_Q5)


###1 QUESTION ONLY, Q12 (we know that the tresholds in Q12 are not identified easily, so here, we expect the WI priors to help)
questions_to_fit <- question_col_names[12]
logit_formula_rslope <- as.formula(paste0(questions_to_fit,
                                          " ~ 1 + alsfrs_dly_mnths*group + (1 + alsfrs_dly_mnths | subject_id)"))
model_formula <- bf(logit_formula_rslope, family = cumulative())

#RESULTS, TABLE
results_Q12 <- get_results(priors, model_formula, args_per_fit, grid, file.path(cache_dir, "small_Q12.rds"))
table_Q12 <- make_table(results_Q12)
print(table_Q12)


###1 QUESTION ONLY, TRY OUT THE REMAINING QUESTIONS
# for(i in 1:11)
# {
#     questions_to_fit <- question_col_names[i]
#     logit_formula_rslope <- as.formula(paste0(questions_to_fit,
#                                               " ~ 1 + alsfrs_dly_mnths*group + (1 + alsfrs_dly_mnths | subject_id)"))
#     model_formula <- bf(logit_formula_rslope, family = cumulative())
#
#     #RESULTS, TABLE
#     results_Q12 <- get_results(priors, model_formula, args_per_fit, grid, file.path(cache_dir, "small_Q12.rds"))
#     table_Q12 <- make_table(results_Q12)
#     print(table_Q12)
#
# }


# ALL QUESTIONS
questions_to_fit <- question_col_names
logit_formula_rslope <- as.formula(paste0(
    "mvbind(", paste0(questions_to_fit, collapse = ", "),
    ") ~ 1 + alsfrs_dly_mnths*group + (1 + alsfrs_dly_mnths | p | subject_id)"))
model_formula <- bf(logit_formula_rslope, family = cumulative()) + set_rescor(FALSE)

#DESIGN
prior_WI <- c(
    do.call(c, lapply(questions_to_fit, function(q) c(
        set_prior("normal(0, 10)", class = "Intercept", resp = q),
        set_prior("normal(0, 1)", class = "b", coef = "alsfrs_dly_mnths", resp = q),
        set_prior("normal(0, 3)", class = "b", coef = "groupTreatment", resp = q),
        set_prior("normal(0, 1)", class = "b", coef = "alsfrs_dly_mnths:groupTreatment", resp = q),
        set_prior("normal(0, 10)", class = "sd", coef = "Intercept", group = "subject_id", resp = q),
        set_prior("normal(0, 1)",  class = "sd", coef = "alsfrs_dly_mnths", group = "subject_id", resp = q)))),
    set_prior("lkj(2)", class = "cor"))

priors <- list(default = NULL, WI = prior_WI)

#RESULTS, TABLE
results_ALL <- get_results(priors, model_formula, args_per_fit, grid, file.path(cache_dir, "small_ALL.rds"))
table_ALL <- make_table(results_ALL)
print(table_ALL)
