# Render the full summary for any s2z / nutpie results file. READ ONLY --
# nothing is refitted and nothing is written. Everything comes from the saved
# .rds, so this runs in seconds against the project brms; it does NOT need the
# s2z library on .libPaths().
#
#   Rscript analyses/s2z_report.R                       # 12 items, mvbind minus |p|
#   Rscript analyses/s2z_report.R single                # Q05 and Q12, one item at a time
#   Rscript analyses/s2z_report.R hier                  # the long-format hier model
#   Rscript analyses/s2z_report.R local_temp_data/x.rds # any file
#
# Handles all three file shapes: arms may or may not carry $item, $par, $tag,
# and the treatment effect may be one coefficient (hier: b_time_treat; single
# item: b_alsfrs_dly_mnths:groupTreatment) or twelve (mvbind: one per response).

suppressMessages({library(dplyr)})
options(width = 200)

arg <- commandArgs(trailingOnly = TRUE)
FILE <- switch(if (length(arg)) arg[1] else "wide",
               wide   = "local_temp_data/s2z_stage_c_wide_results.rds",
               single = "local_temp_data/s2z_stage_c_single_results.rds",
               hier   = "local_temp_data/s2z_stage_c_results.rds",
               arg[1])

r <- readRDS(here::here(FILE)); a <- r$arms
hdr <- function(x) cat("\n", strrep("=", 92), "\n", x, "\n", strrep("=", 92), "\n", sep = "")

# quantities computed in each model's OWN coordinates -- never comparable across
# parameterisations (merged_Intercept holds SHIFTED thresholds under s2z; the
# public reconstructed ones are b_Intercept[...] and those do agree)
EXCL <- function(v) v %in% c("lp__", "lprior") | grepl("^merged_Intercept", v)

item_of <- function(x) if (!is.null(x$item)) x$item else "(all)"
par_of  <- function(x) if (!is.null(x$par)) x$par else sub(" .*$", "", x$label)
tag_of  <- function(x) if (is.null(x$tag) || !nzchar(x$tag)) "" else paste0(" ", x$tag)
armlab  <- function(x) paste0(par_of(x), " ", x$backend, tag_of(x))
is_base <- function(x) grepl("^base", par_of(x))
# treatment effect coefficient(s): hier has one, single item one, mvbind twelve
trt_vars <- function(x) {
    v <- x$summ$variable
    c(grep("^b_time_treat$", v, value = TRUE),
      grep("groupTreatment$", grep("alsfrs_dly_mnths:", v, value = TRUE), value = TRUE))
}

hdr(paste("FILE:", FILE))
cat("run      :", format(r$started), "->", format(r$finished),
    sprintf("(%.1f h)\n", as.numeric(difftime(r$finished, r$started, units = "hours"))))
if (!is.null(r$updated))   cat("updated  :", format(r$updated), "\n")
if (!is.null(r$prior_set)) cat("priors   :", r$prior_set,
    if (identical(r$prior_set, "brms_default")) "(as Sim1.Rmd: class b is FLAT)" else "", "\n")
if (!is.null(r$n_subj_per_group)) cat("n/group  :", r$n_subj_per_group, "\n")
if (!is.null(r$sampler_args)) cat("sampler  :",
    paste(names(r$sampler_args), unlist(r$sampler_args), sep = "=", collapse = "  "), "\n")
if (!is.null(r$brms_sha)) cat("brms     :", substr(r$brms_sha, 1, 8), "\n")
cat("arms     :", paste(names(a), collapse = ", "), "\n")

hdr("1. EFFICIENCY (as saved)")
print(r$table, row.names = FALSE)
cat("\nnutpie arms run 400 warmup draws by default against cmdstan's",
    if (!is.null(r$sampler_args$warmup)) r$sampler_args$warmup else "1000",
    "-- that asymmetry is deliberate\nand is part of any nutpie ratio unless num_warmup was set explicitly.\n")

hdr("2. THE ESTIMAND (conditional DID at 12 months, re_formula = NA)")
for (nm in names(a)) {
    x <- a[[nm]]
    cat(sprintf("  %-26s %7.4f [%7.4f, %7.4f]  sd %6.4f | ESS bulk %5.0f tail %5.0f | div %3d | rhat %.3f\n",
                nm, mean(x$effect), quantile(x$effect, .025), quantile(x$effect, .975), sd(x$effect),
                x$ess_effect_bulk, x$ess_effect_tail, x$divergences, max(x$summ$rhat, na.rm = TRUE)))
}

hdr("3. EQUIVALENCE vs base -- MEANS and SPREADS")
cat("The mean-based z is the usual check. The SD ratio is the one that matters for a\n",
    "coverage study: a sampler can reproduce every posterior mean while systematically\n",
    "narrowing the spreads, and the mean test passes that.\n", sep = "")
for (it in unique(vapply(a, item_of, character(1)))) {
    same <- a[vapply(a, function(x) identical(item_of(x), it), logical(1))]
    bi <- which(vapply(same, is_base, logical(1)))
    if (!length(bi)) { cat("\n(no base arm for", it, ")\n"); next }
    ref <- same[[bi[1]]]
    cat("\n---- item:", it, " | base =", names(same)[bi[1]], "----\n")
    for (nm in setdiff(names(same), names(same)[bi[1]])) {
        x <- same[[nm]]
        cmp <- inner_join(x$summ, ref$summ, by = "variable", suffix = c("_a", "_ref")) |>
            filter(sd_a > 0, sd_ref > 0, !EXCL(variable)) |>
            mutate(z = (mean_a - mean_ref) / sqrt(sd_a^2/ess_bulk_a + sd_ref^2/ess_bulk_ref),
                   sd_ratio = sd_a / sd_ref,
                   sd_se = sqrt(1/(2*ess_bulk_a) + 1/(2*ess_bulk_ref)),
                   sd_z = (sd_ratio - 1) / sd_se)
        tv <- intersect(trt_vars(x), cmp$variable)
        ez <- (mean(x$effect) - mean(ref$effect)) /
              sqrt(sd(x$effect)^2/x$ess_effect_bulk + sd(ref$effect)^2/ref$ess_effect_bulk)
        cat(sprintf("\n  %s  (%d shared vars)\n", nm, nrow(cmp)))
        cat(sprintf("    MEANS  max|z| %5.2f | frac|z|>3 %.4f (chance 0.0027) | %d treatment coef(s) max|z| %5.2f\n",
                    max(abs(cmp$z)), mean(abs(cmp$z) > 3), length(tv),
                    if (length(tv)) max(abs(cmp$z[cmp$variable %in% tv])) else NA))
        cat(sprintf("    SPREAD median sd_ratio %.3f | %d of %d with |sd_z| > 3 | effect sd_ratio %.3f\n",
                    median(cmp$sd_ratio), sum(abs(cmp$sd_z) > 3), nrow(cmp),
                    sd(x$effect) / sd(ref$effect)))
        w <- cmp |> arrange(desc(abs(sd_z))) |> head(3)
        cat("    worst SD:", paste(sprintf("%s(%.2f, %.1f SE)", w$variable, w$sd_ratio, w$sd_z),
                                   collapse = "  "), "\n")
    }
}

hdr("4. TRAJECTORY LENGTH / COST PER ITERATION")
nch <- if (!is.null(r$sampler_args$chains)) r$sampler_args$chains else 3
warm_of <- function(x) {
    if (!identical(x$backend, "nutpier")) return(as.numeric(r$sampler_args$warmup))
    if (!is.null(x$nutpie_args$num_warmup)) as.numeric(x$nutpie_args$num_warmup) else 400
}
cat("Identical model and data, so cost per gradient evaluation is constant across arms;\n",
    "a large difference in s/iteration can only mean different trajectory lengths.\n\n", sep = "")
for (nm in names(a)) {
    x <- a[[nm]]; it <- warm_of(x) + x$n_draws / nch
    lf <- if (!is.null(x$mean_n_leapfrog)) sprintf("%7.1f (measured)", x$mean_n_leapfrog) else "      - (not recorded)"
    cat(sprintf("  %-26s warmup %4.0f | %8.4f s/iter | n_leapfrog %s\n",
                nm, warm_of(x), x$wall_sec / (nch * it), lf))
}

hdr("5. WORST-MIXING VARIABLES PER ARM (is the bad one a nuisance or the estimand?)")
for (nm in names(a)) {
    x <- a[[nm]]
    s <- x$summ |> filter(!is.na(ess_bulk)) |> arrange(ess_bulk) |> head(3)
    cat(sprintf("\n  %s   (rhat>1.01: %d, rhat>1.05: %d)\n", nm,
                sum(x$summ$rhat > 1.01, na.rm = TRUE), sum(x$summ$rhat > 1.05, na.rm = TRUE)))
    print(as.data.frame(s |> mutate(across(where(is.numeric), ~round(.x, 3)))), row.names = FALSE)
}

# per-response treatment coefficients: only meaningful for the mvbind file
tv0 <- trt_vars(a[[1]])
if (length(tv0) > 1) {
    hdr("6. PER-ITEM TREATMENT COEFFICIENTS (12 responses), each arm vs base")
    ref <- a[[which(vapply(a, is_base, logical(1)))[1]]]
    for (nm in setdiff(names(a), names(a)[vapply(a, is_base, logical(1))])) {
        cmp <- inner_join(a[[nm]]$summ, ref$summ, by = "variable", suffix = c("_a", "_ref")) |>
            filter(variable %in% tv0) |>
            mutate(item = sub("^b_(Q[0-9]+)_.*$", "\\1", variable),
                   z = (mean_a - mean_ref) / sqrt(sd_a^2/ess_bulk_a + sd_ref^2/ess_bulk_ref)) |>
            select(item, mean_a, mean_ref, z, ess_a = ess_bulk_a, ess_ref = ess_bulk_ref) |>
            mutate(across(where(is.numeric), ~round(.x, 3)))
        cat("\n---", nm, "vs base ---\n"); print(as.data.frame(cmp), row.names = FALSE)
    }
}

hdr("DONE -- nothing was modified")
