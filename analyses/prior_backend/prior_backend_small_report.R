# prior_backend_small_report.R  --  replicate-aware, item-aware version
#
# READ-ONLY report over the results of analyses/prior_backend_small.R.
# Works for any number of replicates per cell (including one), and for the
# single-item sections (Q5, Q12) and the twelve-item section (ALL) alike.
# Loads local_temp_data/prior_backend_small/small_<section>.rds and prints
# tables. Fits nothing, writes nothing.
#
# Run from the project directory (here::here() needs it):
#   Rscript analyses/prior_backend_small_report.R                # Q5, Q12 and ALL, detail on replicate 1
#   Rscript analyses/prior_backend_small_report.R ALL            # one section
#   Rscript analyses/prior_backend_small_report.R Q5 focus=3     # detail sections on replicate 3
#
# Output per results file
#   0   what was run, model size
#   C   per-fit table: make_table (copied verbatim) + per-fit cost
#   A   cell summaries over replicates: A1 inference, A2 reliability and cost
#   B   paired contrasts over replicates: prior within backend, backend within prior
#   S   structural parameters: paired contrasts over replicates
#       (per parameter for a single-item fit, per item for the twelve-item fit)
#   K   twelve-item fit only: the cross-item correlation matrix (| p |) by block
#   W   worst parameters and class diagnostics over replicates
#   F   one replicate in detail (the focus replicate)
# With one replicate per cell, A, B and S show that replicate's values (no spread)
# and F is the complete single-replicate report.
#
# Summary conventions
#   mean ± SD   over replicates, for posterior quantities
#   gm [range]  geometric mean and range, for ratios and costs
#   k/n         number of replicates satisfying the condition named in the column
#   sign        replicates on the same side of the null value (0 for shifts, 1 for
#               ratios) as the mean; 5/5 has two-sided probability 1/16 under no effect
#   mc_noise    Monte Carlo noise of a DID shift in SD units, sqrt(1/ESS_a + 1/ESS_b);
#               a shift is only meaningful if it clearly exceeds this
#   a vs b      contrasts are always "a relative to b": shift_sd = (mean_a - mean_b) / sd_b,
#               sd_ratio = sd_a / sd_b, cost ratios = a / b. Reference: prior "default",
#               backend "cmdstanr" (independent of the order of the grid in the run script)
#   item tables the twelve-item fit has the same nine own parameters per item as a
#               single-item fit (4 thresholds, time, group, interaction, 2 sds); item
#               tables show one row per item for a chosen parameter (role):
#               tau1 = first threshold, time, group, int = interaction,
#               sd0 = subject intercept sd, sd1 = subject slope sd

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(purrr)
})
options(width = 250)

cache_dir <- here::here("local_temp_data", "prior_backend_small")
cl_args <- commandArgs(trailingOnly = TRUE)
focus <- 1L
if (any(grepl("^focus=", cl_args))) {
    focus <- as.integer(sub("^focus=", "", cl_args[grepl("^focus=", cl_args)][1]))
    cl_args <- cl_args[!grepl("^focus=", cl_args)]
}
sections <- if (length(cl_args) == 0) c("Q5", "Q12", "ALL") else cl_args

## ---- make_table: copied verbatim from analyses/prior_backend_small.R ---------
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

## ---- formatting helpers ------------------------------------------------------
gm  <- function(x) exp(mean(log(x)))
fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)
rng <- function(x, d = 2) if (length(x) == 0) "NA" else sprintf("[%s, %s]", fmt(min(x), d), fmt(max(x), d))
msd <- function(x, d = 3) {
    if (length(x) == 0) return("NA")
    if (length(x) > 1) sprintf("%s ± %s", fmt(mean(x), d), fmt(sd(x), d)) else fmt(x, d)
}
gmr <- function(x, d = 2) {
    if (length(x) == 0 || any(is.na(x))) return("NA")
    if (length(x) > 1) sprintf("%s %s", fmt(gm(x), d), rng(x, d)) else fmt(x, d)
}
gmf <- function(x, d = 2) if (length(x) == 0 || any(is.na(x))) "NA" else fmt(gm(x), d)
same_side <- function(x, center = 0) {
    if (length(x) == 0) return("NA")
    s <- sign(mean(x) - center)
    sprintf("%d/%d", sum(sign(x - center) == s), length(x))
}
kofn <- function(flag) sprintf("%d/%d", sum(flag, na.rm = TRUE), length(flag))
pr <- function(df) print(as.data.frame(df), row.names = FALSE)
hr <- function(title) cat("\n---", title, "---\n")

## ---- parameter classes, items, roles, correlation blocks ---------------------
param_class <- function(v) {
    dplyr::case_when(
        grepl("^b_(Q\\d\\d_)?Intercept\\[", v) ~ "threshold",
        grepl("^b_", v)                        ~ "fixed",
        grepl("^sd_", v)                       ~ "sd",
        grepl("^cor_", v)                      ~ "cor",
        grepl("^r_", v)                        ~ "random",
        TRUE                                   ~ "internal")
}
class_levels <- c("threshold", "fixed", "sd", "cor", "random", "internal")
is_structural <- function(v) grepl("^(b_|sd_|cor_)", v)
is_own <- function(v) grepl("^(b_|sd_)", v)     # an item's own parameters; correlations excluded

## the item an (own) parameter belongs to: "Q01".."Q12" in a multivariate fit, NA otherwise
item_of <- function(v) {
    m <- regexpr("Q\\d\\d", v)
    out <- rep(NA_character_, length(v))
    out[m > 0] <- regmatches(v, m)
    out
}
## the role of an own parameter within its item
role_of <- function(v) {
    dplyr::case_when(
        grepl("^b_(Q\\d\\d_)?Intercept\\[1\\]$", v)                ~ "tau1",
        grepl("^b_(Q\\d\\d_)?alsfrs_dly_mnths$", v)                ~ "time",
        grepl("^b_(Q\\d\\d_)?groupTreatment$", v)                  ~ "group",
        grepl("^b_(Q\\d\\d_)?alsfrs_dly_mnths:groupTreatment$", v) ~ "int",
        grepl("^sd_.*__(Q\\d\\d_)?Intercept$", v)                  ~ "sd0",
        grepl("^sd_.*__(Q\\d\\d_)?alsfrs_dly_mnths$", v)           ~ "sd1",
        TRUE ~ NA_character_)
}
roles <- c(tau1  = "first threshold, b_Intercept[1]",
           time  = "time slope per month, b_alsfrs_dly_mnths",
           group = "group main effect, b_groupTreatment",
           int   = "interaction per month, b_alsfrs_dly_mnths:groupTreatment",
           sd0   = "subject intercept sd, sd_subject_id__Intercept",
           sd1   = "subject slope sd, sd_subject_id__alsfrs_dly_mnths")
short_name <- function(v) sub("Q\\d\\d_", "", v)   # drop the item prefix from a parameter name

## block of a correlation parameter: cor_<group>__<resp1>_<coef1>__<resp2>_<coef2>
cor_block_of <- function(v) {
    out <- rep(NA_character_, length(v))
    idx <- which(grepl("^cor_", v))
    if (length(idx) == 0) return(out)
    parts <- strsplit(v[idx], "__", fixed = TRUE)
    p1 <- vapply(parts, function(p) p[2], "")
    p2 <- vapply(parts, function(p) p[3], "")
    t1 <- ifelse(grepl("Intercept$", p1), "I", "S")
    t2 <- ifelse(grepl("Intercept$", p2), "I", "S")
    i1 <- item_of(p1); i2 <- item_of(p2)
    same <- (is.na(i1) & is.na(i2)) | (!is.na(i1) & !is.na(i2) & i1 == i2)
    out[idx] <- dplyr::case_when(t1 == "I" & t2 == "I" ~ "int-int (cross item)",
                                 t1 == "S" & t2 == "S" ~ "slope-slope (cross item)",
                                 same                  ~ "int-slope (same item)",
                                 TRUE                  ~ "int-slope (cross item)")
    out
}

## ---- per-fit extraction ------------------------------------------------------

## every successful fit as (key, replicate, x); failed cells are skipped
fit_list <- function(results) {
    out <- list()
    for (key in names(results$cells)) {
        cell <- results$cells[[key]]
        if (!is.null(cell$error)) next
        for (j in seq_along(cell)) {
            x <- cell[[j]]
            rep <- if (is.null(x$replicate)) j else x$replicate
            out[[length(out) + 1]] <- list(key = key, replicate = rep, x = x)
        }
    }
    out
}

## one row per fit with the metrics the summaries are built from
fit_metrics_table <- function(results) {
    map_df(fit_list(results), function(f) {
        x <- f$x
        s <- as.data.frame(x$summ)
        st <- s[is_structural(s$variable), ]
        ok <- !is.na(s$ess_bulk)
        grads <- x$n_draws * x$mean_n_leapfrog
        el <- if (is.null(x$elapsed)) NA_real_ else max(rowSums(x$elapsed))
        q <- quantile(x$effect, c(0.025, 0.975))
        nfit <- length(results$cells[[f$key]])
        w <- results$wall_min[[f$key]]
        data.frame(cell = f$key,
                   prior = if (is.null(x$prior)) sub("_[^_]+$", "", f$key) else x$prior,
                   backend = if (is.null(x$backend)) sub("^.*_", "", f$key) else x$backend,
                   replicate = f$replicate,
                   nchains = x$nchains, n_draws = x$n_draws,
                   eff_mean = mean(x$effect), eff_sd = sd(x$effect),
                   eff_lo = q[[1]], eff_hi = q[[2]],
                   excl0 = (q[[1]] > 0 | q[[2]] < 0),
                   ess_eff = x$ess_effect_bulk, ess_eff_tail = x$ess_effect_tail,
                   rhat_eff = x$rhat_effect,
                   max_rhat = max(st$rhat, na.rm = TRUE),
                   min_ess = min(st$ess_bulk, na.rm = TRUE),
                   n_bad_rhat = sum(st$rhat > 1.01, na.rm = TRUE),
                   worst_par = s$variable[ok][which.min(s$ess_bulk[ok])],
                   div = x$divergences,
                   n_leapfrog = x$mean_n_leapfrog, treedepth = x$mean_treedepth,
                   gradients_M = grads / 1e6,
                   ess_eff_per_Mgrad = x$ess_effect_bulk / grads * 1e6,
                   min_ess_per_Mgrad = min(st$ess_bulk, na.rm = TRUE) / grads * 1e6,
                   sampler_s = el,
                   wall_s_per_fit = if (is.null(w)) NA_real_ else w * 60 / nfit,
                   stringsAsFactors = FALSE)
    })
}

## long table: one row per (cell, replicate, variable) with the summarise_draws columns
cells_long <- function(results) {
    map_df(fit_list(results), function(f) {
        as.data.frame(f$x$summ) %>%
            mutate(cell = f$key, replicate = f$replicate,
                   class = param_class(variable), item = item_of(variable),
                   role = role_of(variable), block = cor_block_of(variable),
                   .before = 1)
    })
}

## the contrasts of the design: every prior against the reference prior within each
## backend, every backend against the reference backend within each prior
contrast_pairs <- function(fm) {
    priors <- unique(fm$prior); backends <- unique(fm$backend)
    ref_p <- if ("default" %in% priors) "default" else priors[1]
    ref_b <- if ("cmdstanr" %in% backends) "cmdstanr" else backends[1]
    out <- list()
    for (b in backends) for (p in setdiff(priors, ref_p)) {
        out[[length(out) + 1]] <- list(name = sprintf("%s vs %s | %s", p, ref_p, b),
                                       a = paste(p, b, sep = "_"), b = paste(ref_p, b, sep = "_"))
    }
    for (p in priors) for (b in setdiff(backends, ref_b)) {
        out[[length(out) + 1]] <- list(name = sprintf("%s vs %s | %s", b, ref_b, p),
                                       a = paste(p, b, sep = "_"), b = paste(p, ref_b, sep = "_"))
    }
    out
}

## cell a against reference cell b, per replicate and variable.
## z on the posterior mean uses MCSE ~ sd / sqrt(ess_bulk) on each side;
## constants (sd = 0, e.g. `disc`) are dropped.
struct_contrast <- function(long, cell_a, cell_b) {
    a <- long %>% filter(cell == cell_a) %>%
        select(replicate, variable, class, item, role, block,
               mean_a = mean, sd_a = sd, ess_a = ess_bulk, rhat_a = rhat)
    b <- long %>% filter(cell == cell_b) %>%
        select(replicate, variable, mean_b = mean, sd_b = sd, ess_b = ess_bulk, rhat_b = rhat)
    inner_join(a, b, by = c("replicate", "variable")) %>%
        filter(sd_a > 0, sd_b > 0, !is.na(ess_a), !is.na(ess_b)) %>%
        mutate(z        = (mean_a - mean_b) / sqrt(sd_a^2 / ess_a + sd_b^2 / ess_b),
               shift_sd = (mean_a - mean_b) / sd_b,
               sd_ratio = sd_a / sd_b,
               ess_ratio = ess_a / ess_b)
}

## summary of a struct_contrast over a set of variables: per replicate, then over replicates
overall_summary <- function(ct, label) {
    if (nrow(ct) == 0) return(NULL)
    ct %>% group_by(replicate) %>%
        summarise(mz = max(abs(z)), f3 = mean(abs(z) > 3),
                  ms = median(abs(shift_sd)), mr = median(sd_ratio), me = median(ess_ratio),
                  .groups = "drop") %>%
        summarise(over = label, n_rep = n(), n_vars = nrow(ct) / n(),
                  max_abs_z = fmt(max(mz), 2), frac_abs_z_gt3 = fmt(mean(f3), 4),
                  median_abs_shift_sd = fmt(mean(ms), 3), median_sd_ratio = fmt(mean(mr), 3),
                  median_ess_ratio = fmt(mean(me), 2))
}
overall_rows <- function(ct, multi_item) {
    rows <- list(overall_summary(ct, "all variables"),
                 overall_summary(ct %>% filter(class == "random"), "random effects"),
                 overall_summary(ct %>% filter(is_structural(variable)), "structural"))
    if (multi_item) rows <- c(rows, list(overall_summary(ct %>% filter(class == "cor"), "correlations")))
    bind_rows(rows)
}

## per-parameter contrast summary over replicates (single-item fits)
param_contrast_summary <- function(ct, vars) {
    ct %>% filter(variable %in% vars) %>%
        mutate(variable = factor(variable, levels = vars)) %>%
        group_by(variable) %>%
        summarise(n = n(),
                  shift = msd(shift_sd, 2), shift_rng = rng(shift_sd, 2), shift_sign = same_side(shift_sd),
                  sdr = fmt(gm(sd_ratio), 2), sdr_rng = rng(sd_ratio, 2), sdr_sign = same_side(log(sd_ratio)),
                  essr = fmt(gm(ess_ratio), 2),
                  max_abs_z = fmt(max(abs(z)), 1),
                  .groups = "drop") %>%
        arrange(variable)
}

## per-item contrast summary over replicates (multi-item fits): each item's own parameters
item_contrast_summary <- function(ct) {
    ct %>% filter(is_own(variable), !is.na(item)) %>%
        group_by(item) %>%
        summarise(sd0_shift  = msd(shift_sd[role %in% "sd0"], 2),
                  sd0_sdr    = gmf(sd_ratio[role %in% "sd0"]),
                  tau1_shift = msd(shift_sd[role %in% "tau1"], 2),
                  tau1_sdr   = gmf(sd_ratio[role %in% "tau1"]),
                  time_shift = msd(shift_sd[role %in% "time"], 2),
                  int_shift  = msd(shift_sd[role %in% "int"], 2),
                  int_sdr    = gmf(sd_ratio[role %in% "int"]),
                  max_abs_z  = fmt(max(abs(z)), 1),
                  ess_ratio  = gmf(ess_ratio),
                  .groups = "drop") %>%
        arrange(item)
}

## effect contrast between two fit_metrics rows (one replicate each)
effect_contrast_row <- function(ra, rb) {
    data.frame(mean_a = fmt(ra$eff_mean, 3), mean_b = fmt(rb$eff_mean, 3),
               z = fmt((ra$eff_mean - rb$eff_mean) /
                           sqrt(ra$eff_sd^2 / ra$ess_eff + rb$eff_sd^2 / rb$ess_eff), 2),
               shift_sd = fmt((ra$eff_mean - rb$eff_mean) / rb$eff_sd, 3),
               sd_ratio = fmt(ra$eff_sd / rb$eff_sd, 3))
}

## variables x cells for one replicate, keeping the Stan order of the variables
wide_by_cell <- function(long_one, value, vars, keys) {
    long_one %>%
        filter(variable %in% vars) %>%
        mutate(variable = factor(variable, levels = vars), cell = factor(cell, levels = keys)) %>%
        select(variable, cell, all_of(value)) %>%
        arrange(cell) %>%
        pivot_wider(names_from = cell, values_from = all_of(value)) %>%
        arrange(variable) %>%
        mutate(variable = as.character(variable))
}

## items x cells, one table per role, posterior mean (sd) on one replicate
item_role_tables <- function(long_one, keys) {
    for (r in names(roles)) {
        d <- long_one %>% filter(role %in% r) %>%
            mutate(v = sprintf("%.2f (%.2f)", mean, sd), cell = factor(cell, levels = keys)) %>%
            select(item, cell, v) %>% arrange(cell) %>%
            pivot_wider(names_from = cell, values_from = v) %>% arrange(item)
        cat("\n  ", r, ":", roles[[r]], "\n")
        pr(d)
        if (r == "tau1") cat("  a first threshold with posterior sd above ~4 under the default prior means the lowest\n",
                             " category is (nearly) unobserved; confirm with table(sims[[1]]$Qxx)\n")
    }
}

## items x cells: min ess_bulk / max rhat over the item's own parameters (worst parameter);
## over all replicates present in long_x
item_diag_table <- function(long_x, keys) {
    long_x %>% filter(is_own(variable), !is.na(item), !is.na(ess_bulk)) %>%
        group_by(item, cell) %>%
        summarise(v = sprintf("%d / %.3f (%s)", round(min(ess_bulk)), max(rhat),
                              short_name(variable[which.min(ess_bulk)])), .groups = "drop") %>%
        mutate(cell = factor(cell, levels = keys)) %>% arrange(cell) %>%
        pivot_wider(names_from = cell, values_from = v) %>% arrange(item)
}

## the correlation matrix by block, per cell, over variables and replicates
cor_summary <- function(long, keys) {
    long %>% filter(class == "cor") %>%
        group_by(cell, block) %>%
        summarise(n_cor = n_distinct(variable),
                  mean_of_means = fmt(mean(.data$mean), 3),
                  range_of_means = rng(.data$mean, 2),
                  mean_post_sd = fmt(mean(.data$sd), 3),
                  min_ess = round(min(ess_bulk, na.rm = TRUE)),
                  max_rhat = fmt(max(rhat, na.rm = TRUE), 3),
                  .groups = "drop") %>%
        mutate(cell = factor(cell, levels = keys)) %>% arrange(block, cell)
}
## correlation blocks under a contrast; mean_abs_change < 0 = correlations pulled toward 0 in a
cor_contrast_summary <- function(ct) {
    ct %>% filter(class == "cor") %>%
        group_by(block) %>%
        summarise(n_cor = n_distinct(variable),
                  mean_abs_change = fmt(mean(abs(mean_a) - abs(mean_b)), 3),
                  median_abs_shift_sd = fmt(median(abs(shift_sd)), 3),
                  median_sd_ratio = fmt(median(sd_ratio), 3),
                  max_abs_z = fmt(max(abs(z)), 2),
                  frac_abs_z_gt3 = fmt(mean(abs(z) > 3), 3),
                  .groups = "drop")
}

## ---- main --------------------------------------------------------------------
for (sec in sections) {
    f <- file.path(cache_dir, paste0("small_", sec, ".rds"))
    cat("\n\n############################################################################\n")
    cat("###", sec, "  ", f, "\n")
    cat("############################################################################\n")
    if (!file.exists(f)) { cat("   (file missing)\n"); next }
    results <- readRDS(f)
    keys <- names(results$cells)
    failed <- keys[vapply(keys, function(k) !is.null(results$cells[[k]]$error), logical(1))]
    fm <- fit_metrics_table(results)
    if (nrow(fm) == 0) { cat("   no successful cells\n"); next }
    fm <- fm %>% mutate(cell = factor(cell, levels = keys)) %>% arrange(cell, replicate) %>%
        mutate(cell = as.character(cell))
    long <- cells_long(results)
    reps <- sort(unique(fm$replicate))
    pairs <- contrast_pairs(fm)
    pairs <- pairs[vapply(pairs, function(p) all(c(p$a, p$b) %in% fm$cell), logical(1))]
    struct_vars <- unique(long$variable[is_structural(long$variable)])

    ## single item or several: items are read from the parameter names
    items <- sort(unique(long$item[is_own(long$variable) & !is.na(long$item)]))
    multi_item <- length(items) > 1
    if (!multi_item) long$item[is.na(long$item)] <- sec

    ## 0 ------------------------------------------------------------------------
    hr("0  what was run")
    cat("sim_seed:", results$sim_seed, "  nsims:", results$nsims,
        "  replicates present:", paste(reps, collapse = ","), "  focus replicate:", focus, "\n")
    cat("args_sh: ", paste(names(results$args_sh), unlist(results$args_sh), sep = " = ", collapse = ", "), "\n")
    cat("cells:   ", paste(keys, collapse = ", "), "\n")
    for (k in failed) cat("FAILED", k, ":", results$cells[[k]]$error, "\n")
    cat("model:   ", if (multi_item) sprintf("%d items (%s)", length(items), paste(items, collapse = ", "))
                     else "single item", "\n")
    cat("variables per class in one fit:\n")
    pr(long %>% filter(cell == keys[keys %in% fm$cell][1], replicate == reps[1]) %>%
           count(class) %>% mutate(class = factor(class, levels = class_levels)) %>% arrange(class) %>%
           pivot_wider(names_from = class, values_from = n))
    cat("WI prior block in force (identical lines across responses collapsed; n_lines = how many):\n")
    if (!is.null(results$priors$WI)) {
        wi <- as.data.frame(results$priors$WI) %>% mutate(o = row_number())
        pr(wi %>% group_by(prior, class, coef, group) %>%
               summarise(n_lines = n(), o = min(o), .groups = "drop") %>% arrange(o) %>% select(-o))
    }
    cat("wall clock per cell (min), the whole brm_parallel call, all replicates of the cell:\n")
    pr(fm %>% distinct(cell, wall_s_per_fit) %>%
           left_join(fm %>% count(cell, name = "n_fits"), by = "cell") %>%
           mutate(wall_min_cell = round(wall_s_per_fit * n_fits / 60, 2),
                  wall_s_per_fit = round(wall_s_per_fit)) %>%
           select(cell, n_fits, wall_min_cell, wall_s_per_fit))

    ## C ------------------------------------------------------------------------
    hr("C  per-fit table (make_table)")
    pr(make_table(results))
    hr("C2 per-fit cost (post-warm-up gradients only; nutpie 400 vs cmdstan 1000 warm-up draws by design)")
    pr(fm %>% transmute(cell, replicate, nchains, n_draws,
                        n_leapfrog = round(n_leapfrog, 1), treedepth = round(treedepth, 2),
                        gradients_M = round(gradients_M, 3),
                        ess_eff_per_Mgrad = round(ess_eff_per_Mgrad),
                        min_ess_per_Mgrad = round(min_ess_per_Mgrad),
                        worst_par))

    ## A ------------------------------------------------------------------------
    hr("A1 inference on the estimand, over replicates (mean ± SD; k/n = replicates whose 95% interval excludes 0)")
    pr(fm %>% group_by(cell) %>%
           summarise(n_fits = n(),
                     DID_mean = msd(eff_mean, 3),
                     DID_sd = msd(eff_sd, 3),
                     sd_of_means_over_mean_sd = if (n() > 1) fmt(sd(eff_mean) / mean(eff_sd), 2) else "n=1",
                     intervals_excl_0 = kofn(excl0),
                     ess_eff_gm = gmr(ess_eff, 0),
                     ess_eff_tail_min = round(min(ess_eff_tail)),
                     rhat_eff_max = fmt(max(rhat_eff), 3),
                     .groups = "drop") %>%
           mutate(cell = factor(cell, levels = keys)) %>% arrange(cell))
    cat("  DID = 12-month difference in differences in points on the scale of the fitted response(s):\n",
        " a single item spans 0-4, the twelve-item sum 0-48\n",
        " sd_of_means_over_mean_sd: empirical SD of the posterior mean across replicates / mean posterior SD;\n",
        " near 1 if the posterior SD is honest; only gross departures are detectable with few replicates\n")

    hr("A2 reliability (worst case over replicates, k/n = fits failing) and cost (gm [range])")
    pr(fm %>% group_by(cell) %>%
           summarise(n_fits = n(),
                     div_max = max(div),
                     fits_with_div = kofn(div > 0),
                     max_rhat_worst = fmt(max(max_rhat), 3),
                     fits_rhat_gt_1.01 = kofn(max_rhat > 1.01),
                     min_ess_worst = round(min(min_ess)),
                     fits_ess_lt_100_per_chain = kofn(min_ess < 100 * nchains),
                     leapfrog_gm = gmr(n_leapfrog, 1),
                     ess_eff_per_Mgrad_gm = gmr(ess_eff_per_Mgrad, 0),
                     min_ess_per_Mgrad_gm = gmr(min_ess_per_Mgrad, 0),
                     sampler_s_gm = gmr(sampler_s, 0),
                     wall_s_per_fit = fmt(first(wall_s_per_fit), 0),
                     .groups = "drop") %>%
           mutate(cell = factor(cell, levels = keys)) %>% arrange(cell))

    ## B ------------------------------------------------------------------------
    hr("B  paired contrasts over replicates (a vs b; shift in units of b's posterior SD; ratios a/b)")
    B <- map_df(pairs, function(p) {
        a <- fm %>% filter(cell == p$a) %>%
            select(replicate, eff_mean_a = eff_mean, eff_sd_a = eff_sd, ess_a = ess_eff,
                   cost_a = ess_eff_per_Mgrad, mcost_a = min_ess_per_Mgrad, div_a = div, rhat_a = max_rhat)
        b <- fm %>% filter(cell == p$b) %>%
            select(replicate, eff_mean_b = eff_mean, eff_sd_b = eff_sd, ess_b = ess_eff,
                   cost_b = ess_eff_per_Mgrad, mcost_b = min_ess_per_Mgrad, div_b = div, rhat_b = max_rhat)
        d <- inner_join(a, b, by = "replicate") %>%
            mutate(shift_sd = (eff_mean_a - eff_mean_b) / eff_sd_b,
                   mc_noise = sqrt(1 / ess_a + 1 / ess_b),
                   width_ratio = eff_sd_a / eff_sd_b,
                   cost_ratio = cost_a / cost_b,
                   mcost_ratio = mcost_a / mcost_b,
                   div_diff = div_a - div_b,
                   rhat_diff = rhat_a - rhat_b)
        if (nrow(d) == 0) return(NULL)
        data.frame(contrast = p$name, n = nrow(d),
                   DID_shift_sd = msd(d$shift_sd, 3),
                   shift_rng = rng(d$shift_sd, 3),
                   shift_sign = same_side(d$shift_sd),
                   mc_noise = fmt(mean(d$mc_noise), 3),
                   DID_width_ratio = fmt(gm(d$width_ratio), 3),
                   width_rng = rng(d$width_ratio, 3),
                   width_sign = same_side(log(d$width_ratio)),
                   ess_per_grad_ratio = fmt(gm(d$cost_ratio), 2),
                   ess_ratio_rng = rng(d$cost_ratio, 2),
                   min_ess_per_grad_ratio = fmt(gm(d$mcost_ratio), 2),
                   div_diff_rng = rng(d$div_diff, 0),
                   fewer_div = kofn(d$div_diff < 0),
                   max_rhat_diff_rng = rng(d$rhat_diff, 3),
                   stringsAsFactors = FALSE)
    })
    pr(B)

    ## S ------------------------------------------------------------------------
    cts <- lapply(pairs, function(p) struct_contrast(long, p$a, p$b))
    hr(if (multi_item) "S  structural parameters per item: paired contrasts over replicates (each item's own parameters)"
       else "S  structural parameters: paired contrasts over replicates")
    for (i in seq_along(pairs)) {
        cat("\n  ", pairs[[i]]$name, "\n")
        if (multi_item) pr(item_contrast_summary(cts[[i]]))
        else pr(param_contrast_summary(cts[[i]], struct_vars))
        pr(overall_rows(cts[[i]], multi_item))
    }
    cat("  shift = mean shift over replicates in units of b's posterior SD; sdr = sd_a/sd_b; ess_ratio = ess_a/ess_b;\n",
        " max_abs_z = largest Monte-Carlo z over replicates (differences beyond MC noise are expected between priors)\n")

    ## K ------------------------------------------------------------------------
    if (multi_item) {
        hr("K1 cross-item correlation matrix (| p |) by block, per cell, over variables and replicates")
        pr(cor_summary(long, keys))
        hr("K2 correlation blocks under each contrast (mean_abs_change < 0 = correlations pulled toward 0 in a)")
        for (i in seq_along(pairs)) {
            cat("\n  ", pairs[[i]]$name, "\n")
            pr(cor_contrast_summary(cts[[i]]))
        }
    }

    ## W ------------------------------------------------------------------------
    if (multi_item) {
        hr("W0 per item x cell: min ess_bulk / max rhat over the item's own parameters and all replicates (worst parameter)")
        pr(item_diag_table(long, keys))
    }
    hr("W1 the parameter with the lowest ess_bulk in each fit, counted over replicates")
    pr(fm %>% count(cell, worst_par, name = "times_worst") %>%
           mutate(cell = factor(cell, levels = keys)) %>% arrange(cell, desc(times_worst)))

    hr("W2 eight lowest parameters by min ess_bulk over replicates, per cell:  name [min ess_bulk, max rhat]")
    w2 <- long %>% filter(!is.na(ess_bulk)) %>%
        group_by(cell, variable) %>%
        summarise(min_ess = min(ess_bulk), worst_rhat = max(rhat), .groups = "drop") %>%
        mutate(cell = factor(cell, levels = keys)) %>%
        group_by(cell) %>% slice_min(min_ess, n = 8, with_ties = FALSE) %>%
        summarise(worst = paste(sprintf("%s [%d, %.3f]", variable, round(min_ess), worst_rhat), collapse = "; "),
                  .groups = "drop") %>% arrange(cell)
    for (i in seq_len(nrow(w2))) cat("   ", as.character(w2$cell[i]), ":", w2$worst[i], "\n")

    hr("W3 diagnostics by parameter class, worst over replicates (threshold = b_Intercept[k]; fixed = other b_; random = r_; internal = z_, L_, Intercept[k], disc, lp__, lprior)")
    pr(long %>% group_by(cell, replicate, class) %>%
           summarise(ce = min(ess_bulk, na.rm = TRUE), cr = max(rhat, na.rm = TRUE),
                     nb = sum(rhat > 1.01, na.rm = TRUE), .groups = "drop") %>%
           group_by(cell, class) %>%
           summarise(n_fits = n(), min_ess_worst = round(min(ce)), max_rhat_worst = fmt(max(cr), 3),
                     n_rhat_gt_1.01_mean = fmt(mean(nb), 1), .groups = "drop") %>%
           mutate(class = factor(class, levels = class_levels), cell = factor(cell, levels = keys)) %>%
           arrange(class, cell))

    ## F ------------------------------------------------------------------------
    focus_use <- if (focus %in% reps) focus else reps[1]
    hr(sprintf("F  replicate %d in detail%s", focus_use,
               if (focus_use != focus) sprintf(" (focus=%d not present)", focus) else ""))
    long_f <- long %>% filter(replicate == focus_use)
    fm_f <- fm %>% filter(replicate == focus_use)
    fits_f <- list()
    for (fx in fit_list(results)) if (fx$replicate == focus_use) fits_f[[fx$key]] <- fx$x

    if (multi_item) {
        hr("F1 items x cells: posterior mean (sd), one table per parameter role")
        item_role_tables(long_f, keys)
        hr("F2 items x cells: min ess_bulk / max rhat over the item's own parameters (worst parameter)")
        pr(item_diag_table(long_f, keys))
    } else {
        hr("F1 structural parameters: posterior mean (sd) per cell")
        pr(wide_by_cell(long_f %>% mutate(v = sprintf("%.3f (%.3f)", mean, sd)), "v", struct_vars, keys))
        hr("F2 structural parameters: ess_bulk / ess_tail / rhat per cell")
        pr(wide_by_cell(long_f %>% mutate(v = sprintf("%d / %d / %.3f", round(ess_bulk), round(ess_tail), rhat)),
                        "v", struct_vars, keys))
    }

    hr("F3 contrasts on this replicate (a vs b)")
    for (p in pairs) {
        ct <- struct_contrast(long_f, p$a, p$b)
        cat("\n  ", p$name, "\n")
        if (multi_item) {
            pr(item_contrast_summary(ct))
        } else {
            pr(ct %>% filter(variable %in% struct_vars) %>%
                   mutate(variable = factor(variable, levels = struct_vars)) %>% arrange(variable) %>%
                   transmute(variable = as.character(variable),
                             mean_b = round(mean_b, 3), mean_a = round(mean_a, 3),
                             shift_sd = round(shift_sd, 3), sd_ratio = round(sd_ratio, 3), z = round(z, 2),
                             ess_b = round(ess_b), ess_a = round(ess_a),
                             rhat_b = round(rhat_b, 3), rhat_a = round(rhat_a, 3)))
        }
        pr(overall_rows(ct, multi_item))
        cat("  ten largest movers (|shift_sd|) over all parameters:\n")
        pr(ct %>% arrange(desc(abs(shift_sd))) %>% slice_head(n = 10) %>%
               transmute(variable, class, mean_b = round(mean_b, 3), mean_a = round(mean_a, 3),
                         shift_sd = round(shift_sd, 3), sd_ratio = round(sd_ratio, 3),
                         ess_ratio = round(ess_ratio, 3)))
        cat("  effect:\n")
        pr(effect_contrast_row(fm_f %>% filter(cell == p$a), fm_f %>% filter(cell == p$b)))
    }

    hr("F4 diagnostics by parameter class on this replicate")
    pr(long_f %>% group_by(cell, class) %>%
           summarise(n = n(),
                     min_ess_bulk = round(min(ess_bulk, na.rm = TRUE)),
                     med_ess_bulk = round(median(ess_bulk, na.rm = TRUE)),
                     min_ess_tail = round(min(ess_tail, na.rm = TRUE)),
                     max_rhat = round(max(rhat, na.rm = TRUE), 3),
                     n_rhat_gt_1.01 = sum(rhat > 1.01, na.rm = TRUE),
                     .groups = "drop") %>%
           mutate(class = factor(class, levels = class_levels), cell = factor(cell, levels = keys)) %>%
           arrange(class, cell))

    cat("\n  eight worst parameters by ess_bulk, per cell:  name [ess_bulk, rhat]\n")
    w <- long_f %>% filter(!is.na(ess_bulk)) %>% mutate(cell = factor(cell, levels = keys)) %>%
        group_by(cell) %>% slice_min(ess_bulk, n = 8, with_ties = FALSE) %>%
        summarise(worst = paste(sprintf("%s [%d, %.3f]", variable, round(ess_bulk), rhat), collapse = "; "),
                  .groups = "drop") %>% arrange(cell)
    for (i in seq_len(nrow(w))) cat("   ", as.character(w$cell[i]), ":", w$worst[i], "\n")

    cat("\n  eight worst parameters by rhat, per cell:  name [rhat, ess_bulk]\n")
    w <- long_f %>% filter(!is.na(rhat)) %>% mutate(cell = factor(cell, levels = keys)) %>%
        group_by(cell) %>% slice_max(rhat, n = 8, with_ties = FALSE) %>%
        summarise(worst = paste(sprintf("%s [%.3f, %d]", variable, rhat, round(ess_bulk)), collapse = "; "),
                  .groups = "drop") %>% arrange(cell)
    for (i in seq_len(nrow(w))) cat("   ", as.character(w$cell[i]), ":", w$worst[i], "\n")

    hr("F5 the estimand (12-month DID in points, conditional u = 0) on this replicate")
    pr(map_df(names(fits_f), function(k) {
        x <- fits_f[[k]]
        q <- quantile(x$effect, c(0.025, 0.5, 0.975))
        data.frame(cell = k, mean = round(mean(x$effect), 3), sd = round(sd(x$effect), 3),
                   q2.5 = round(q[[1]], 3), q50 = round(q[[2]], 3), q97.5 = round(q[[3]], 3),
                   ess_bulk = round(x$ess_effect_bulk), ess_tail = round(x$ess_effect_tail),
                   rhat = round(x$rhat_effect, 3))
    }))

    hr("F6 timing on this replicate: cmdstan sampler clock per chain (s); nutpie fits carry no clock")
    pr(map_df(names(fits_f), function(k) {
        x <- fits_f[[k]]
        if (is.null(x$elapsed)) return(data.frame(cell = k, chain = "-", warmup_s = NA_real_, sample_s = NA_real_))
        el <- x$elapsed
        data.frame(cell = k, chain = rownames(el), warmup_s = round(el[, 1], 1), sample_s = round(el[, 2], 1))
    }))
}
cat("\nDONE_REPORT (read-only)\n")
