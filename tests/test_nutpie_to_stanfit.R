# End-to-end test of nutpie_draws_to_stanfit(): fits the same brms-generated
# ordinal mixed model via (A) the current cmdstanr -> read_csv_as_stanfit path
# and (B) nutpieR -> nutpie_draws_to_stanfit, pushes BOTH through the full brms
# pipeline (emptyfit + rename_pars), and compares the results.
# Rerun this whenever nutpieR, brms, or rstan versions change.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({
    library(brms); library(posterior); library(dplyr)
})
library(nutpieR)
source(here::here("R", "nutpie_to_stanfit.R"))
options(brms.backend = "cmdstanr", width = 200)

## ---- small ordinal mixed model with the features that matter ----
## (item-style thresholds, a random effect, brms's intercept-centering GQ)
set.seed(2026)
n_g <- 60; n_per <- 5
d <- data.frame(g = rep(seq_len(n_g), each = n_per),
                x = rnorm(n_g * n_per))
re <- rnorm(n_g, 0, 0.8)
eta <- 0.8 * d$x + re[d$g] + rlogis(nrow(d))
d$score <- cut(eta, c(-Inf, -1.5, -0.5, 0.5, 1.5, Inf), labels = FALSE)

f <- bf(score | thres(4) ~ x + (1 | g), family = cumulative("logit"))
pr <- c(prior(normal(0, 2), class = "Intercept"),
        prior(normal(0, 1), class = "b"),
        prior(student_t(3, 0, 2.5), class = "sd"))

stancode <- make_stancode(f, data = d, prior = pr)
sdata <- make_standata(f, data = d, prior = pr)
class(sdata) <- NULL
emptyfit <- brm(f, data = d, prior = pr, empty = TRUE)

## ---- path A: current pipeline (cmdstanr -> read_csv_as_stanfit) ----
mod_cs <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(as.character(stancode)))
t_cs <- system.time(
    fit_cs_raw <- mod_cs$sample(data = sdata, chains = 2, parallel_chains = 2,
                                iter_warmup = 500, iter_sampling = 500,
                                init = 0.1, seed = 11, refresh = 0)
)
stanfit_cs <- brms::read_csv_as_stanfit(fit_cs_raw$output_files())
bfit_cs <- emptyfit; bfit_cs$fit <- stanfit_cs; bfit_cs <- rename_pars(bfit_cs)

## ---- path B: nutpie -> nutpie_draws_to_stanfit ----
mod_np <- nutpie_compile_model(code = as.character(stancode))
t_np <- system.time(
    draws_np <- nutpie_sample(mod_np, data = sdata, num_chains = 2, cores = 2,
                              num_draws = 500, seed = 11,
                              store_divergences = TRUE)
)
stanfit_np <- nutpie_draws_to_stanfit(draws_np, seed = 11)
bfit_np <- emptyfit; bfit_np$fit <- stanfit_np; bfit_np <- rename_pars(bfit_np)

## ---- comparisons ----
va <- variables(as_draws(bfit_cs)); vb <- variables(as_draws(bfit_np))
cat("CHECK variables identical (incl lp__):", identical(sort(va), sort(vb)),
    "| n_cs:", length(va), " n_np:", length(vb), "\n")

s_cs <- summarise_draws(as_draws(bfit_cs), "mean", "sd", "ess_bulk")
s_np <- summarise_draws(as_draws(bfit_np), "mean", "sd", "ess_bulk")
cmp <- inner_join(s_np, s_cs, by = "variable", suffix = c("_np", "_cs")) |>
    # drop constants (disc: sd=0 -> ESS=NA) and lp__ (samplers may differ by
    # an additive normalizing constant)
    filter(sd_np > 0 | sd_cs > 0, variable != "lp__") |>
    mutate(z = (mean_np - mean_cs) /
               sqrt(sd_np^2 / ess_bulk_np + sd_cs^2 / ess_bulk_cs))
cat("CHECK join covers all:", nrow(cmp) == length(vb) - 2,
    "| max |z|:", round(max(abs(cmp$z)), 2),
    "| frac |z|>3:", round(mean(abs(cmp$z) > 3), 4), "\n")

div_np <- sum(rstan::get_divergent_iterations(bfit_np$fit))
cat("CHECK get_divergent_iterations works, count:", div_np, "\n")

nd <- data.frame(x = c(-1, 0, 1))
pe_cs <- posterior_epred(bfit_cs, newdata = nd, re_formula = NA)
pe_np <- posterior_epred(bfit_np, newdata = nd, re_formula = NA)
dp <- max(abs(apply(pe_cs, c(2, 3), mean) - apply(pe_np, c(2, 3), mean)))
cat("CHECK posterior_epred runs both | max |diff| of mean category probs:",
    round(dp, 4), "(should be MC-noise, ~<0.02)\n")

cat("\nsummary() smoke test (nutpie-backed brmsfit):\n")
print(summary(bfit_np), digits = 2)

cat(sprintf("\nwall: cmdstanr %.1fs (500+500) vs nutpie %.1fs (400+500 default warmup)\n",
            t_cs["elapsed"], t_np["elapsed"]))
cat("DONE_TEST\n")

## ---- part 2: multivariate (mvbind) coverage -- the production model shape ----
set.seed(7)
dm <- data.frame(g = rep(1:50, each = 4), x = rnorm(200))
rem <- rnorm(50, 0, 0.8)
mk <- function(shift) cut(shift + 0.6 * dm$x + rem[dm$g] + rlogis(200),
                          c(-Inf, -1.5, -0.5, 0.5, 1.5, Inf), labels = FALSE)
dm$Q01 <- mk(0); dm$Q02 <- mk(0.4); dm$Q03 <- mk(-0.4)

fm <- bf(mvbind(Q01, Q02, Q03) ~ 1 + x + (1 | p | g),
         family = cumulative()) + set_rescor(FALSE)
scm <- make_stancode(fm, data = dm)
sdm <- make_standata(fm, data = dm); class(sdm) <- NULL
emptym <- brm(fm, data = dm, empty = TRUE)

mod_cs2 <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(as.character(scm)))
fit_cs2 <- mod_cs2$sample(data = sdm, chains = 2, parallel_chains = 2,
                          iter_warmup = 400, iter_sampling = 400, init = 0.1,
                          seed = 12, refresh = 0)
b_cs2 <- emptym
b_cs2$fit <- brms::read_csv_as_stanfit(fit_cs2$output_files())
b_cs2 <- rename_pars(b_cs2)

mod_np2 <- nutpie_compile_model(code = as.character(scm))
dr2 <- nutpie_sample(mod_np2, data = sdm, num_chains = 2, cores = 2,
                     num_draws = 400, seed = 12, store_divergences = TRUE)
b_np2 <- emptym
b_np2$fit <- nutpie_draws_to_stanfit(dr2, seed = 12)
b_np2 <- rename_pars(b_np2)

va2 <- variables(as_draws(b_cs2)); vb2 <- variables(as_draws(b_np2))
cat("CHECK-MV variables identical:", identical(sort(va2), sort(vb2)),
    "| n:", length(va2), "/", length(vb2), "\n")
s1 <- summarise_draws(as_draws(b_np2), "mean", "sd", "ess_bulk")
s2 <- summarise_draws(as_draws(b_cs2), "mean", "sd", "ess_bulk")
cmp2 <- inner_join(s1, s2, by = "variable", suffix = c("_np", "_cs")) |>
    filter(sd_np > 0 | sd_cs > 0, variable != "lp__") |>
    mutate(z = (mean_np - mean_cs) /
               sqrt(sd_np^2 / ess_bulk_np + sd_cs^2 / ess_bulk_cs))
if (anyNA(cmp2$z)) {
    cat("NOTE-MV z is NA for:",
        paste(head(cmp2$variable[is.na(cmp2$z)], 10), collapse = ", "), "\n")
    cmp2 <- filter(cmp2, !is.na(z))
}
pe1 <- posterior_epred(b_np2, newdata = data.frame(x = 0), re_formula = NA)
pe2 <- posterior_epred(b_cs2, newdata = data.frame(x = 0), re_formula = NA)
cat("CHECK-MV max |z|:", round(max(abs(cmp2$z)), 2),
    "| epred dims equal:", identical(dim(pe1), dim(pe2)),
    "| max epred diff:",
    round(max(abs(apply(pe1, c(2, 3), mean) - apply(pe2, c(2, 3), mean))), 4), "\n")
cat("CHECK-MV divergences count:",
    sum(rstan::get_divergent_iterations(b_np2$fit)), "\n")
cat("DONE_TEST_MV\n")
