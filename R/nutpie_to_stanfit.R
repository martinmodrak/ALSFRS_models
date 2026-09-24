# Convert nutpieR draws into an rstan `stanfit`, so the standard brm_parallel
# downstream (emptyfit$fit <- fit; brms::rename_pars(); posterior_epred();
# rstan::get_divergent_iterations()) works unchanged with a nutpie backend.
#
# Modeled on brms::read_csv_as_stanfit (brms 2.23.1), which parses CmdStan CSVs
# and then assembles the stanfit from plain R objects; here the same structure
# is assembled directly from nutpie's in-memory results, skipping the CSV
# round-trip. Like the brms original, the resulting stanfit carries a null
# stanmodel (no compiled DSO), so draws-based methods (epred, loo, pp_check)
# work while log_prob-based ones (bridgesampling) do not --> identical to the
# current cmdstanr conversion path.
#
# draws: posterior::draws_array returned by nutpieR::nutpie_sample()
#        (post-warmup only; sampler stats are pulled via nutpie_nuts_params,
#        lp__ is recovered from nutpie_diagnostics()$logp).
#
# metric: label only, recorded in the stanfit's sampler metadata so a
#         low-rank-adapted fit is distinguishable from a diagonal one after the
#         fact ("diag_e" for adaptation = "diag", "low_rank" for "low_rank").
nutpie_draws_to_stanfit <- function(draws,
                                    model_name = "nutpie_model",
                                    seed = NA_integer_,
                                    adapt_delta = NA_real_,
                                    max_treedepth = NA_integer_,
                                    metric = "diag_e") {
    stopifnot(inherits(draws, "draws_array"))

    n_iter <- posterior::niterations(draws)
    n_chains <- posterior::nchains(draws)

    ## draws as one data.frame with bracket-style names, then split by chain
    dd <- as.data.frame(posterior::as_draws_df(draws))
    chain_ids <- dd$.chain
    dd[c(".chain", ".iteration", ".draw")] <- NULL

    ## lp__ is not a draws variable in nutpie, but the per-draw log density is
    ## available from the diagnostics; append it as the final column (CmdStan
    ## convention). Note: may differ from CmdStan's lp__ by an additive
    ## constant (dropped normalizing terms), so compare across samplers only
    ## up to a constant.
    dg <- nutpieR::nutpie_diagnostics(draws)
    if (!is.null(dg$logp) && !is.null(dg$chain) && !is.null(dg$draw)) {
        keep <- if (!is.null(dg$tuning)) !dg$tuning else rep(TRUE, length(dg$logp))
        lp <- dg$logp[keep][order(dg$chain[keep], dg$draw[keep])]
        if (length(lp) == nrow(dd)) {
            dd$lp__ <- lp
        }
    }

    ## sampler diagnostics: long (Chain, Iteration, Parameter, Value) -> wide,
    ## in the column order rstan expects; missing ones are simply omitted
    nuts <- nutpieR::nutpie_nuts_params(draws)
    rstan_diagn_order <- c("accept_stat__", "treedepth__", "stepsize__",
                           "divergent__", "n_leapfrog__", "energy__")
    diagnostics <- lapply(seq_len(n_chains), function(ch) {
        x <- nuts[nuts$Chain == ch, , drop = FALSE]
        cols <- lapply(rstan_diagn_order, function(p) {
            xi <- x[x$Parameter == p, , drop = FALSE]
            if (nrow(xi) == 0) return(NULL)
            xi$Value[order(xi$Iteration)]
        })
        names(cols) <- rstan_diagn_order
        cols <- cols[!vapply(cols, is.null, logical(1))]
        as.data.frame(cols, check.names = FALSE)
    })

    ## parameter names and dimensions, parsed from the flat names
    fnames <- colnames(dd)
    base <- sub("\\[.*\\]$", "", fnames)
    model_pars <- unique(base)
    par_dims <- lapply(model_pars, function(p) {
        idx <- fnames[base == p]
        if (!grepl("[", idx[1], fixed = TRUE)) return(integer(0))
        ind <- sub("^.*\\[", "", sub("\\]$", "", idx))
        im <- do.call(rbind, lapply(strsplit(ind, ","), as.integer))
        as.integer(apply(im, 2, max))
    })
    names(par_dims) <- model_pars

    ## per-chain sample data.frames with the attributes rstan reads
    samples <- split(dd, chain_ids)
    names(samples) <- NULL
    step_sizes <- rep(NA_real_, n_chains)
    for (i in seq_len(n_chains)) {
        rownames(samples[[i]]) <- seq_len(nrow(samples[[i]]))
        attr(samples[[i]], "sampler_params") <- diagnostics[[i]]
        rownames(attr(samples[[i]], "sampler_params")) <-
            seq_len(nrow(diagnostics[[i]]))
        attr(samples[[i]], "adaptation_info") <- character(0)
        attr(samples[[i]], "args") <- list(sampler_t = paste0("NUTS(", metric, ")"),
                                           chain_id = i)
        m <- colMeans(samples[[i]])
        if ("lp__" %in% names(m)) {
            attr(samples[[i]], "mean_pars") <- m[-length(m)]
            attr(samples[[i]], "mean_lp__") <- m[["lp__"]]
        } else {
            attr(samples[[i]], "mean_pars") <- m
            attr(samples[[i]], "mean_lp__") <- NA_real_
        }
        if ("stepsize__" %in% names(diagnostics[[i]])) {
            step_sizes[i] <- utils::tail(diagnostics[[i]]$stepsize__, 1)
        }
    }

    sim <- list(
        samples = samples,
        iter = n_iter,
        thin = 1L,
        warmup = 0L,
        chains = n_chains,
        n_save = rep(n_iter, n_chains),
        warmup2 = rep(0L, n_chains),
        permutation = lapply(seq_len(n_chains), function(i) sample.int(n_iter)),
        pars_oi = model_pars,
        dims_oi = par_dims,
        fnames_oi = fnames,
        n_flatnames = length(fnames)
    )

    sargs <- list(
        stan_version_major = "2", stan_version_minor = "0",
        stan_version_patch = "0",
        model = model_name, method = "sample",
        iter = n_iter, warmup = 0L, save_warmup = 0L, thin = 1L,
        algorithm = "hmc", engine = "nuts", metric = metric,
        stepsize = NA_real_, chain_id = NA_integer_,
        seed = as.character(seed),
        sampler_t = paste0("NUTS(", metric, ")"),
        # rstan looks in @stan_args$control for some metadata
        control = list(adapt_delta = adapt_delta,
                       max_treedepth = max_treedepth)
    )
    sargs_rep <- replicate(n_chains, sargs, simplify = FALSE)
    for (i in seq_len(n_chains)) {
        sargs_rep[[i]]$chain_id <- i
        sargs_rep[[i]]$stepsize <- step_sizes[i]
    }

    # rstan's S4 classes must be registered for new("stanmodel")/new("stanfit")
    loadNamespace("rstan")
    cxxdso_class <- "cxxdso"
    attr(cxxdso_class, "package") <- "rstan"
    null_dso <- methods::new(
        cxxdso_class, sig = list(character(0)), dso_saved = FALSE,
        dso_filename = character(0), modulename = character(0),
        system = R.version$system, cxxflags = character(0),
        .CXXDSOMISC = new.env(parent = emptyenv())
    )
    null_sm <- methods::new(
        "stanmodel", model_name = model_name, model_code = character(0),
        model_cpp = list(), dso = null_dso
    )

    methods::new(
        "stanfit",
        model_name = model_name,
        model_pars = model_pars,
        par_dims = par_dims,
        mode = 0L,
        sim = sim,
        inits = list(),
        stan_args = sargs_rep,
        stanmodel = null_sm,
        date = format(Sys.time(), "%a %b %d %X %Y"),
        .MISC = new.env(parent = emptyenv())
    )
}
