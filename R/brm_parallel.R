brm_parallel <- function(args_shared, args_per_fit,
                         summarise_fun = NULL,
                         backend = options("brms.backend"),
                         ...) {


    n_fits <- length(args_per_fit)

    brms_argnames <- c("formula", "family", "data", "prior", "autocor",
                   "data2", "cov_ranef", "sample_prior", "stanvars",
                   "threads", "normalize")

    processed_args_per_fit <- list()
    model_codes <- array(NA_character_, n_fits)


    for(i in 1:n_fits) {
        for(a in brms_argnames) {
            if(a %in% names(args_shared) && a %in% names(args_per_fit)) {
                stop(paste0("Both args_shared and args_per_fit[[", i, "]] provide '", a, "'"))
            }
        }
        brms_args <- c(args_shared[intersect(names(args_shared), brms_argnames)],
                               args_per_fit[[i]][intersect(names(args_per_fit[[i]]), brms_argnames)])

        stan_args <- c(args_shared[setdiff(names(args_shared), brms_argnames)],
                               args_per_fit[[i]][setdiff(names(args_per_fit[[i]]), brms_argnames)])


        stancode <- do.call(brms::make_stancode, args = brms_args)
        standata <- do.call(brms::make_standata, args = brms_args)
        emptyfit <- do.call(brms::brm, args = c(brms_args, list(empty = TRUE)))

        class(standata) <- NULL
        stan_args$data <- standata

        same_model_indices <- which(model_codes == as.character(stancode))
        if(length(same_model_indices) > 0) {
            stan_args$model <- processed_args_per_fit[[same_model_indices[1]]]$model
        } else {
            if(backend == "rstan") {
                stan_args$model <- rstan::stan_model(model_code = stancode)
            } else if(backend == "cmdstanr"){
                stan_args$model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
            } else {
                stop("Unrecognized backend")
            }
        }

        model_codes[i] <- stancode

        stan_args$summarise_fun_args <- c(list(emptyfit = emptyfit), stan_args$summarise_fun_args)

        processed_args_per_fit[[i]] <- stan_args
    }

    brms_summarise_fun <- function(fit, emptyfit, ...) {
        brmsfit <- emptyfit
        brmsfit$fit <- fit
        brmsfit <- rename_pars(brmsfit)
        if(is.null(summarise_fun)) {
            brmsfit
        } else {
            summarise_fun(brmsfit, ...)
        }
    }

    sampling_parallel(args_shared = list(), args_per_fit = processed_args_per_fit,
                   convert_cmdstan_fits_to_rstan = TRUE,
                   summarise_fun = brms_summarise_fun,
                   ...)
}
