#' @param summarise_fun a function to process each fit. This function is run in parallel.
#'   Note that this function must be runnable in new RStudio sessions. You may
#'   use the  `summarise_fun_dependencies` parameter to ensure libraries are loaded (but you should ideally just be explicit and use `::`).
#' @param summarise_fun_dependencies a list of package names that need to be loaded for
#'   `summarise_fun` to run. IMPORTANT: when developing packages, you need to install
#'   the latest version, `devtools::load_all()` won't be enough
#'    (the packages are loaded from the default library)
#' @param cache_dir if not NULL, fits will be cached in this directory
#' @param nutpie_args extra arguments for the nutpieR backend, e.g.
#'   `list(adaptation = "low_rank")` or `list(num_warmup = 800)`. Ignored for the
#'   rstan/cmdstanr backends. NOTE: `warmup`/`iter_warmup` is deliberately NOT translated to
#'   `num_warmup` -- nuts-rs adapts within its own budget (400 draws for
#'   `adaptation = "diag"`, 800 for `"low_rank"`), and forcing cmdstan's warmup
#'   would discard the adaptation difference this backend exists to exploit.
#'   Set `nutpie_args$num_warmup` to override.
#' @param extra_cache_key optional object folded into the SUMMARY cache key only
#'   (not the fit key -- a fit does not depend on how it is summarised).
#'   `brm_parallel` uses it to pass the source of `summarise_fun`, so that
#'   editing the summary function invalidates cached summaries but not the
#'   expensive fits behind them.
#' @param nutpie_converter function turning nutpieR draws into a `stanfit`.
#'   Defaults to `nutpie_draws_to_stanfit()` from R/nutpie_to_stanfit.R, which
#'   must be sourced when a nutpie model is used.
#' @return A list of length `length(data)` containing the result of applying
#'   `summarise_fun` to each fit.
sampling_parallel <- function(args_shared, args_per_fit,
                            cores_per_fit = NULL,
                            convert_cmdstan_fits_to_rstan = FALSE,
                            fits_in_parallel = NULL,
                            summarise_fun = NULL,
                            cache_fits = FALSE,
                            cache_summaries = FALSE,
                            cache_dir = NULL,
                            extra_cache_key = NULL,
                            nutpie_args = list(),
                            nutpie_converter = NULL,
                            future.chunk.size = 1,
                            future.globals = TRUE,
                            future.stdout = FALSE
                            ) {

  if(!is.list(args_shared)) {
    stop("args_shared must be a list")
  }
  if(!is.list(args_per_fit) || length(args_per_fit) <= 0) {
    stop("args_per_fit must be a non-empty list")
  }

  if(!is.null(cache_dir) && !dir.exists(cache_dir)) {
    stop(paste0("Cache dir '", cache_dir,"'  does not exist"))
  }

  if((cache_fits || cache_summaries) && is.null(cache_dir)){
    stop("Caching turned on but cache_dir not given")
  }

  if(cache_summaries && is.null(summarise_fun)) {
    stop("cache_summaries can only be used if summarise_fun is not null")
  }

  total_cores <- future::nbrOfWorkers()

  n_fits <- length(args_per_fit)

  if("cores" %in% names(args_shared) || "num_cores" %in% names(args_shared)) {
    stop("args_shared must not specify cores or num_cores")
  }

  uses_rstan <- FALSE
  uses_cmdstan <- FALSE
  uses_nutpie <- FALSE
  model_in_shared_args <- FALSE
  data_in_shared_args <- "data" %in% names(args_shared)

  if("model" %in% names(args_shared)) {
    if(inherits(args_shared$model, "stanmodel")) {
      uses_rstan <- TRUE
    } else if(inherits(args_shared$model, "CmdStanModel")) {
      uses_cmdstan <- TRUE
    } else if(inherits(args_shared$model, "nutpie_model")) {
      uses_nutpie <- TRUE
    } else {
      stop("Model in shared args is not of class 'stanmodel', 'CmdStanModel' or 'nutpie_model'")
    }
    model_in_shared_args <- TRUE
  }

  for(i in 1:n_fits) {
    if(!is.list(args_per_fit[[i]])) {
      stop("All elements of args_per_fit have to be lists")
    }

    if(length(intersect(names(args_shared), names(args_per_fit[[i]]))) > 0) {
      stop(paste0("No parameters provided in args_per_fit can be given in args_shared.\n
                 Found intersection at index ", i, "."))
    }

    if("model" %in% names(args_per_fit[[i]])) {
      if(inherits(args_per_fit[[i]]$model, "stanmodel")) {
        uses_rstan <- TRUE
      } else if(inherits(args_per_fit[[i]]$model, "CmdStanModel")) {
        uses_cmdstan <- TRUE
      } else if(inherits(args_per_fit[[i]]$model, "nutpie_model")) {
        uses_nutpie <- TRUE
      } else {
        stop(paste0("Model for fit id ", i,
                    " is not of class 'stanmodel', 'CmdStanModel' or 'nutpie_model'"))
      }
    } else if(!model_in_shared_args) {
      stop(paste0("No model argument in shared_args and fit id ", i, " does not provide model"))
    }

    if(!data_in_shared_args && !("data" %in% names(args_per_fit[[i]]))) {
      stop(paste0("No data argument in shared_args and fit id ", i, " does not provide data"))
    }


    if("cores" %in% names(args_per_fit[[i]]) || "num_cores" %in% names(args_per_fit[[i]])) {
      stop(paste0("args_per_fit[[", i, "]] must not specify cores or num_cores"))
    }
  }

  ## ---- nutpieR backend setup (all resolved in the master process) ----
  nutpie_tag <- NULL
  nutpie_metric <- NULL
  if(uses_nutpie) {
    if(!requireNamespace("nutpieR", quietly = TRUE)) {
      stop("A nutpie_model was supplied but the nutpieR package is not available")
    }
    if(is.null(nutpie_converter)) {
      nutpie_converter <- get0("nutpie_draws_to_stanfit", mode = "function")
      if(is.null(nutpie_converter)) {
        stop("nutpie_draws_to_stanfit() not found: source(here::here('R', 'nutpie_to_stanfit.R'))",
             " or pass nutpie_converter explicitly")
      }
    }
    if(!is.list(nutpie_args)) {
      stop("nutpie_args must be a list")
    }
    unknown <- setdiff(names(nutpie_args), names(formals(nutpieR::nutpie_sample)))
    if(length(unknown) > 0) {
      stop(paste0("nutpie_args contains arguments nutpie_sample() does not take: ",
                  paste(unknown, collapse = ", ")))
    }
    adaptation <- if(is.null(nutpie_args$adaptation)) "diag" else nutpie_args$adaptation
    adaptation <- match.arg(adaptation, c("diag", "low_rank", "low-rank"))
    # cache keys and the stanfit metric label must distinguish the two
    # adaptations: they are different samplers on the same model + data
    nutpie_tag <- if(identical(adaptation, "diag")) "nutpie-diag" else "nutpie-lowrank"
    nutpie_metric <- if(identical(adaptation, "diag")) "diag_e" else "low_rank"
    default_warmup <- if(identical(adaptation, "diag")) 400L else 800L
    message("sampling_parallel: nutpieR backend, adaptation = '", adaptation,
            "', num_warmup = ",
            if(is.null(nutpie_args$num_warmup)) {
              paste0(default_warmup, " (nuts-rs default; warmup/iter_warmup is not copied over)")
            } else {
              nutpie_args$num_warmup
            },
            ", cache tag = '", nutpie_tag, "'")
  }

  if(is.null(fits_in_parallel)) {
    if(2 * n_fits <= total_cores) {
      fits_in_parallel <- n_fits
    } else {
      fits_in_parallel <- min(c(total_cores, n_fits))
    }
  }

  if(is.null(cores_per_fit)) {
    if(2 * n_fits <= total_cores) {
      cores_per_fit <- floor(total_cores / n_fits)
    } else {
      cores_per_fit <- 1
    }
  }

  fit_fun <- function(args, args_shared, summarise_fun,
                      convert_cmdstan_fits_to_rstan,
                      cores_per_fit,
                      cache_dir, cache_fits, cache_summaries,
                      cmdstan_fit_dir, summary_key_extra,
                      nutpie_args, nutpie_converter, nutpie_tag, nutpie_metric, #nutpie adds 5 parameters: summary_key_extra - nutpie_metric
                      progressor) {
    all_args <- c(args_shared, args)
    all_args$cores <- cores_per_fit

    model <- all_args$model
    all_args$model <- NULL

    summarise_fun_args <- all_args$summarise_fun_args
    all_args$summarise_fun_args <- NULL

    cache_hashes <- all_args$.cache_hashes
    all_args$.cache_hashes <- NULL



    if(inherits(model, "stanmodel")) {
      model_code <- model@model_code
      backend <- "rstan"
    } else if(inherits(model, "CmdStanModel")) {
      model_code <- model$code()
      backend <- "cmdstanr"
    } else if(inherits(model, "nutpie_model")) {
      # the staged .stan file in nutpie's content-hashed compile cache is an
      # exact, stable copy of the code the shared library was built from
      model_code <- paste(readLines(model$staged_source), collapse = "\n")
      backend <- "nutpier"
    } else {
      stop("Invalid model")
    }
    is_rstan <- backend == "rstan"


    if(cache_fits || cache_summaries) {
      if(is.null(cache_hashes)) {
        stop("Cache hashes missing -- they are computed in the master process")
      }
      code_hash <- cache_hashes[["code"]]
      data_hash <- cache_hashes[["data"]]
      # nutpie fits are stanfits too, so without a tag of their own an untagged
      # key would let the backends silently hand each other's fits back
      cache_tag <- if(backend == "nutpier") paste0(nutpie_tag, "_") else ""
    }

    summary_cached <-  FALSE
    if(!is.null(summarise_fun) && cache_summaries) {
      summary_cache_file <- paste0(cache_dir, "/summary_", cache_tag, code_hash, "_",
                                   data_hash, "_", summary_key_extra, ".rds")
      if(file.exists(summary_cache_file)) {
        result <- readRDS(summary_cache_file)
        summary_cached <- TRUE
      }
    }

    if(!summary_cached) {
      fit_cached <- FALSE
      if(cache_fits) {
        fit_cache_file <- paste0(cache_dir, "/fit_", cache_tag, code_hash, "_", data_hash, ".rds")
        if(file.exists(fit_cache_file)) {
          fit_from_file <- readRDS(fit_cache_file)
          if((is_rstan && inherits(fit_from_file, "stanfit"))
             || (backend == "nutpier" && inherits(fit_from_file, "stanfit"))
             || (backend == "cmdstanr" && !convert_cmdstan_fits_to_rstan && inherits(fit_from_file, "CmdStanMCMC"))
             || (backend == "cmdstanr" && convert_cmdstan_fits_to_rstan && inherits(fit_from_file, "stanfit"))
             ) {
            fit <- fit_from_file
            fit_cached <- TRUE
          }
        }
      }


      if(!fit_cached) {
        if(backend == "rstan") { #rewritten from if(inherits(model, "stanmodel"))
          all_args_ordered <- c(list(model), all_args)
          fit <- do.call(rstan::sampling, args = all_args_ordered)
          if(!is.null(cache_dir)) {
            saveRDS(fit, fit_cache_file)
          }
        } else if(backend == "nutpier") {
          # brms/cmdstanr-shaped sampling arguments -> nutpie_sample() arguments
          np <- list(data = all_args$data,
                     cores = all_args$cores,
                     store_divergences = TRUE,  # the converter needs the NUTS stats
                     progress = "none",
                     refresh = 0)

          if(!is.null(all_args$chains)) {
            np$num_chains <- all_args$chains
          }
          if(!is.null(all_args$iter_sampling)) {
            np$num_draws <- all_args$iter_sampling
          } else if(!is.null(all_args$iter)) {
            if(!is.null(all_args$warmup)) {
              np$num_draws <- all_args$iter - all_args$warmup
            } else {
              np$num_draws <- all_args$iter / 2
            }
          }
          # warmup/iter_warmup is intentionally NOT mapped to num_warmup; see
          # the nutpie_args documentation above.
          if(!is.null(all_args$adapt_delta)) {
            np$target_accept <- all_args$adapt_delta
          }
          if(!is.null(all_args$control$adapt_delta)) {
            np$target_accept <- all_args$control$adapt_delta
          }
          if(!is.null(all_args$max_treedepth)) {
            np$max_treedepth <- all_args$max_treedepth
          }
          if(!is.null(all_args$control$max_treedepth)) {
            np$max_treedepth <- all_args$control$max_treedepth
          }
          # `init` means the same thing in both backends: a scalar x starts each
          # chain from Uniform(-x, x) on the unconstrained scale
          if(!is.null(all_args$init)) {
            np$init <- all_args$init
          }
          if(!is.null(all_args$seed)) {
            np$seed <- all_args$seed
          }
          # anything else the caller passed that nutpie_sample() understands;
          # cmdstanr-only arguments (parallel_chains, iter_warmup, ...) are dropped
          passthrough <- setdiff(intersect(names(all_args),
                                           names(formals(nutpieR::nutpie_sample))),
                                 c(names(np), "model"))
          for(nm in passthrough) {
            np[[nm]] <- all_args[[nm]]
          }
          np[names(nutpie_args)] <- nutpie_args

          draws <- do.call(nutpieR::nutpie_sample, args = c(list(model), np))
          fit <- nutpie_converter(
            draws,
            seed = if(is.null(np$seed)) NA_integer_ else np$seed,
            adapt_delta = if(is.null(np$target_accept)) NA_real_ else np$target_accept,
            max_treedepth = if(is.null(np$max_treedepth)) NA_integer_ else np$max_treedepth,
            metric = nutpie_metric)
          if(!is.null(cache_dir) && cache_fits) {
            saveRDS(fit, fit_cache_file)
          }
        } else {
          translated_args <- list()
          for(old in names(all_args)) {
            if(old == "chains") {
              translated_args$num_chains = all_args$chains
            } else if(old == "cores") {
              translated_args$parallel_chains = all_args$cores
            } else if(old == "control") {
              if(!is.null(all_args$control$adapt_delta)) {
                translated_args$adapt_delta = all_args$control$adapt_delta
              }
              if(!is.null(all_args$control$max_treedepth)) {
                translated_args$max_depth = all_args$control$max_treedepth
              }
            } else if(old == "warmup") {
              translated_args$iter_warmup = all_args$warmup
            } else if(old == "iter") {
              if("warmup" %in% names(all_args)) {
                  translated_args$iter_sampling = all_args$iter - all_args$warmup
              } else {
                  translated_args$iter_warmup = all_args$iter / 2
                  translated_args$iter_sampling = all_args$iter/ 2
              }
            } else if(old == "warmup") {
                translated_args$iter_warmup = all_args$warmup
            } else {
              translated_args[[old]] = all_args[[old]]
            }
          }
          fit <- do.call(model$sample, args = translated_args)
          if(convert_cmdstan_fits_to_rstan) {
            fit <- brms::read_csv_as_stanfit(fit$output_files())
            if(!is.null(cache_dir) && cache_fits) {
              saveRDS(fit, fit_cache_file)
            }
          } else {
            fit$save_output_files(cmdstan_fit_dir)
            if(!is.null(cache_dir) && cache_fits) {
              fit$save_object(fit_cache_file)
            }
          }
        }
      } # End - if(!fit_cached)

      if(!is.null(summarise_fun)) {
        result <- do.call(summarise_fun, args = c(list(fit), summarise_fun_args))
        if(!is.null(cache_dir) && cache_summaries) {
          saveRDS(result, summary_cache_file)
        }
      } else {
        result <- fit
      }

    } # End - if(!summary_cached)

    if(!is.null(progressor)) {
        progressor()
    }
    result
  }

  if(requireNamespace("progressr", quietly = TRUE)) {
      progressor <- progressr::progressor(n_fits)
  } else {
      progressor <- NULL
  }

  ## ---- cache keys ----
  ## These MUST be computed here, on the un-shared objects. rlang::hash() of a
  ## mori::share()d object is not reproducible: the shared wrapper serialises
  ## differently on every share (while identical() still returns TRUE), so
  ## hashing inside the worker handed every fit a fresh key and the cache never
  ## hit - silently, since a miss just refits.
  ##
  ## The key is (model code, data), as originally intended. It does NOT cover
  ## the sampler arguments: two runs differing only in iter/warmup/adapt_delta/
  ## seed share a key and the first one wins, so use a separate cache_dir per
  ## sampler configuration.
  summary_key_extra <- substr(rlang::hash(extra_cache_key), 1, 8)
  if(cache_fits || cache_summaries) {
    model_code_of <- function(model) {
      if(inherits(model, "stanmodel")) {
        model@model_code
      } else if(inherits(model, "CmdStanModel")) {
        model$code()
      } else {
        paste(readLines(model$staged_source), collapse = "\n")
      }
    }
    for(i in 1:n_fits) {
      all_args_i <- c(args_shared, args_per_fit[[i]])
      args_per_fit[[i]]$.cache_hashes <- c(
        code = rlang::hash(model_code_of(all_args_i$model)),
        data = rlang::hash(all_args_i$data))
    }
  }

  shared_args_per_fit <- mori::share(args_per_fit)
  shared_args_shared <- mori::share(args_shared)

  results <- future.apply::future_lapply(
      X = shared_args_per_fit,
      FUN = fit_fun,
      args_shared = shared_args_shared,
      summarise_fun = summarise_fun,
      convert_cmdstan_fits_to_rstan = convert_cmdstan_fits_to_rstan,
      cores_per_fit = cores_per_fit,
      cache_dir = cache_dir,
      cmdstan_fit_dir = tempdir(),
      cache_fits = cache_fits,
      cache_summaries = cache_summaries,
      summary_key_extra = summary_key_extra,
      nutpie_args = nutpie_args,
      nutpie_converter = nutpie_converter,
      nutpie_tag = nutpie_tag,
      nutpie_metric = nutpie_metric,
      progressor = progressor,
      future.seed = TRUE,
      future.chunk.size = future.chunk.size,
      future.globals = future.globals,
      future.stdout = future.stdout
  )

  results
}

