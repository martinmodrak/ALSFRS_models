#' @param summarise_fun a function to process each fit. This function is run in parallel.
#'   Note that this function must be runnable in new RStudio sessions. You may
#'   use the  `summarise_fun_dependencies` parameter to ensure libraries are loaded (but you should ideally just be explicit and use `::`).
#' @param summarise_fun_dependencies a list of package names that need to be loaded for
#'   `summarise_fun` to run. IMPORTANT: when developing packages, you need to install
#'   the latest version, `devtools::load_all()` won't be enough
#'    (the packages are loaded from the default library)
#' @param cache_dir if not NULL, fits will be cached in this directory
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
                            future.chunk.size = 1,
                            future.globals = NULL
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
  model_in_shared_args <- FALSE
  data_in_shared_args <- "data" %in% names(args_shared)

  if("model" %in% names(args_shared)) {
    if(inherits(args_shared$model, "stanmodel")) {
      uses_rstan <- TRUE
    } else if(inherits(args_shared$model, "CmdStanModel")) {
      uses_cmdstan <- TRUE
    } else {
      stop("Model in shared args is not of class 'stanmodel' or 'CmdStanModel'")
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
      } else {
        stop(paste0("Model for fit id ", i," is not of class 'stanmodel' or 'CmdStanModel'"))
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
                      cmdstan_fit_dir,
                      progressor) {
    all_args <- c(args_shared, args)
    all_args$cores <- cores_per_fit

    model <- all_args$model
    all_args$model <- NULL

    summarise_fun_args <- all_args$summarise_fun_args
    all_args$summarise_fun_args <- NULL



    if(inherits(model, "stanmodel")) {
      model_code <- model@model_code
      is_rstan <- TRUE
    } else if(inherits(model, "CmdStanModel")) {
      model_code <- model$code()
      is_rstan <- FALSE
    } else {
      stop("Invalid model")
    }




    if(cache_fits || cache_summaries) {
      data <- all_args$data
      code_hash <- rlang::hash(model_code)
      data_hash <- rlang::hash(data)
    }

    summary_cached <-  FALSE
    if(!is.null(summarise_fun) && cache_summaries) {
      summary_cache_file <- paste0(cache_dir, "/summary_", code_hash, "_", data_hash, ".rds")
      if(file.exists(summary_cache_file)) {
        result <- readRDS(summary_cache_file)
        summary_cached <- TRUE
      }
    }

    if(!summary_cached) {
      fit_cached <- FALSE
      if(cache_fits) {
        fit_cache_file <- paste0(cache_dir, "/fit_", code_hash, "_", data_hash, ".rds")
        if(file.exists(fit_cache_file)) {
          fit_from_file <- readRDS(fit_cache_file)
          if((is_rstan && inherits(fit_from_file, "stanfit"))
             || (!is_rstan && !convert_cmdstan_fits_to_rstan && inherits(fit_from_file, "CmdStanMCMC"))
             || (!is_rstan && convert_cmdstan_fits_to_rstan && inherits(fit_from_file, "stanfit"))
             ) {
            fit <- fit_from_file
            fit_cached <- TRUE
          }
        }
      }


      if(!fit_cached) {
        if(inherits(model,"stanmodel")) {
          all_args_ordered <- c(list(model), all_args)
          fit <- do.call(rstan::sampling, args = all_args_ordered)
          if(!is.null(cache_dir)) {
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

  results <- future.apply::future_lapply(
      X = args_per_fit,
      FUN = fit_fun,
      args_shared = args_shared,
      summarise_fun = summarise_fun,
      convert_cmdstan_fits_to_rstan = convert_cmdstan_fits_to_rstan,
      cores_per_fit = cores_per_fit,
      cache_dir = cache_dir,
      cmdstan_fit_dir = tempdir(),
      cache_fits = cache_fits,
      cache_summaries = cache_summaries,
      future.seed = TRUE,
      future.chunk.size = future.chunk.size,
      future.globals = future.globals,
      progressor = progressor
  )

  results
}

