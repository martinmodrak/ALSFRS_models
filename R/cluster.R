setup_cluster <- function(max_local_workers = Inf, max_server_workers = 140, force = FALSE) {
  is_local <- parallelly::availableCores() < 80
  if(is_local) {
    n_workers <- min(parallelly::availableCores(), max_local_workers)
    # future::plan(multisession, workers = n_workers)
  } else {
    n_workers <- min(parallelly::availableCores(), max_server_workers)
    # mirai::daemons(n_workers, dispatcher = FALSE, force = force)
    # future::plan(future.mirai::mirai_cluster)
  }
  future::plan(future::multisession, workers = n_workers)
}
