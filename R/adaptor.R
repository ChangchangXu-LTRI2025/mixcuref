#' Adapter for mixcure.meta: fit from a cure_spec
#'
#' @param spec A mc_spec object
#' @param control A list of engine-specific controls
#' @return An engine-specific fitted object
#' @export
mixcuref_fit <- function(spec, control = list()) {
  dat <- spec$data

  # Replace `fit_mixture_cure` with YOUR existing exported function name
  mixcuref::fit_mixture_cure(
    incidence = spec$incidence,
    latency   = spec$latency,
    time      = dat[[spec$time]],
    status    = dat[[spec$status]],
    data      = dat,
    control   = control
  )
}
