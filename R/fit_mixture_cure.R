#' Fit a Weibull mixture cure model (mixcuref backend)
#'
#' A thin wrapper around \code{\link{mixcure.penal.est}} that accepts separate
#' incidence and latency formulas (as used by the meta-layer package) and
#' constructs a single survival regression formula for the mixcuref estimator.
#'
#' The current mixcuref engine uses a single covariate design matrix for both the
#' cure (incidence) and survival (latency) components. Therefore, this wrapper
#' takes the union of covariates appearing on the right-hand sides of
#' \code{incidence} and \code{latency}.
#'
#' @param incidence Formula for the incidence (cure) part. Only the right-hand
#'   side is used; any left-hand side is ignored.
#' @param latency Formula for the latency (susceptible survival) part. Only the
#'   right-hand side is used; any left-hand side is ignored.
#' @param time Numeric vector of follow-up times.
#' @param status Event indicator (1 = event, 0 = censored).
#' @param data Data frame containing covariates referenced by the formulas.
#' @param control List of engine-specific control options. Recognized elements:
#'   \itemize{
#'     \item \code{pl}: logical; use Firth-type penalized likelihood (default
#'       \code{FALSE}).
#'     \item \code{iterlim}: integer; iteration limit passed to
#'       \code{mixcure.penal.est} (default \code{200}).
#'     \item \code{init}: numeric vector of initial parameter values. If not
#'       supplied, a default vector of zeros is used with Weibull shape
#'       initialized at 1.
#'   }
#'
#' @return An object of class \code{"mixcure"} (a list) returned by
#'   \code{\link{mixcure.penal.est}}.
#' @export
fit_mixture_cure <- function(incidence, latency, time, status, data, control = list()) {

  stopifnot(inherits(incidence, "formula"))
  stopifnot(inherits(latency, "formula"))
  stopifnot(is.data.frame(data))

  if (length(time) != nrow(data)) {
    stop("'time' must have length equal to nrow(data).", call. = FALSE)
  }
  if (length(status) != nrow(data)) {
    stop("'status' must have length equal to nrow(data).", call. = FALSE)
  }

  # Coerce status to 0/1 integer where possible
  status <- as.integer(status)
  if (any(!status %in% c(0L, 1L), na.rm = TRUE)) {
    stop("'status' must be coded 0/1 (0=censored, 1=event).", call. = FALSE)
  }

  # Extract RHS covariates; ignore any LHS in incidence/latency.
  # Use terms(..., data=...) so '.' can expand when provided.
  rhs_inc <- attr(stats::terms(incidence, data = data), "term.labels")
  rhs_lat <- attr(stats::terms(latency, data = data), "term.labels")
  rhs_all <- unique(c(rhs_inc, rhs_lat))

  rhs <- if (length(rhs_all) == 0) "1" else paste(rhs_all, collapse = " + ")

  # Build a single Surv(...) ~ RHS formula for the engine
  surv_formula <- stats::as.formula(paste0("survival::Surv(time, status) ~ ", rhs))

  # Build a local data frame that contains time/status as named columns
  dat <- data
  dat$time <- time
  dat$status <- status

  # Control defaults
  pl <- isTRUE(control$pl)
  iterlim <- if (!is.null(control$iterlim)) as.integer(control$iterlim) else 200L

  # Compute default init if not provided
  init <- control$init
  if (is.null(init)) {
    mf <- stats::model.frame(surv_formula, data = dat, na.action = stats::na.omit)
    X <- stats::model.matrix(surv_formula, data = mf)
    k <- ncol(X)
    init <- c(rep(0, 2 * k), 1)  # last element is Weibull shape > 0
  }

  mixcure.penal.est(
    formula = surv_formula,
    data = dat,
    init = init,
    pl = pl,
    iterlim = iterlim
  )
}
