#' Compute k and b Burst Parameters from Vector Inputs
#'
#' This function computes the burst size parameter \code{k} and burst frequency
#' parameter \code{b} from vector inputs of negative binomial estimates
#' (\code{r_est}, \code{p_est}). Error (RMSE) values are optional; if omitted, 
#' the function returns \code{NA} for RMSE outputs.
#'
#' The transformations applied are:
#' \itemize{
#'   \item \code{k_est  = r_est}
#'   \item \code{k_rmse = r_rmse}
#'   \item \code{b_est  = (1 - p_est) / p_est}
#'   \item \code{b_rmse = p_rmse / (p_est^2)}
#' }
#'
#' @param r_est Numeric vector of NB size estimates.
#' @param p_est Numeric vector of NB probability estimates.
#' @param r_rmse Optional numeric vector of RMSE values for \code{r_est}.
#' @param p_rmse Optional numeric vector of RMSE values for \code{p_est}.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{k_est}
#'     \item \code{k_rmse}
#'     \item \code{b_est}
#'     \item \code{b_rmse}
#'   }
#'
#' @export
kb_calc <- function(r_est,
                    p_est,
                    r_rmse = NULL,
                    p_rmse = NULL) {
  
  # ---- Validate required inputs ----
  if (!is.numeric(r_est)) stop("r_est must be numeric.")
  if (!is.numeric(p_est)) stop("p_est must be numeric.")
  
  if (length(r_est) != length(p_est)) {
    stop("r_est and p_est must have the same length.")
  }
  
  n <- length(r_est)
  
  # ---- Validate optional RMSE inputs ----
  if (!is.null(r_rmse)) {
    if (!is.numeric(r_rmse)) stop("r_rmse must be numeric.")
    if (length(r_rmse) != n) stop("r_rmse must match length of r_est.")
  } else {
    r_rmse <- rep(NA_real_, n)
  }
  
  if (!is.null(p_rmse)) {
    if (!is.numeric(p_rmse)) stop("p_rmse must be numeric.")
    if (length(p_rmse) != n) stop("p_rmse must match length of p_est.")
  } else {
    p_rmse <- rep(NA_real_, n)
  }
  
  # ---- Compute k parameters ----
  k_est  <- r_est
  k_rmse <- r_rmse
  
  # ---- Compute b parameters safely ----
  b_est  <- (1 - p_est) / p_est
  b_rmse <- p_rmse / (p_est^2)
  
  # ---- Return tidy data frame ----
  list(
    k_est  = k_est,
    k_rmse = k_rmse,
    b_est  = b_est,
    b_rmse = b_rmse
  )
}
