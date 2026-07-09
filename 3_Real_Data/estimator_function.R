nb_params <- function(mat, seq_depth = NULL, max_depth = 0.3,
                      rmse = FALSE,
                      size_rmse_file = "error_fit_size.rds",
                      prob_rmse_file = "error_fit_prob.rds") {
  
  # ---- Validate inputs ----
  if (inherits(mat, "Matrix")) {
    if (!is.numeric(mat@x))
      stop("mat must contain numeric values")
  } else {
    if (!is.numeric(mat))
      stop("mat must contain numeric values")
  }
  if (any(mat < 0, na.rm = TRUE))
    stop("mat contains negative values")
  
  # --- Validate max_depth ---
  if (!is.numeric(max_depth) || length(max_depth) != 1)
    stop("max_depth must be a single numeric value")
  if (max_depth <= 0 || max_depth > 1)
    stop("max_depth must be >0 and <=1")
  
  # ---- Validate or generate seq_depth ----
  if (!is.null(seq_depth)) {
    if (!is.numeric(seq_depth))
      stop("seq_depth must be numeric")
    if (length(seq_depth) != ncol(mat))
      stop("seq_depth must have length equal to number of matrix columns")
    if (any(seq_depth <= 0))
      stop("seq_depth values must be >0")
    
    seq_depth <- max_depth * (seq_depth / max(seq_depth))
    
  } else {
    seq_sum <- colSums(mat)
    seq_depth <- seq_sum / mean(seq_sum)
    seq_depth <- max_depth * (seq_depth / max(seq_depth))
  }
  
  # Calculate proportions of zeroes per gene
  zero_prop_rows <- Matrix::rowSums(mat == 0) / ncol(mat)
  
  # ---- NB estimator ----
  nb_mle <- function(vals, depths) {
    
    if (var(vals) <= mean(vals)) {
      return(list(r_est = NA, p_est = 1, lambda = mean(vals)))
    }
    
    negloglik <- function(par, y, q) {
      r <- par[1]
      p <- par[2]
      p_prime <- p / (p + q * (1 - p))
      -sum(dnbinom(y, size = r, prob = p_prime, log = TRUE))
    }
    
    fit <- optim(
      par    = c(1, 0.5),
      fn     = negloglik,
      y      = vals,
      q      = depths,
      method = "L-BFGS-B",
      lower  = c(1e-6, 1e-6),
      upper  = c(1e3, 1 - 1e-6)
    )
    
    list(r_est = fit$par[1], p_est = fit$par[2])
  }
  
  # ---- Define progress printer ----
  
  progress_printer <- local({
    thresholds <- NULL
    
    function(i, total) {
      
      if (total < 10) {
        if (i == total) message("100% done")
        return(invisible())
      }
      
      if (is.null(thresholds)) {
        thresholds <<- ceiling(seq(0.1, 1, by = 0.1) * total)
      }
      
      if (i %in% thresholds) {
        pct <- which(thresholds == i) * 10
        message(pct, "% done (", i, "/", total, ")")
      }
    }
  })
  
  
  # ---- Estimate parameters ----
  G <- nrow(mat)
  
  r_vals <- numeric(G)
  p_vals <- numeric(G)
  
  for (i in seq_len(G)) {
    fit <- nb_mle(mat[i, ], seq_depth)
    r_vals[i] <- fit$r_est
    p_vals[i] <- fit$p_est
    
    progress_printer(i, G)
  }
  
  # ---- Compute Error ----
  if (rmse) {
    
    size_exists <- file.exists(size_rmse_file)
    prob_exists <- file.exists(prob_rmse_file)
    
    if (!size_exists || !prob_exists) {
      message("RMSE requested but error fit files not found: ",
              paste(c(size_rmse_file, prob_rmse_file)[!c(size_exists, prob_exists)],
                    collapse = ", "))
      
    } else {
      
      fit_rmse_size <- readRDS(size_rmse_file)
      fit_rmse_prob <- readRDS(prob_rmse_file)
      
      newdat <- data.frame(
        pred_sizes_corr = r_vals,
        pred_probs_corr = p_vals,
        zero_prop_corr  = zero_prop_rows,
        cells           = ncol(mat)
      )
      
      size_err <- abs(predict(fit_rmse_size, newdata = newdat))
      prob_err <- abs(predict(fit_rmse_prob, newdata = newdat))
    }
  }
  
  # ---- Return ----
  list(
    r_est      = r_vals,
    p_est      = p_vals,
    r_rmse     = unname(size_err),
    p_rmse     = unname(prob_err),
    zero_prop  = unname(zero_prop_rows),
    seq_depth  = unname(seq_depth)
  )
}

