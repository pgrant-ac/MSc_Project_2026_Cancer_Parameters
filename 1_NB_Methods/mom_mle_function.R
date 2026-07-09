nb_mle <- function(vals, depths = NULL){
  #' Function for estimating size(r) and probability(p) for count data 
  #' assuming a negative binomial distribution.
  #' Uses method of moments to calculate estimates.
  
  # If no depths supplied assume relative sequencing depths are equal
  if (is.null(depths)) {
    depths <- rep(1, length(vals))
  }
  
  # Ensure count values input is integer vector
  if (!(is.numeric(vals) &&
        is.atomic(vals) &&
        length(vals) >= 1 &&
        all(is.na(vals) | vals == as.integer(vals)))){
    stop("Input values must be a vector containing integer values")
  }
  # Check depths are correct numeric
  if (!(is.numeric(depths) && is.atomic(depths))) {
    stop("Depths must be a numeric vector")
  }
  
  # Ensure lengths are equal before cleaning
  if (length(vals) != length(depths)) {
    stop("Values and depths must be vectors of the same length")
  }
  
  # Remove NA in either vector
  na_idx <- is.na(vals) | is.na(depths)
  if (any(na_idx)) {
    warning(glue("{sum(na_idx)} NA depth(s)/value(s) removed from data"))
    vals  <- vals[!na_idx]
    depths <- depths[!na_idx]
  }
  
  # Remove negative values
  neg_idx <- vals < 0
  if (any(neg_idx)) {
    warning(glue("{sum(neg_idx)} negative value(s) removed"))
    vals  <- vals[!neg_idx]
    depths <- depths[!neg_idx]
  }
  
  # Remove depths <= 0
  zero_idx <- depths <= 0
  if (any(zero_idx)) {
    warning(glue("{sum(zero_idx)} depth(s) removed for being zero or less"))
    vals  <- vals[!zero_idx]
    depths <- depths[!zero_idx]
  }
  
  # Check still have data
  if (length(vals) == 0) {
    stop("No valid values remain after cleaning")
  }
  
  # Check if values' variance is equal or less than mean
  if (var(vals) <= mean(vals)) {
    warning("Distribution is underdispersed. Returning Poisson fit")
    return(list(r_est = NA, p_est = 1, lambda = mean(vals)))
  }
  
  # Method of moments likelihood
  gmm_nb <- function(par, y, q) {
    r <- exp(par[1])
    p <- plogis(par[2])
    
    # thinning-corrected probability
    p_prime <- p / (p + q * (1 - p))
    
    # theoretical moments
    mu  <- r * (1 - p_prime) / p_prime
    var <- mu + mu^2 / r
    
    # sample moments
    m1 <- mean(y)
    m2 <- var(y)
    
    # moment conditions
    g1 <- m1 - mean(mu)
    g2 <- m2 - mean(var)
    
    # GMM objective
    g1^2 + g2^2
  }
  
  fit <- optim(
    par    = c(0, 0),
    fn     = gmm_nb,
    y      = vals,
    q      = depths,
    method = "BFGS"
  )
  
  r_est <- exp(fit$par[1])
  p_est <- plogis(fit$par[2])
  
  list(r_est = r_est, p_est = p_est)
}