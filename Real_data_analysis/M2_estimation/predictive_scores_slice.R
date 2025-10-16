### RMSPE

RMSPE <- function(Y.pred.samples, Y.true, pred.iters){
  
  out <- 0
  
  q <- ncol(Y.true)
  
  for(i in 1:pred.iters){
    
    for(j in 1:q){
      
      out <- out + sum((Y.pred.samples[[i]][,q] - Y.true[,q])^2)
    }
    
  }
  
  return(sqrt(out/pred.iters))
}



# Function to compute Predictive Coverage

pred_coverage <- function(Y_true, Y_pred_samples, alpha = 0.05) {
  N_pred <- nrow(Y_true)
  q <- ncol(Y_true)
  pred_iters <- length(Y_pred_samples)
  
  coverage <- matrix(NA, nrow = N_pred, ncol = q)
  
  for (j in 1:q) {
    pred_samples_j <- do.call(cbind, lapply(Y_pred_samples, function(x) x[, j]))
    lower_bound <- apply(pred_samples_j, 1, quantile, probs = alpha / 2)
    upper_bound <- apply(pred_samples_j, 1, quantile, probs = 1 - alpha / 2)
    
    coverage[, j] <- (Y_true[, j] >= lower_bound) & (Y_true[, j] <= upper_bound)
  }
  
  return(mean(coverage))
}


# Function to compute Energy Score
energy_score <- function(Y_true, Y_pred_samples) {
  N_pred <- nrow(Y_true)
  q <- ncol(Y_true)
  M <- length(Y_pred_samples)
  
  energy_scores <- numeric(N_pred)
  
  for (i in 1:N_pred) {
    pred_samples_i <- sapply(Y_pred_samples, function(x) x[i, ], simplify = "array")
    mean_term <- mean(apply(pred_samples_i, 2, function(x) sqrt(sum((Y_true[i, ] - x)^2))))
    pairwise_term <- mean(sapply(1:(M-1), function(m) sqrt(sum((pred_samples_i[, m] - pred_samples_i[, m+1])^2))))
    energy_scores[i] <- mean_term - 0.5 * pairwise_term
  }
  
  return(mean(energy_scores))
}

#### Poisson scores ####

crps_pois <- function(y, lambda) {
  c1 <- (y - lambda) * (2 * ppois(y, lambda) - 1)
  c2 <- 2 * dpois(floor(y), lambda) -
    besselI(2 * lambda, 0, expon.scaled = TRUE) -
    besselI(2 * lambda, 1, expon.scaled = TRUE)
  return(c1 + lambda * c2)
}


logs_pois <- function(y, lambda) {
  -dpois(y, lambda, log = TRUE)
}


dss_pois <- function(y, lambda) {
  lambda[lambda <= 0] <- NaN
  (y - lambda)^2 / lambda + log(lambda)
}

#### Normal scores ####

crps_norm <- function(y, mean = 0, sd = 1, location = mean, scale = sd) {
  if (!missing(mean) && !missing(location))
    stop("specify 'mean' or 'location' but not both")
  if (!missing(sd) && !missing(scale))
    stop("specify 'sd' or 'scale' but not both")
  if (!identical(location, 0)) y <- y - location
  if (identical(scale, 1)) {
    y * (2 * pnorm(y) - 1) + (sqrt(2) * exp(-0.5 * y^2) - 1) / sqrt(pi)
  } else {
    z <- y / scale
    z[y == 0 & scale == 0] <- 0
    y * (2 * pnorm(y, sd = scale) - 1) +
      scale * (sqrt(2) * exp(-0.5 * z^2) - 1) / sqrt(pi)
  }
}


logs_norm <- function(y, mean = 0, sd = 1, location = mean, scale = sd) {
  if (!missing(mean) && !missing(location))
    stop("specify 'mean' or 'location' but not both")
  if (!missing(sd) && !missing(scale))
    stop("specify 'sd' or 'scale' but not both")
  -dnorm(y, location, scale, log = TRUE)
}

dss_norm <- function(y, mean = 0, sd = 1, location = mean, scale = sd) {
  if (!missing(mean) && !missing(location))
    stop("specify 'mean' or 'location' but not both")
  if (!missing(sd) && !missing(scale))
    stop("specify 'sd' or 'scale' but not both")
  if (!identical(location, 0)) y <- y - location
  if (identical(scale, 1)) {
    y^2
  } else {
    scale[scale <= 0] <- NaN
    (y / scale)^2 + 2*log(scale)
  }
}


###### Gamma scores

crps_gamma <- function(y, shape, rate = 1, scale = 1/rate) {
  if (!missing(rate) && !missing(scale))
    stop("specify 'rate' or 'scale' but not both")
  p1 <- pgamma(y, shape, scale = scale)
  p2 <- pgamma(y, shape + 1, scale = scale)
  y * (2 * p1 - 1) - scale * (shape * (2 * p2 - 1) + 1 / beta(.5, shape))
}


logs_gamma <- function(y, shape, rate = 1, scale = 1/rate) {
  if (!missing(rate) && !missing(scale))
    stop("specify 'rate' or 'scale' but not both")
  -dgamma(y, shape, scale = scale, log = TRUE)
}


dss_gamma <- function(y, shape, rate = 1, scale = 1/rate) {
  if (!missing(rate) && !missing(scale))
    stop("specify 'rate' or 'scale' but not both")
  ms <- sqrt(shape)
  scale[scale <= 0] <- NaN
  s <- ms * scale
  (y / s - ms)^2 + 2*log(s)
}


### Binomial scores ### 

crps_binom <- function(y, size, prob) {
  n_param <- max(length(size), length(prob))
  n_y <- length(y)
  size <- rep(size, length.out = n_param)
  prob <- rep(prob, length.out = n_param)
  
  if (n_y <= n_param) {
    y <- rep(y, length.out = n_param)
    sapply(
      seq_along(y),
      function(i) {
        y <- y[i]
        size <- size[i]
        prob <- prob[i]
        if (anyNA(c(y, size, prob))) return(y * size * prob)
        size_rounded <- round(size)
        tol <- .Machine$double.eps^0.5
        size <- if (abs(size - size_rounded) < tol) {
          size_rounded
        } else {
          warning(sprintf("non-integer n = %.6f", size))
          return(NaN)
        }
        if (size >= 0) {
          x <- seq.int(0, size, 1)
          w <- dbinom(x, size, prob)
          a <- pbinom(x, size, prob) - 0.5 * w
          2 * sum(w * ((y < x) - a) * (x - y))
        } else {
          NaN
        }
      }
    )
  } else {
    list_param <- lapply(
      seq_along(size),
      function(i) {
        size <- size[i]
        prob <- prob[i]
        if (anyNA(c(size, prob))) {
          typeNA <- size * prob
          return(list(x = typeNA, w = typeNA, a = typeNA))
        }
        size_rounded <- round(size)
        tol <- .Machine$double.eps^0.5
        listNaN <- list(x = NaN, w = NaN, a = NaN)
        size <- if (abs(size - size_rounded) < tol) {
          size_rounded
        } else {
          warning(sprintf("non-integer n = %.6f", size))
          return(listNaN)
        }
        if (size >= 0) {
          x <- seq.int(0, size, 1)
          w <- dbinom(x, size, prob)
          a <- pbinom(x, size, prob) - 0.5 * w
          list(x = x, w = w, a = a)
        } else {
          listNaN
        }
      }
    )
    list_param <- rep(list_param, length.out = n_y)
    sapply(
      seq_along(y),
      function(i) {
        with(list_param[[i]], 2 * sum(w * ((y[i] < x) - a) * (x - y[i])))
      }
    )
  }
}


logs_binom <- function(y, size, prob) {
  -dbinom(y, size, prob, log = TRUE)
}



dss_binom <- function(y, size, prob) {
  # Calculate the mean (mu_P) and standard deviation (sigma_P)
  mu_P <- size * prob
  sigma_P <- sqrt(size * prob * (1 - prob))
  
  # Calculate the Dawid and Sebastiani score
  ds_score <- ((y - mu_P) / sigma_P)^2 + 2 * log(sigma_P)
  
  return(ds_score)
}


#### Negative Binomial scores ######

crps_nbinom <- function(y, size, prob, mu) {
  # check from stats::pnbinom
  if (!missing(mu)) {
    if (!missing(prob))
      stop("specify 'prob' or 'mu' but not both")
    prob <- size / (size + mu)
  }
  if (!requireNamespace("hypergeo", quietly = TRUE)) {
    stop(paste(
      "Calculations require an implementation of the gaussian hypergeometric function.",
      "Please install the following package: hypergeo (>= 1.0)",
      sep = "\n"))
  }
  c1 <- y * (2 * pnbinom(y, size, prob) - 1)
  c2 <- (1 - prob) / prob ^ 2
  c3 <- (prob * (2 * pnbinom(y - 1, size + 1, prob) - 1)
         #+ Re(hypergeo::hypergeo(size + 1, 0.5, 2, -4 * c2)))
         + hypergeo_0.5_2(size + 1, -4 * c2))
  return(c1 - size * c2 * c3)
}

logs_nbinom <- function(y, size, prob, mu) {
  if (!missing(prob)) {
    if (!missing(mu))
      stop("specify 'prob' or 'mu' but not both")
    -dnbinom(y, size, prob, log = TRUE)
  } else {
    -dnbinom(y, size, mu = mu, log = TRUE)
  }
}


dss_nbinom <- function(y, size, prob, mu) {
  size[size <= 0] <- NaN
  if (!missing(prob)) {
    if (!missing(mu))
      stop("specify 'prob' or 'mu' but not both")
    prob[prob <= 0] <- NaN
    prob[prob >= 1] <- NaN
    mu <- size * (1 - prob) / prob
    v <- mu / prob
  } else {
    mu[mu < 0] <- NaN
    v <- mu * (1 + mu / size)
  }
  (y - mu)^2 / v + log(v)
}


compute_crps <- function(Y_true, Y_pred_samples, family) {
  N.pred <- nrow(Y_true)
  q <- ncol(Y_true)
  
  # Function to compute CRPS for a single observation and its predictions
  crps_func <- function(y_true, y_pred_samples, fam) {
    if (fam == "Gaussian") {
      return(mean(sapply(y_pred_samples, function(mu) crps_norm(y_true, mean = mu, sd = 1))))
    } else if (fam == "Poisson") {
      return(mean(sapply(y_pred_samples, function(lambda) crps_pois(y_true, lambda = lambda))))
    } else if (fam == "Binomial") {
      return(mean(sapply(y_pred_samples, function(prob) crps_binom(y_true, size = 1, prob = prob))))
    } else if (fam == "Negative-Binomial") {
      return(mean(sapply(y_pred_samples, function(mu) crps_nbinom(y_true, size = 1, mu = mu))))
    } else if (fam == "Gamma") {
      return(mean(sapply(y_pred_samples, function(shape) crps_gamma(y_true, shape = shape, rate = 1))))
    } else {
      return(NA)  # Unsupported family
    }
  }
  
  # Compute CRPS for each column (response variable)
  crps_scores <- sapply(1:q, function(j) {
    y_true_col <- Y_true[, j]  # Extract column j
    y_pred_col <- lapply(Y_pred_samples, function(sample) sample[, j])  # Extract corresponding predictions
    sapply(1:N.pred, function(i) crps_func(y_true_col[i], sapply(y_pred_col, "[", i), family[j]))
  })
  
  return(mean(crps_scores, na.rm = TRUE))
}

compute_dss <- function(Y_true, Y_pred_samples, family) {
  N.pred <- nrow(Y_true)
  q <- ncol(Y_true)
  
  dss_scores <- sapply(1:q, function(j) {
    y_true_col <- Y_true[, j]
    y_pred_col <- lapply(Y_pred_samples, function(sample) sample[, j])
    
    sapply(1:N.pred, function(i) {
      y <- y_true_col[i]
      preds <- sapply(y_pred_col, "[", i)
      
      if (family[j] == "Gaussian") {
        return(mean(sapply(preds, function(mu) dss_norm(y, mean = mu, sd = 1))))
      } else if (family[j] == "Poisson") {
        return(mean(sapply(preds, function(lambda) dss_pois(y, lambda = lambda))))
      } else if (family[j] == "Binomial") {
        return(mean(sapply(preds, function(prob) dss_binom(y, size = 1, prob = prob))))
      } else if (family[j] == "Negative-Binomial") {
        return(mean(sapply(preds, function(mu) dss_nbinom(y, size = 1, mu = mu))))
      } else if (family[j] == "Gamma") {
        return(mean(sapply(preds, function(shape) dss_gamma(y, shape = shape, rate = 1))))
      } else {
        return(NA)
      }
    })
  })
  
  return(mean(dss_scores, na.rm = TRUE))
}




compute_logs <- function(Y_true, Y_pred_samples, family) {
  N.pred <- nrow(Y_true)
  q <- ncol(Y_true)
  
  logs_scores <- sapply(1:q, function(j) {
    y_true_col <- Y_true[, j]
    y_pred_col <- lapply(Y_pred_samples, function(sample) sample[, j])
    
    sapply(1:N.pred, function(i) {
      y <- y_true_col[i]
      preds <- sapply(y_pred_col, "[", i)
      
      # Compute log score per family
      score <- if (family[j] == "Gaussian") {
        mean(sapply(preds, function(mu) logs_norm(y, mean = mu, sd = 1)))
      } else if (family[j] == "Poisson") {
        mean(sapply(preds, function(lambda) logs_pois(y, lambda = pmax(lambda, .Machine$double.eps))))
      } else if (family[j] == "Binomial") {
        mean(sapply(preds, function(prob) logs_binom(y, size = 1, prob = pmin(pmax(prob, .Machine$double.eps), 1 - .Machine$double.eps))))
      } else if (family[j] == "Negative-Binomial") {
        mean(sapply(preds, function(mu) logs_nbinom(y, size = 1, mu = pmax(mu, .Machine$double.eps))))
      } else if (family[j] == "Gamma") {
        mean(sapply(preds, function(shape) logs_gamma(y, shape = pmax(shape, .Machine$double.eps), rate = 1)))
      } else {
        NA
      }
      
      # Ensure finite values only
      if (!is.finite(score)) return(NA)
      return(score)
    })
  })
  
  return(mean(logs_scores, na.rm = TRUE))
}

