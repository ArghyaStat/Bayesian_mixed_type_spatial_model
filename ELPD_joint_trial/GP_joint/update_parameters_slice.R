# Column-wise Elliptical Slice Sampler for W
update.W <- function(W.obs.ord, beta, chol.Sigma.inv, chol.K.obs.ord.inv,
                     Y.obs.ord, X.obs.ord, N.obs, log.like, family){

  q <- ncol(W.obs.ord)

  ## Precompute XB
  W.prior.mean.full <- X.obs.ord %*% beta

  ## Current log-likelihood
  post.W <- log.likelihood(W.obs.ord, Y.obs.ord, family)

  ## Recover Sigma from its inverse Cholesky (cheap for small q)
  Sigma.inv <- crossprod(chol.Sigma.inv)
  Sigma <- solve(Sigma.inv)

  for (j in seq_len(q)) {

    ind.minus <- setdiff(seq_len(q), j)

    ## ---- conditional mean ----
    if (length(ind.minus) > 0) {

      Sigma.jm <- Sigma[j, ind.minus, drop = FALSE]
      Sigma.mm <- Sigma[ind.minus, ind.minus, drop = FALSE]
      Sigma.mm.inv <- solve(Sigma.mm)

      mu.j <- W.prior.mean.full[, j] +
        (W.obs.ord[, ind.minus, drop = FALSE] -
           W.prior.mean.full[, ind.minus, drop = FALSE]) %*%
        t(Sigma.mm.inv %*% t(Sigma.jm))

      Sigma.cond <- as.numeric(
        Sigma[j, j] -
          Sigma.jm %*% Sigma.mm.inv %*% t(Sigma.jm)
      )

    } else {
      mu.j <- W.prior.mean.full[, j]
      Sigma.cond <- Sigma[j, j]
    }

    ## ---- auxiliary direction ----
    z <- rnorm(N.obs)
    W.prior.j <- sqrt(Sigma.cond) *
      backsolve.spam(chol.K.obs.ord.inv, z)

    ## ---- ESS setup (same structure as your code) ----
    gamma <- runif(1, 0, 2*pi)
    gamma_min <- gamma - 2*pi
    gamma_max <- gamma

    logy <- post.W + log(runif(1))

    while (TRUE) {

      can.W.obs.ord <- W.obs.ord
      can.W.obs.ord[, j] <-
        mu.j +
        (W.obs.ord[, j] - mu.j) * cos(gamma) +
        W.prior.j * sin(gamma)

      post.can.W <- log.likelihood(can.W.obs.ord, Y.obs.ord, family)

      log.r.W <- post.can.W - post.W

      if (log(runif(1)) < log.r.W) {

        W.obs.ord <- can.W.obs.ord
        post.W <- post.can.W
        break

      } else {

        if (gamma < 0) gamma_min <- gamma else gamma_max <- gamma
        gamma <- runif(1, gamma_min, gamma_max)

      }
    }
  }

  list(W.obs.ord = W.obs.ord,
       log.like = post.W)
}



# Update of phi with truncated normal proposal 

update.phi <- function(phi, beta, chol.Sigma.inv, nu,
                       W.obs.ord, distobs.nn, NNarray.obs,
                       chol.K.obs.ord.inv, acc.phi, tuning.phi, b_phi){
  
  # Calculate the posterior for the current phi
  post.phi <- dmatnorm.sgv(X = W.obs.ord,
                           M = X.obs.ord %*% beta,
                           chol.U.prec = chol.K.obs.ord.inv,
                           chol.V.prec = chol.Sigma.inv,
                           log = TRUE) + dunif(phi, 0, b_phi, log = TRUE)
  
  # Generate a candidate phi using a truncated Gaussian distribution
  can.phi <- rtruncnorm(1, a = 0, b = b_phi, 
                        mean = phi, sd = sqrt(tuning.phi))
    
    # correlation matrix details for the candidate phi
    
    can.chol.K.obs.ord.inv <- U.sgv(distobs.nn, NNarray.obs, 
                                    can.phi, nu, m)
    
    # the posterior for the candidate phi
    post.can.phi <- dmatnorm.sgv(W.obs.ord,
                                 M = X.obs.ord %*% beta,
                                 chol.U.prec = can.chol.K.obs.ord.inv,
                                 chol.V.prec = chol.Sigma.inv,
                                 log = TRUE) + dunif(can.phi, 0, b_phi, log = TRUE)
    
    # Compute the truncated normal density at the current and candidate phi
    q_curr_given_can <- dtruncnorm(phi, a = 0, b = b_phi, mean = can.phi, 
                                   sd = sqrt(tuning.phi))
    q_can_given_curr <- dtruncnorm(can.phi, a = 0, b = b_phi, mean = phi, 
                                   sd = sqrt(tuning.phi))
    
    # the log acceptance ratio
    log.r.phi <- (post.can.phi - post.phi) + 
      log(q_curr_given_can) - log(q_can_given_curr)
    
    # Metropolis-Hastings acceptance step
    if(log(runif(1)) < log.r.phi){
      phi <- can.phi
      chol.K.obs.ord.inv <- can.chol.K.obs.ord.inv
      acc.phi <- acc.phi + 1
    }
  
  
  results <- list(phi = phi, 
                  chol.cormat.obs.inv = chol.K.obs.ord.inv, 
                  acc.phi = acc.phi)
  
  
  return(results)
}

#### Sigma-Beta MNIW block Gibbs update

# Sigma update

update.Sigma <- function(W.obs.ord, X.obs.ord, chol.K.obs.ord.inv, NNarray,
                         M.prior, chol.prec.V.prior, M.tilde.part, chol.V.tilde, 
                         S.prior, df.prior, N.obs){
  
  rss.W.part <- chol.K.obs.ord.inv %*% W.obs.ord
  
  rss.W.obs.ord <- crossprod(rss.W.part) 
  
  rss.beta.prior.part <- chol.prec.V.prior %*% M.prior
  
  rss.beta.prior <- crossprod(rss.beta.prior.part) 
  
  rss.beta.posterior.part <- chol.V.tilde %*% M.tilde.part
    
  rss.beta.posterior <- crossprod(rss.beta.posterior.part)
  
  S.tilde <- S.prior + rss.W.obs.ord + rss.beta.prior - rss.beta.posterior
  df.tilde <- df.prior + N.obs 
  
  Sigma <- riwish(v = df.tilde, S = S.tilde)
  
  chol.Sigma <- chol(Sigma)
  Sigma.inv <- chol2inv(chol.Sigma)
  chol.Sigma.inv <- chol(Sigma.inv)
  
  results <- list(Sigma = Sigma,
                  chol.Sigma.inv = chol.Sigma.inv,
                  chol.Sigma = chol.Sigma)
  
  return(results)
}


# beta update

update.beta <- function(M.tilde.part, chol.Sigma, chol.V.tilde){

  beta.samp <- (matrix(rnorm(p*q, 0, 1), nrow = p , ncol = q) %*% chol.Sigma)

  beta.temp <- (chol.V.tilde %*% M.tilde.part + beta.samp)

  beta <- crossprod(chol.V.tilde, beta.temp)
  
  beta
  
}
