# Elliptic slice sampler for W

update.W <- function(W.obs.ord, beta, chol.Sigma.inv, chol.K.obs.ord.inv,
                     Y.obs.ord, X.obs.ord, N.obs, log.like, family){

  q <- dim(chol.Sigma.inv)[1]

  gamma <- runif(1, min = 0, max = 2*pi)

  gamma_min <- gamma - 2*pi
  gamma_max <- gamma
  
  
  Z <- matrix(rnorm(N.obs*q, 0, 1), nrow = N.obs , ncol = q)
  
  Z_new <- Z %*% backsolve(chol.Sigma.inv, diag(q), transpose = TRUE)

  W.prior <- backsolve.spam(chol.K.obs.ord.inv, Z_new)
 
  W.prior.mean <- X.obs.ord %*% beta

  post.W <- log.likelihood(W.obs.ord, Y.obs.ord, family) 
 

  while (TRUE){

    can.W.obs.ord <- W.prior.mean + (W.obs.ord - W.prior.mean)*cos(gamma) + W.prior*sin(gamma)
    
    post.can.W <- log.likelihood(can.W.obs.ord, Y.obs.ord, family)
    
    log.like <- log.like + post.can.W
    

    log.r.W <- post.can.W - post.W

    if(log(runif(1)) < log.r.W){

      W.obs.ord <- can.W.obs.ord

      break

    }else{
      if(gamma < 0) gamma_min <- gamma else gamma_max <- gamma
      gamma <- runif(1, gamma_min, gamma_max)
    }
  }

  
  # log.like <- log.like + log.likelihood(W.obs.ord, Y.obs.ord, family)
  
  results <- list(W.obs.ord = W.obs.ord,
                  log.like = post.can.W)
  
  
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
  
  q <- dim(chol.Sigma)[1]

  beta.samp <- (matrix(rnorm(p*q, 0, 1), nrow = p , ncol = q) %*% chol.Sigma)

  beta.temp <- (chol.V.tilde %*% M.tilde.part + beta.samp)

  beta <- crossprod(chol.V.tilde, beta.temp)
  
  beta
  
}
