# Elliptical Slice Sampler for W

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
  
  W.obs.ord
  
}


update.phi <- function(phi, beta, chol.Sigma.inv, nu,
                       W.obs.ord, distobs.nn, NNarray.obs,
                       chol.K.obs.ord.inv, acc.phi, tuning.phi, b_phi){

 
  # Transform current phi -> theta (logit scale)

  theta <- log(phi / (b_phi - phi))

  # Current posterior (include Jacobian)
 
  log.jacobian.curr <- log(b_phi) + theta -
    2 * log(1 + exp(theta))

  post.phi <- dmatnorm.sgv(X = W.obs.ord,
                           M = X.obs.ord %*% beta,
                           chol.U.prec = chol.K.obs.ord.inv,
                           chol.V.prec = chol.Sigma.inv,
                           log = TRUE) + 
    dunif(phi, 0, b_phi, log = TRUE) + 
    log.jacobian.curr

  
  # Propose new theta (Gaussian RW)

  can.theta <- theta + rnorm(1, 0, sqrt(tuning.phi))

  # Transform back to phi
  can.phi <- b_phi * exp(can.theta) / (1 + exp(can.theta))

  can.chol.K.obs.ord.inv <- U.sgv(distobs.nn, NNarray.obs,
                                  can.phi, nu, m)

  # Candidate posterior (include Jacobian)

  log.jacobian.can <- log(b_phi) + can.theta -
    2 * log(1 + exp(can.theta))

  post.can.phi <- dmatnorm.sgv(W.obs.ord,
                               M = X.obs.ord %*% beta,
                               chol.U.prec = can.chol.K.obs.ord.inv,
                               chol.V.prec = chol.Sigma.inv,
                               log = TRUE) +
    dunif(can.phi, 0, b_phi, log = TRUE) +
    log.jacobian.can


  log.r.phi <- post.can.phi - post.phi


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


update.Sigma <- function(W.obs.ord, chol.K.obs.ord.inv,
                         M.prior, chol.prec.V.prior,
                         M.tilde.part, chol.V.tilde,
                         S.prior, df.prior, N.obs,
                         family,
                         fixed_variances) {
  
  q <- ncol(W.obs.ord)
  
  rss.W.part <- chol.K.obs.ord.inv %*% W.obs.ord
  rss.W.obs  <- crossprod(rss.W.part)
  
  rss.beta.prior <- crossprod(chol.prec.V.prior %*% M.prior)
  rss.beta.post  <- crossprod(chol.V.tilde %*% M.tilde.part)
  
  S.tilde <- S.prior + rss.W.obs + rss.beta.prior - rss.beta.post
  df.tilde <- df.prior + N.obs
  
  Sigma.star <- riwish(df.tilde, S.tilde)
  
  d.star <- sqrt(diag(Sigma.star))
  R.star <- diag(1/d.star) %*% Sigma.star %*% diag(1/d.star)
  
  D.new <- diag(q)
  
  for (j in 1:q) {
    if (family[j] == "Binomial") {
      D.new[j,j] <- sqrt(fixed_variances[j])
    } else {
      D.new[j,j] <- d.star[j]
    }
  }
  
  Sigma <- D.new %*% R.star %*% D.new
  
  chol.Sigma <- chol(Sigma)
  
  return(list(
    Sigma = Sigma,
    chol.Sigma = chol.Sigma,
    chol.Sigma.inv = chol(chol2inv(chol.Sigma))
  ))
}

# beta update

update.beta <- function(M.tilde.part, chol.Sigma, chol.V.tilde){
  
  q <- dim(chol.Sigma)[1]
  
  beta.samp <- (matrix(rnorm(p*q, 0, 1), nrow = p , ncol = q) %*% chol.Sigma)
  
  beta.temp <- (chol.V.tilde %*% M.tilde.part + beta.samp)
  
  beta <- crossprod(chol.V.tilde, beta.temp)
  
  beta
  
}