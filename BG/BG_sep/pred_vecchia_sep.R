W.pred.ord <- function(W.obs.ord, beta, chol.Sigma, phi, nu, m,
                       X.obs.ord, X.pred.ord, U.obs.col, U.pred.col,
                       chol.W.pred.var, N.obs, N.pred, q){
  
  W.pred.mean.part <- U.obs.col %*% (W.obs.ord - X.obs.ord %*% beta)
  W.pred.mean.tilde <- crossprod.spam(U.pred.col,  W.pred.mean.part)
  
  W.pred.samp <- (matrix(rnorm(N.pred*q, 0, 1), nrow = N.pred , ncol = q) %*% chol.Sigma)
  
  W.pred.temp <- chol.W.pred.var %*% W.pred.mean.tilde +  W.pred.samp
  
  W.pred.ord <- X.pred.ord %*% beta -  crossprod(chol.W.pred.var,  W.pred.temp)
  
  W.pred.ord
  
}

Y.pred.ord <- function(W.pred.ord, N.pred, q, family, family.par = NULL, link = NULL) {
  
  # Checks
  
  if (!all(dim(W.pred.ord) == c(N.pred, q)))
    stop("W.pred.ord must be N.pred x q.")
  
  if (length(family) != q)
    stop("Length of 'family' must equal q.")
  
  if (is.null(family.par))
    family.par <- vector("list", q)
  
  if (length(family.par) != q)
    stop("family.par must be list of length q.")
  
  if (is.null(link))
    link <- rep("canonical", q)
  
  if (length(link) != q)
    stop("Length of 'link' must equal q.")
  
  # Link inverse
  
  linkinv <- function(eta, fam, linkname) {
    
    if (linkname == "canonical") {
      
      switch(fam,
             "Gaussian" = eta,
             "Poisson"  = exp(eta),
             "Binomial" = plogis(eta),
             "Gamma"    = exp(eta),
             "Negative-Binomial" = exp(eta),
             stop(paste("Unsupported family:", fam)))
      
    } else if (linkname == "logit") {
      
      plogis(eta)
      
    } else if (linkname == "probit") {
      
      pnorm(eta)
      
    } else if (linkname == "cloglog") {
      
      1 - exp(-exp(eta))
      
    } else if (linkname == "identity") {
      
      eta
      
    } else if (linkname == "log") {
      
      exp(eta)
      
    } else {
      
      stop("Unsupported link.")
    }
  }
  
  # Generate predictive responses
  
  Y.mat <- matrix(NA, N.pred, q)
  
  for (j in seq_len(q)) {
    
    fam  <- family[j]
    pars <- family.par[[j]]
    eta  <- W.pred.ord[, j]
    
    mu <- linkinv(eta, fam, link[j])
    
    if (fam == "Gaussian") {
      
      dispersion <- ifelse(is.null(pars$dispersion), 1, pars$dispersion)
      
      Y.mat[, j] <- rnorm(
        N.pred,
        mean = mu,
        sd   = sqrt(dispersion)
      )
      
    } else if (fam == "Poisson") {
      
      Y.mat[, j] <- rpois(
        N.pred,
        lambda = mu
      )
      
    } else if (fam == "Binomial") {
      
      size <- ifelse(is.null(pars$size), 1, pars$size)
      
      Y.mat[, j] <- rbinom(
        N.pred,
        size = size,
        prob = mu
      )
      
    } else if (fam == "Gamma") {
      
      shape <- ifelse(is.null(pars$shape), 1, pars$shape)
      rate  <- shape / mu
      
      Y.mat[, j] <- rgamma(
        N.pred,
        shape = shape,
        rate  = rate
      )
      
    } else if (fam == "Negative-Binomial") {
      
      theta <- ifelse(is.null(pars$theta), 1, pars$theta)
      
      Y.mat[, j] <- rnbinom(
        N.pred,
        mu   = mu,
        size = theta
      )
      
    } else {
      
      stop(paste("Unsupported family:", fam))
    }
  }
  
  Y.mat
}
