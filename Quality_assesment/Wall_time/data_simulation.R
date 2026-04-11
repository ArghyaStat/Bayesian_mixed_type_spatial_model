sim.data <- function(
    q,
    N,
    family,
    true.beta,
    true.Sigma,
    true.phi,
    true.nu,
    family.par = NULL,
    link = NULL,
    pred.prop,
    grid = TRUE,
    locations = NULL
) {
  
 
  
  if (length(family) != q)
    stop("Length of 'family' must equal q.")
  
  if (is.null(family.par))
    family.par <- vector("list", q)
  
  if (length(family.par) != q)
    stop("family.par must be a list of length q.")
  
  if (is.null(link))
    link <- rep("canonical", q)
  
  if (length(link) != q)
    stop("Length of 'link' must equal q.")
  
  
  if (!is.null(locations)) {
    
    if (!is.matrix(locations) || nrow(locations) != N)
      stop("locations must be an N x d matrix.")
    
  } else {
    
    if (grid) {
      
      m_x <- ceiling(sqrt(N))
      m_y <- ceiling(N / m_x)
      
      grid_x <- seq(0, 1, length.out = m_x)
      grid_y <- seq(0, 1, length.out = m_y)
      
      full_grid <- as.matrix(expand.grid(grid_x, grid_y))
      locations <- full_grid[1:N, , drop = FALSE]
      
    } else {
      
      locations <- cbind(runif(N), runif(N))
    }
  }
  
  d <- ncol(locations)
 
  # 2. Design matrix
  
  X <- cbind(1, locations)
  p <- ncol(X)
  
  mu.true <- X %*% true.beta   # N x q
  
  # 3. Spatial covariance (Matern)
  
  # distance
  distmat <- fields::rdist(locations)
  
  K.true <- fields::Matern(
    distmat,
    range = true.phi,
    smoothness = true.nu
  )
  
  chol.K.true <- chol(K.true)
  chol.Sigma.true <- chol(true.Sigma)
  
  W.error <- matrix(rnorm(N*q), N, q)
  
  # avoid large intermediate
  W.error <- W.error %*% chol.Sigma.true
  
  true.W <- mu.true + t(chol.K.true) %*% W.error
  
  # -------------------------------------------------
  # 4. Link inverse (family driven)
  # -------------------------------------------------
  
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
  
  # -------------------------------------------------
  # 5. Generate responses
  # -------------------------------------------------
  
  Y <- matrix(NA, N, q)
  
  for (j in 1:q) {
    
    fam  <- family[j]
    pars <- family.par[[j]]
    eta  <- true.W[, j]
    
    mu <- linkinv(eta, fam, link[j])
    
    if (fam == "Gaussian") {
      
      dispersion <- ifelse(is.null(pars$dispersion), 1, pars$dispersion)
      Y[, j] <- rnorm(N, mean = mu, sd = sqrt(dispersion))
      
    } else if (fam == "Poisson") {
      
      Y[, j] <- rpois(N, lambda = mu)
      
    } else if (fam == "Binomial") {
      
      size <- ifelse(is.null(pars$size), 1, pars$size)
      Y[, j] <- rbinom(N, size = size, prob = mu)
      
    } else if (fam == "Gamma") {
      
      shape <- ifelse(is.null(pars$shape), 1, pars$shape)
      rate  <- shape / mu
      Y[, j] <- rgamma(N, shape = shape, rate = rate)
      
    } else if (fam == "Negative-Binomial") {
      
      theta <- ifelse(is.null(pars$theta), 1, pars$theta)
      Y[, j] <- MASS::rnegbin(N, mu = mu, theta = theta)
      
    } else {
      
      stop(paste("Unsupported family:", fam))
    }
  }
  
  # -------------------------------------------------
  # 6. Train/Test split
  # -------------------------------------------------
  
  N.pred <- ceiling(pred.prop * N)
  N.obs  <- N - N.pred
  
  pred.indices <- sample(seq_len(N), N.pred)
  
  pred.locs <- locations[pred.indices, , drop = FALSE]
  obs.locs  <- locations[-pred.indices, , drop = FALSE]
  
  X.obs  <- X[-pred.indices, , drop = FALSE]
  X.pred <- X[pred.indices, , drop = FALSE]
  
  Y.obs       <- Y[-pred.indices, , drop = FALSE]
  Y.pred.true <- Y[pred.indices, , drop = FALSE]
  
  true.W.obs <- true.W[-pred.indices, , drop = FALSE]
  
  distmat.obs <- as.matrix(dist(obs.locs))
  diameter    <- max(distmat.obs)
  
  
  list(
    N = N,
    N.obs = N.obs,
    N.pred = N.pred,
    p = p,
    q = q,
    locations = locations,
    joint.locs = rbind(pred.locs, obs.locs),
    obs.locs = obs.locs,
    pred.locs = pred.locs,
    X.obs = X.obs,
    X.pred = X.pred,
    diameter = diameter,
    Y.obs = Y.obs,
    Y.pred.true = Y.pred.true,
    true.W.obs = true.W.obs,
    true.beta = true.beta,
    true.Sigma = true.Sigma,
    true.phi = true.phi,
    true.nu = true.nu,
    family = family,
    family.par = family.par,
    link = link
  )
}
