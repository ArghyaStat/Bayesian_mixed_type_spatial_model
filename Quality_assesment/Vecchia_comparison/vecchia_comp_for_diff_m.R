### Binomial Gaussian Poisson phi 0.1 joint dep

rm(list = ls())

libraries <- c("spam", "this.path", "fields", "fBasics", "MCMCpack", "truncnorm", 
               "rlist", "foreach","doParallel", "FNN", "mcmcse")


# Load libraries
invisible(lapply(libraries, library, character.only = TRUE))

mydir <- this.path::here()
setwd(mydir)


source("data_simulation.R")
source("vecchia.R")
source("aux_functions_joint.R")
source("update_parameters_joint.R")
source("mixed_workflow_joint.R")
source("predictive_scores.R")


  set.seed(35)
  
  p <- 3
  q <- 3
  
  
  true.beta <- matrix(
    c(1.0, -0.5,  0.8,
      3.0,  1.5, -2.0,
      -1.2,  0.0, 0.7),
    nrow = p, ncol = q, byrow = TRUE
  )
  
  
  true.Sigma <- matrix(
    c(9, 5, 4,
      5, 3, 9/4,
      4, 9/4, 2),
    nrow = q, ncol = q, byrow = TRUE)
  
  
  true.phi <- 0.3
  true.nu <- 0.5
  pred.prop <- 0.2
  
  family <- c("Binomial", "Gaussian", "Poisson")
  
  data <- sim.data(q = q, N = 5e2, 
                   family = family,
                   true.beta = true.beta,
                   true.Sigma = true.Sigma, 
                   true.phi = true.phi,
                   true.nu = true.nu,
                   pred.prop = 0.2)
  
  
  ### Vecchia specifications and pre-computations ###
  
  obs.locs <- data$obs.locs
  N.obs <- data$N.obs
  Y.obs <- data$Y.obs
  X.obs <- data$X.obs
  
  
  obs.ord <- max_min(obs.locs) 
  # Assume max_min returns max-min order indices
  obs.locs.ord <- obs.locs[obs.ord, , drop = FALSE]  # Reorder locations based on max-min ordering
  distobs.ord <- rdist(obs.locs.ord)
  diameter <-  max(distobs.ord)
  b_phi <- diameter/3
  
  
  
  
  #### Prediction set up ####
  
  pred.locs <- data$pred.locs
  N.pred <- data$N.pred
  pred.ord <- max_min(pred.locs)
  ord <- c(obs.ord, pred.ord + nrow(obs.locs))
  locs.ord <- rbind(obs.locs, pred.locs)[ord, , drop = FALSE]
  distlocs.ord <- rdist(locs.ord)
  
  
  
  phi <- data$true.phi
  nu <- data$true.nu
  
  
  m.values <- c(10, 20, 50, 500)
  
  cov.err.weak  <- numeric(length(m.values))
  cov.err.strong <- numeric(length(m.values))
  
  phi.values <- c(0.1, 0.3)
  
  for(phi_idx in 1:2){
    
    phi <- phi.values[phi_idx]
    
    # TRUE covariance (ALL locations)
    C.true <- Matern(distlocs.ord, range = phi, smoothness = nu)
    
    for(k in seq_along(m.values)){
      
      m <- m.values[k]
      
      cat("phi =", phi, " m =", m, "\n")  # progress
      
      NNarray.all <- neighbor_matrix(locs.ord, m)
      distall.nn  <- dist.nn(locs.ord, NNarray.all)
      
      coeff.all <- vecchia_coeff(distall.nn, NNarray.all, phi, nu)
      
      C.vecchia <- vecchia_cov_from_coeff(coeff.all, NNarray.all)
      
      err <- norm(C.true - C.vecchia, "F") / norm(C.true, "F")
      
      if(!is.finite(err)) err <- NA
      
      if(phi_idx == 1){
        cov.err.weak[k] <- err
      } else {
        cov.err.strong[k] <- err
      }
    }
  }
  
  ### Figure 6 (left panel)
  
  
  pdf("vecchia_elbow_plot.pdf", width = 9, height = 7.5)
  
  par(mar = c(3.5, 4.5, 1, 1), mgp = c(2.5,0.8,0), tcl = -0.3)
  
  valid.idx <- which(is.finite(cov.err.weak) & is.finite(cov.err.strong))
  
  m.plot <- m.values[valid.idx]
  weak.plot <- cov.err.weak[valid.idx]
  strong.plot <- cov.err.strong[valid.idx]
  
  plot(m.plot, weak.plot,
       type = "b",
       pch = 19,
       lwd = 2,
       col = "blue",
       log = "x",                 
       xaxt = "n",
       ylim = c(0, max(c(weak.plot, strong.plot))),
       xlab = "Neighbor size (m, log scale)",
       ylab = "Relative Frobenius norm",
       cex.lab = 2, cex.axis = 1.4)
  
  # Custom ticks at actual m values
  axis(1, at = m.plot, labels = m.plot, cex.axis = 1.6)
  
  lines(m.plot, strong.plot,
        type = "b",
        pch = 17,
        lwd = 2,
        lty = 2,
        col = "red")
  
  points(m.plot, weak.plot, pch = 19, col = "blue")
  points(m.plot, strong.plot, pch = 17, col = "red")
  
  legend("topright",
         legend = c(expression("Weak ("*phi==0.1*")"),
                    expression("Strong ("*phi==0.3*")")),
         col = c("blue", "red"),
         lty = c(1, 2),
         pch = c(19, 17),
         cex = 1.6,
         lwd = 2,
         bty = "n")
  
  box()
  
  dev.off()
  
  #### Prior specifications ####
  
  M.prior <- matrix(0, p, q)
  V.prior <-  1e2*diag(p)
  S.prior <-  diag(q)
  df.prior <- q + 1
  
  
  X.pred <- data$X.pred
  Y.pred.true <- data$Y.pred.true
  
  X.pred.ord <- X.pred[pred.ord, , drop = FALSE]
  Y.pred.ord.true <- Y.pred.true[pred.ord, , drop = FALSE]
  
  
  
  # Reorder Y.obs, X.obs, and W.obs accordingly
  Y.obs.ord <- Y.obs[obs.ord, , drop = FALSE]
  X.obs.ord <- X.obs[obs.ord, , drop = FALSE]
  
  
  # Sample initial parameters from MLE
  
  beta <- data$true.beta
  Sigma <- data$true.Sigma
  W.obs.ord <- data$true.W.obs[obs.ord, , drop = FALSE]
  
  
  tuning.phi <- 1e-3
  
  # Number of iterations
  niters <- 1e4
  
 
  # Average wall time reporting for different values of m
  
  m.values <- c(10, 20, 50, 500)
  
  est.time   <- numeric(length(m.values))
  pred.time  <- numeric(length(m.values))
  total.time <- numeric(length(m.values))
  
  for(k in seq_along(m.values)){
    
    m <- m.values[k]
    
    cat("Running for m =", m, "\n")
    
    # Recompute Vecchia structures (IMPORTANT)
    NNarray.all <- neighbor_matrix(locs.ord, m)
    distall.nn  <- dist.nn(locs.ord, NNarray.all)
    coeff.all   <- vecchia_coeff(distall.nn, NNarray.all, phi, nu)
    
    out <- mixed.workflow.joint(
      Y.obs.ord, X.obs.ord, obs.locs.ord, m, N.obs, p, q,
      nu, W.obs.ord, beta, Sigma, phi,
      M.prior, V.prior, S.prior, df.prior, b_phi,
      niters, tuning.phi,
      Y.pred.ord.true, X.pred.ord, N.pred,
      distall.nn, coeff.all, NNarray.all,
      family, family.par = NULL, link = NULL
    )
    
    est.time[k]   <- out$est.time
    pred.time[k]  <- out$pred.time
    total.time[k] <- out$total.time
  }
  
  walltime <- cbind(est.time, pred.time, total.time)
  
  save(walltime, file = "bgp_walltime500.RData")
  
  ### Figure 6 (right panel)
  
  pdf("vecchia_runtime_plot_secs.pdf", width = 9, height = 7.5)
  
  par(mar = c(3.5, 4.5, 1, 1), mgp = c(2.5,0.8,0), tcl = -0.3)
  
  # Convert ALL to seconds
  est.sec   <- 60 * est.time
  pred.sec  <- 60 * pred.time
  total.sec <- 60 * total.time
  
  # Ensure positivity for log scale
  eps <- 1e-8
  est.sec   <- pmax(est.sec, eps)
  pred.sec  <- pmax(pred.sec, eps)
  total.sec <- pmax(total.sec, eps)
  
  plot(m.values, est.sec,
       type = "b",
       pch = 19,
       lwd = 2,
       col = "blue",
       log = "xy",
       xaxt = "n",
       ylim = range(c(est.sec, pred.sec, total.sec)),
       xlab = expression("Neighbor size ("*m*", log scale)"),
       ylab = "Runtime (seconds, log scale)",
       cex.lab = 2, cex.axis = 1.6)
  
  axis(1, at = m.values, labels = m.values, cex.axis = 1.6)
  
  
  lines(m.values, pred.sec,
        type = "b",
        pch = 17,
        lwd = 2,
        lty = 2,
        col = "red")
  
  lines(m.values, total.sec,
        type = "b",
        pch = 15,
        lwd = 2,
        lty = 3,
        col = "darkgreen")
  
  # Overlay points (NOW consistent)
  points(m.values, est.sec,   pch = 19, col = "blue")
  points(m.values, pred.sec,  pch = 17, col = "red")
  points(m.values, total.sec, pch = 15, col = "darkgreen")
  
  legend("topleft",
         legend = c("Estimation time",
                    "Prediction time",
                    "Total time"),
         col = c("blue", "red", "darkgreen"),
         lty = c(1, 2, 3),
         pch = c(19, 17, 15),
         cex = 1.6,
         lwd = 2,
         bty = "n")
  
  box()
  
  dev.off()
