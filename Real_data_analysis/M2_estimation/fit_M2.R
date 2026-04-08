### Real data analysis M1

rm(list = ls())

libraries <- c("spam", "this.path", "fields", "fBasics", "MCMCpack", "truncnorm", 
               "rlist", "foreach","doParallel", "FNN", "mcmcse", "sp", "gstat")


# Load libraries
invisible(lapply(libraries, library, character.only = TRUE))

mydir <- this.path::here()
setwd(mydir)

source("vecchia.R")
source("aux_functions_sep.R")
source("update_parameters_sep.R")
source("mixed_workflow_sep.R")
source("predictive_scores.R")

load(file.path("data_full.RData"), eva.data <- new.env())

full.data.DF <- eva.data$data_DF
full.data <- as.list(eva.data)

N.t <- 161
N.months <- 7
N.years <- 23
N.s <- 3503

t.index <- 140

spatial.subset <- seq(N.s * (t.index - 1) + 1, N.s * t.index)
data <- full.data.DF[spatial.subset, ]


k.folds <- 10
set.seed(35)

N <- nrow(data)
fold_id <- sample(rep(1:k.folds, length.out = N))

n.cores <- 10
cl <- makeCluster(n.cores)
registerDoParallel(cl)

results <- foreach(k = 1:k.folds,
                   .packages = libraries) %dopar% {
                     
                     ### =========================
                     ### 1. Split train / test
                     ### =========================
                     
                     train_idx <- which(fold_id != k)
                     test_idx  <- which(fold_id == k)
                     
                     ### Locations
                     locations <- cbind(data$lon, data$lat)
                     
                     obs.locs  <- locations[train_idx, , drop = FALSE]
                     pred.locs <- locations[test_idx, , drop = FALSE]
                     
                     ### Design matrix
                     X <- cbind(1, locations)
                     X.obs  <- X[train_idx, , drop = FALSE]
                     X.pred <- X[test_idx, , drop = FALSE]
                     
                     ### Responses
                     Y_BA  <- log(1 + data$BA)
                     Y_CNT <- data$CNT
                     Y <- cbind(Y_BA, Y_CNT)
                     
                     Y.obs       <- Y[train_idx, , drop = FALSE]
                     Y.pred.true <- Y[test_idx, , drop = FALSE]
                     
                     N.obs  <- nrow(obs.locs)
                     N.pred <- nrow(pred.locs)
                     
                     p <- ncol(X)
                     q <- ncol(Y)
                     
                     family <- c("Gaussian", "Poisson")
                     
                     ### =========================
                     ### 2. Vecchia ordering (TRAIN)
                     ### =========================
                     
                     m <- 20
                     
                     obs.ord <- max_min(obs.locs)
                     obs.locs.ord <- obs.locs[obs.ord, , drop = FALSE]
                     
                     distobs.ord <- rdist(obs.locs.ord)
                     diameter <- max(distobs.ord)
                     b_phi <- diameter/ 3
                     # b_phi <- median(distobs.ord) / 3
                     
                     Y.obs.ord <- Y.obs[obs.ord, , drop = FALSE]
                     X.obs.ord <- X.obs[obs.ord, , drop = FALSE]
                     
                     ### =========================
                     ### 3. Joint ordering (TRAIN + TEST)
                     ### =========================
                     
                     pred.ord <- max_min(pred.locs)
                     
                     ord <- c(obs.ord, pred.ord + N.obs)
                     locs.all <- rbind(obs.locs, pred.locs)
                     locs.ord <- locs.all[ord, , drop = FALSE]
                     
                     NNarray.all <- neighbor_matrix(locs.ord, m)
                     distall.nn  <- dist.nn(locs.ord, NNarray.all)
                     
                     
                     
                     
                     ### Reorder prediction
                     X.pred.ord <- X.pred[pred.ord, , drop = FALSE]
                     Y.pred.ord.true <- Y.pred.true[pred.ord, , drop = FALSE]
                     
                     ### =========================
                     ### 4. Priors & Initialization
                     ### =========================
                     
                     M.prior <- matrix(0, p, q)
                     V.prior <- 1e2 * diag(p)
                     a.prior <- 1e-2
                     b.prior <- 1e-2
                     
                     W.obs.ord <- cbind(Y.obs.ord[,1], log(1 + Y.obs.ord[,2]))
                     beta <- solve(crossprod(X.obs.ord)) %*% crossprod(X.obs.ord, W.obs.ord)
                     Sigma_full <- cov(W.obs.ord)
                     R <- chol(Sigma_full)
              
        
                     # 2. Whiten responses (remove cross-correlation)
                     W.white <- W.obs.ord %*% solve(R)
                     
                     # 3. Stack responses into a single field
                     values <- as.vector(W.white)
                     coords.rep <- obs.locs.ord[rep(1:nrow(obs.locs.ord), 
                                                    times = ncol(W.white)), ]
                     
                     # 4. Create spatial object
                     df.joint <- data.frame(
                       lon = coords.rep[,1],
                       lat = coords.rep[,2],
                       values = values
                     )
                     coordinates(df.joint) <- ~ lon + lat
                     
                     # 5. Empirical variogram
                     emp.vario.joint <- variogram(values ~ 1, df.joint)
                     
                     init.vgm <- vgm(
                       model = "Mat",
                       psill = var(values),
                       range = median(emp.vario.joint$dist),
                       nugget = 0
                     )
                     
                     fit.vgm <- fit.variogram(emp.vario.joint, model = init.vgm, fit.sills = TRUE)
                     
                     # 8. Extract phi (range parameter)
                     phi <- fit.vgm$range[2]
                     
                     # (optional) smoothness if needed
                     nu <- fit.vgm$kappa[2]
                     Sigma <- diag(diag(cov(W.obs.ord)))
                     
                     coeff.all <- vecchia_coeff(distall.nn, NNarray.all, phi, nu)
                     
                     tuning.phi <- 1e-3
                     niters <- 1e4
                     burnin <- 0
                     
                     
                     out <- mixed.workflow.sep(Y.obs.ord, X.obs.ord, obs.locs.ord, m, N.obs, p, q,
                                               nu, W.obs.ord, beta, Sigma, phi,
                                               M.prior, V.prior, a.prior, b.prior, b_phi,
                                               niters, tuning.phi,
                                               Y.pred.ord.true, X.pred.ord, N.pred,
                                               distall.nn, coeff.all, NNarray.all,
                                               family, family.par = NULL, link = NULL,
                                               burnin)
                     
                     out
                     
                     
                   }

stopCluster(cl)


# # Save results
save(results, file = "results_M2.RData")
