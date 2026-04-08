rm(list = ls())

libraries <- c("spam", "this.path", "fields", "fBasics", "MCMCpack", "truncnorm", 
               "rlist", "FNN", "mcmcse")

invisible(lapply(libraries, library, character.only = TRUE))

mydir <- this.path::here()
setwd(mydir)

source("data_simulation.R")
source("vecchia.R")
source("aux_functions_joint.R")
source("update_parameters_ess_joint.R")
source("mixed_workflow_joint.R")

set.seed(35)

### Model
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
  nrow = q, byrow = TRUE)

true.phi <- 0.1
true.nu <- 0.5

family <- c("Binomial", "Gaussian", "Poisson")

### Data (NO prediction)
data <- sim.data(q = q, N = 5e2, 
                 family = family,
                 true.beta = true.beta,
                 true.Sigma = true.Sigma, 
                 true.phi = true.phi,
                 true.nu = true.nu,
                 pred.prop = 0)

### Priors
M.prior <- matrix(0, p, q)
V.prior <- 1e2 * diag(p)
S.prior <- diag(q)
df.prior <- q + 1

obs.locs <- data$obs.locs
N.obs <- nrow(obs.locs)
Y.obs <- data$Y.obs
X.obs <- data$X.obs

m <- 20
  
# min(20, N.obs - 1)

obs.ord <- seq_len(N.obs)

obs.locs.ord <- obs.locs[obs.ord, , drop = FALSE]

## neighbor structure
NNarray.obs <- neighbor_matrix(obs.locs.ord, m)

## FULL distance matrix (NOT list)
distobs.ord <- rdist(obs.locs.ord)

distobs.nn <- dist.nn(obs.locs.ord, NNarray.obs)


## Vecchia coeff
coeff.obs <- vecchia_coeff(distobs.nn, NNarray.obs, true.phi, true.nu)

## reorder data
Y.obs.ord <- Y.obs[obs.ord, , drop = FALSE]
X.obs.ord <- X.obs[obs.ord, , drop = FALSE]

## phi bound
b_phi <- max(distobs.ord) / 3



### Initial
beta <- true.beta
Sigma <- true.Sigma
phi <- true.phi
nu <- true.nu
W.obs.ord <- data$true.W.obs[obs.ord, , drop = FALSE]



### MCMC
niters <- 1e4
tuning.phi <- 5e-4



out.joint.ess <- mixed.workflow.joint(
  Y.obs.ord, X.obs.ord, obs.locs.ord, m, N.obs, p, q,
  nu, W.obs.ord, beta, Sigma, phi,
  M.prior, V.prior, S.prior, df.prior, b_phi,
  niters, tuning.phi,
  distobs.nn, coeff.obs, NNarray.obs,
  family, family.par = NULL, link = NULL
)

save(out.joint.ess, file = "bgp_jointESS_N500.RData")