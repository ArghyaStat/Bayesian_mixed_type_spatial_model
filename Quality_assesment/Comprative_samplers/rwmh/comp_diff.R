## ===============================
## TRUE LOG-POSTERIOR (benchmark)
## ===============================

## Use TRUE values
beta.true  <- true.beta
Sigma.true <- true.Sigma
phi.true   <- true.phi
W.true     <- data$true.W.obs[obs.ord, , drop = FALSE]

## --- Likelihood ---
log.like.true <- log.likelihood(W.true, Y.obs.ord, family,
                                family.par = NULL,
                                link = NULL)$like.total

## --- Sigma prior ---
R.Sigma <- chol(Sigma.true)
Sigma.inv <- chol2inv(R.Sigma)
log.det.Sigma <- 2 * sum(log(diag(R.Sigma)))

log.prior.Sigma.true <- -(df.prior + q + 1)/2 * log.det.Sigma -
  0.5 * sum(diag(S.prior %*% Sigma.inv))

## --- Vecchia prior for W ---
log.prior.W.true <- vecchia_ll(
  W.ord = W.true,
  M.ord = X.obs.ord %*% beta.true,
  NNarray = NNarray.obs,
  coeff = coeff.obs,
  chol.V = R.Sigma,
  log = TRUE
)

## --- Beta prior (matrix normal) ---
E <- beta.true - M.prior

chol.V.prior <- chol(V.prior)
tmp <- forwardsolve(chol.V.prior, E)

quad.beta <- sum((tmp %*% Sigma.inv) * tmp)

log.prior.beta.true <- -0.5 * quad.beta

## --- Phi prior ---
log.prior.phi.true <- ifelse(phi.true > 0 & phi.true < b_phi, 0, -Inf)

## --- Final ---
log.post.W.true <- log.like.true + log.prior.W.true

log.post.full.true <- log.post.W.true +
  log.prior.beta.true +
  log.prior.Sigma.true +
  log.prior.phi.true

cat("True log-posterior (W only):", log.post.W.true, "\n")
cat("True log-posterior (full):", log.post.full.true, "\n")