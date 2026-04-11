rm(list = ls())

#### ---------------- LIBRARIES ---------------- ####
libraries <- c("spam", "this.path", "fields", "fBasics", "MCMCpack", 
               "truncnorm", "rlist", "foreach","doParallel", "FNN", "mcmcse",
               "xtable")

invisible(lapply(libraries, library, character.only = TRUE))

mydir <- this.path::here()
setwd(mydir)

#### ---------------- SOURCES ---------------- ####
source("data_simulation.R")
source("vecchia.R")
source("aux_functions_joint.R")
source("update_parameters_joint.R")
source("mixed_workflow_joint.R")
source("predictive_scores.R")

#### ---------------- SETUP ---------------- ####
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
  nrow = q, ncol = q, byrow = TRUE
)

true.phi <- 0.3
true.nu  <- 0.5
pred.prop <- 0.2

family <- c("Binomial", "Gaussian", "Poisson")

#### ---------------- PRIOR ---------------- ####
M.prior <- matrix(0, p, q)
V.prior <- 1e2 * diag(p)
S.prior <- diag(q)
df.prior <- q + 1

#### ---------------- SETTINGS ---------------- ####
N <- c(100, 500, 2500)

niters <- 1e2
tuning.phi <- 1e-3

#### ---------------- STORAGE ---------------- ####
runtime.df <- data.frame()

#### ---------------- MAIN LOOP ---------------- ####
for(j in seq_along(N)){
  
  current.N <- N[j]
  
  # Define m values
  if(current.N == 2500){
    m.values <- c(20)   # only Vecchia
  } else {
    m.values <- c(20, current.N)  # Vecchia + Full
  }
  
  for(m in m.values){
    
    cat("Running N =", current.N, " m =", m, "\n")
    
    #### ---------- PRE-COMPUTATION ---------- ####
    start_pre <- Sys.time()
    
    data <- sim.data(
      q = q, N = current.N,
      family = family,
      true.beta = true.beta,
      true.Sigma = true.Sigma,
      true.phi = true.phi,
      true.nu = true.nu,
      pred.prop = pred.prop
    )
    
    obs.locs <- data$obs.locs
    N.obs    <- data$N.obs
    Y.obs    <- data$Y.obs
    X.obs    <- data$X.obs
    
    phi <- data$true.phi
    nu  <- data$true.nu
    
    obs.ord <- max_min(obs.locs)
    obs.locs.ord <- obs.locs[obs.ord, , drop = FALSE]
    
    distobs.ord <- rdist(obs.locs.ord)
    diameter <- max(distobs.ord)
    b_phi <- diameter / 3
    
    #### Prediction setup ####
    pred.locs <- data$pred.locs
    N.pred    <- data$N.pred
    
    pred.ord <- max_min(pred.locs)
    ord <- c(obs.ord, pred.ord + nrow(obs.locs))
    
    locs.ord <- rbind(obs.locs, pred.locs)[ord, , drop = FALSE]
    
    NNarray.all <- neighbor_matrix(locs.ord, m)
    distall.nn  <- dist.nn(locs.ord, NNarray.all)
    coeff.all   <- vecchia_coeff(distall.nn, NNarray.all, phi, nu)
    
    X.pred <- data$X.pred
    Y.pred.true <- data$Y.pred.true
    
    X.pred.ord <- X.pred[pred.ord, , drop = FALSE]
    Y.pred.ord.true <- Y.pred.true[pred.ord, , drop = FALSE]
    
    Y.obs.ord <- Y.obs[obs.ord, , drop = FALSE]
    X.obs.ord <- X.obs[obs.ord, , drop = FALSE]
    
    beta  <- data$true.beta
    Sigma <- data$true.Sigma
    W.obs.ord <- data$true.W.obs[obs.ord, , drop = FALSE]
    
    end_pre <- Sys.time()
    
    pre_time <- as.numeric(difftime(end_pre, start_pre, units = "secs"))
    
    #### ---------- MCMC ---------- ####
    out <- mixed.workflow.joint(
      Y.obs.ord, X.obs.ord, obs.locs.ord, m, N.obs, p, q,
      nu, W.obs.ord, beta, Sigma, phi,
      M.prior, V.prior, S.prior, df.prior, b_phi,
      niters, tuning.phi,
      Y.pred.ord.true, X.pred.ord, N.pred,
      distall.nn, coeff.all, NNarray.all,
      family, family.par = NULL, link = NULL
    )
    
    #### ---------- STORE ---------- ####
    runtime.df <- rbind(runtime.df, data.frame(
      N = current.N,
      m = m,
      Method = ifelse(m == current.N, "Full", "Vecchia"),
      Precompute = pre_time,
      Estimation = out$est.time,
      Prediction = out$pred.time,
      Total = pre_time + out$est.time + out$pred.time
    ))
  }
  
  #### ---------- ADD NA ROW FOR FULL (N = 2500) ---------- ####
  if(current.N == 2500){
    
    runtime.df <- rbind(runtime.df, data.frame(
      N = 2500,
      m = 2500,
      Method = "Full",
      Precompute = NA,
      Estimation = NA,
      Prediction = NA,
      Total = NA
    ))
  }
}

save(runtime.df, file = "runtime_results.RData")

#### ---------------- RESHAPE FOR TABLE ---------------- ####

get_val <- function(method, Nval, var){
  val <- runtime.df[runtime.df$Method == method & runtime.df$N == Nval, var]
  if(length(val) == 0 || is.na(val)) return(NA)
  return(round(val, 4))
}

Ns <- c(100, 500, 2500)

runtime.wide <- data.frame(
  N = Ns,
  
  Pre_Vecchia = sapply(Ns, function(n) get_val("Vecchia", n, "Precompute")),
  Pre_Full    = sapply(Ns, function(n) if(n == 2500) NA else get_val("Full", n, "Precompute")),
  
  Est_Vecchia = sapply(Ns, function(n) get_val("Vecchia", n, "Estimation")),
  Est_Full    = sapply(Ns, function(n) if(n == 2500) NA else get_val("Full", n, "Estimation")),
  
  Pred_Vecchia = sapply(Ns, function(n) get_val("Vecchia", n, "Prediction")),
  Pred_Full    = sapply(Ns, function(n) if(n == 2500) NA else get_val("Full", n, "Prediction"))
)

#### ---------------- XTABLE ---------------- ####
library(xtable)

xt <- xtable(
  runtime.wide,
  caption = paste(
    "Wall-clock times (seconds) including pre-computation, estimation, and prediction",
    "for 100 MCMC iterations of our Algorithm 1. Full GP results for $N=2500$ are not feasible due to memory constraints."
  ),
  label = "tab:wall_time",
  align = c("l","c","c","c","c","c","c","c")
)

addtorow <- list()
addtorow$pos <- list(0)

addtorow$command <- paste0(
  "\\hline\n",
  " & \\multicolumn{2}{c}{\\text{Pre-computation}}",
  " & \\multicolumn{2}{c}{\\text{Estimation}}",
  " & \\multicolumn{2}{c}{\\text{Prediction}} \\\\\n",
  "\\cline{2-3} \\cline{4-5} \\cline{6-7}\n",
  " \\text{$N$}",
  " & \\text{Vecchia} & \\text{Full}",
  " & \\text{Vecchia} & \\text{Full}",
  " & \\text{Vecchia} & \\text{Full} \\\\\n",
  "\\hline\n"
)

print(xt,
      include.rownames = FALSE,
      include.colnames = FALSE,  
      sanitize.text.function = identity,
      add.to.row = addtorow,
      hline.after = c(nrow(runtime.wide)))