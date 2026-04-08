# Traceplots

par(mfrow = c(1,1))

pdf("traceplot_phi.pdf", width=5 ,height=5)
trace.phi <- plot.ts(out$phi.samples, ylab = "phi", main = "Traceplot of phi")
dev.off()

# acfplots

pdf("acf_phi.pdf", width=5 ,height=5)
acf.phi <- acf(out$phi.samples, main = "ACF plot of phi", lag.max = 200)
dev.off()




labels.beta <- list()
for (i in 1:p) {
  for (j in 1:q) {
    labels.beta <- c(labels.beta, 
                     list(paste0('beta (', i, ',', j, ')')))
  }
}


pdf("traceplot_beta.pdf", width=7.5 ,height=6)
par(mfrow=c(p,q))

for (i in 1:p) {
  for (j in 1:q) {
    plot(1:niters, sapply(out$beta.samples, function(x) x[i, j]), 
         type='l', col=1, xlab='Iteration', ylab='beta',
         main=labels.beta[(i-1)*q + j])
    
  }
}
dev.off()

pdf("acf_beta.pdf", width= 7.5 ,height=6)
par(mfrow=c(p,q))

for (i in 1:p) {
  for (j in 1:q) {
    
    acf(sapply(out$beta.samples, function(x) x[i, j]), 
        main = labels.beta[(i-1)*q + j], lag.max = 100)
    
  }
}
dev.off()

# Output Analysis of Sigma

labels.Sigma <- list()
for (i in 1:q) {
  for (j in 1:q) {
    labels.Sigma <- c(labels.Sigma, 
                      list(paste0('Sigma (', i, ',', j, ')')))
  }
}


pdf("traceplot_Sigma.pdf", width = 7.5 ,height = 6)
par(mfrow=c(q,q))

for (i in 1:q) {
  for (j in 1:q) {
    plot(1:niters, sapply(out$Sigma.samples, function(x) x[i, j]), 
         type='l', col=1, xlab='Iteration', ylab='Sigma',
         main=labels.Sigma[(i-1)*q + j])
    
  }
}
dev.off()


pdf("acf_Sigma.pdf", width = 7.5 ,height = 6)
par(mfrow=c(q,q))

for (i in 1:q) {
  for (j in 1:q) {
    
    acf(sapply(out$Sigma.samples, function(x) x[i, j]), 
        main=labels.Sigma[(i-1)*q + j],
        lag.max = 200)
    
  }
}
dev.off()