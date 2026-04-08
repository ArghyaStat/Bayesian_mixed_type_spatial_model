rm(list = ls())

library(rlist)
library(this.path)

setwd(this.path::here())

load("bgp_jointESS_N500.RData")
load("bgp_compESS_N500.RData")
load("bgp_pCN_N500.RData")
load("bgp_rwmh_N500.RData")

chain.joint.ess <- out.joint.ess$log.post.full
chain.comp.ess  <- out.comp.ess$log.post.full
chain.pCN       <- out.pCN$log.post.full
chain.rwmh      <- out.rwmh$log.post.full

essW.joint.ess <- out.joint.ess$W.obs.ord.ess
essW.comp.ess  <- out.comp.ess$W.obs.ord.ess
essW.pCN       <- out.pCN$W.obs.ord.ess
essW.rwmh      <- out.rwmh$W.obs.ord.ess


### Figure 7 (left panel)

cols <- c("#CC79A7", "#0072B2", "#009E73", "#D55E00")


pdf("traceplot_posterior.pdf", width = 9, height = 7.5)

par(mar = c(3.5, 4.5, 1, 1), mgp = c(2.5,0.8,0), tcl = -0.3)

# -------- global y-limits --------
ymin <- min(chain.joint.ess, chain.comp.ess,
            chain.pCN, chain.rwmh)

ymax <- max(chain.joint.ess, chain.comp.ess,
            chain.pCN, chain.rwmh)

# -------- colors with transparency --------
cols.alpha <- c(
  adjustcolor(cols[1], alpha.f = 0.5),
  adjustcolor(cols[2], alpha.f = 0.5),
  adjustcolor(cols[3], alpha.f = 0.5),
  adjustcolor(cols[4], alpha.f = 0.5)
)

# -------- plot RAW chains --------
plot(chain.joint.ess, type = "l", col = cols.alpha[1], lwd = 1,
     ylim = c(ymin, ymax),
     xlab = "Iterations", ylab = "Log-posterior (unnormalized)",
     cex.lab = 2, cex.axis = 1.4)

lines(chain.comp.ess, col = cols.alpha[2], lwd = 1)
lines(chain.pCN,      col = cols.alpha[3], lwd = 1)
lines(chain.rwmh,     col = cols.alpha[4], lwd = 1)

# -------- legend --------
legend("bottomleft",
       legend = c("Joint Elliptical","Comp Elliptical","Comp pCN","RWMH"),
       col = cols,
       lwd = 2, cex = 1.6, bty = "n")

dev.off()


ess.joint <- as.vector(essW.joint.ess)
ess.comp  <- as.vector(essW.comp.ess)
ess.pCN   <- as.vector(essW.pCN)
ess.rwmh  <- as.vector(essW.rwmh)

### Figure 7 (right panel)

pdf("ess_boxplot_log.pdf", width = 9, height = 7.5)

par(mar = c(3.5, 4.5, 1, 1), mgp = c(2.5,0.8,0), tcl = -0.3)

boxplot(log(ess.joint), log(ess.comp),
        log(ess.pCN), log(ess.rwmh),
        col = cols,
        names = c("Joint Elliptical","Comp Elliptical","Comp pCN","RWMH"),
        ylab = expression(log*"(ESS)"),
        cex.lab = 2, cex.axis = 1.4)

dev.off()