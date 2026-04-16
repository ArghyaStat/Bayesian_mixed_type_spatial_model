rm(list = ls())

############################
# Libraries
############################
library(this.path)
library(sp)
library(gstat)

setwd(this.path::here())

############################
# Load data
############################
load("data_full.RData", eva.data <- new.env())
full.data.DF <- eva.data$data_DF

############################
# Extract time slice
############################
N.s <- 3503
t.index <- 140

idx <- seq(N.s * (t.index - 1) + 1, N.s * t.index)
data <- full.data.DF[idx, ]

############################
# Locations
############################
obs.locs <- cbind(data$lon, data$lat)

############################
# Responses (IMPORTANT: keep original form)
############################
Y.obs <- cbind(data$BA, data$CNT)

############################
# Residuals (UNCHANGED MODEL)
############################
lm.BA  <- lm(log(1 + Y.obs[,1]) ~ obs.locs)
lm.CNT <- lm(log(1 + Y.obs[,2]) ~ obs.locs)

resid.BA  <- lm.BA$residuals
resid.CNT <- lm.CNT$residuals

############################
# Spatial objects
############################
df.BA <- data.frame(lon = obs.locs[,1],
                    lat = obs.locs[,2],
                    values = resid.BA)
coordinates(df.BA) <- ~ lon + lat

df.CNT <- data.frame(lon = obs.locs[,1],
                     lat = obs.locs[,2],
                     values = resid.CNT)
coordinates(df.CNT) <- ~ lon + lat

############################
# Empirical variograms
############################
emp.BA  <- variogram(values ~ 1, df.BA)
emp.CNT <- variogram(values ~ 1, df.CNT)

############################
# ORIGINAL INITIALIZATION (KEEP THIS)
############################
init.BA <- vgm(model = "Mat",
               psill  = 0.5 * var(resid.BA),
               range  = 0.2 * max(emp.BA$dist),
               nugget = 0.1 * var(resid.BA),
               kappa  = 0.5)

init.CNT <- vgm(model = "Mat",
                psill  = 0.5 * var(resid.CNT),
                range  = 0.2 * max(emp.CNT$dist),
                nugget = 0.1 * var(resid.CNT),
                kappa  = 0.5)

############################
# Fit (UNCHANGED)
############################
fit.BA  <- fit.variogram(emp.BA,  model = init.BA,  fit.kappa = TRUE)
fit.CNT <- fit.variogram(emp.CNT, model = init.CNT, fit.kappa = TRUE)

############################
# Extract range
############################
phi.BA  <- fit.BA$range[2]
phi.CNT <- fit.CNT$range[2]

pdf("Variogram_resid_combined.pdf", width = 10, height = 4)

par(mfrow = c(1,2),
    mar = c(3.2, 3.2, 0.5, 0.5),
    oma = c(0,0,0,0),
    mgp = c(2, 0.7, 0),
    tcl = -0.3)

############################
# BA
############################
plot(emp.BA$dist, emp.BA$gamma,
     xlab = "Distance",
     ylab = "Semivariance",
     pch = 16,
     col = "black",
     cex.lab = 1.2,
     cex.axis = 1.2)

lines(variogramLine(fit.BA, maxdist = max(emp.BA$dist)),
      lwd = 2,
      col = "black")

# ---- ADD THIS (range inside panel) ----
usr <- par("usr")
text(x = usr[2] - 0.05*(usr[2]-usr[1]),
     y = usr[3] + 0.05*(usr[4]-usr[3]),
     labels = bquote(phi == .(round(phi.BA, 4))),
     adj = c(1,0),
     cex = 1)

############################
# CNT
############################
plot(emp.CNT$dist, emp.CNT$gamma,
     xlab = "Distance",
     ylab = "Semivariance",
     pch = 16,
     col = "black",
     cex.lab = 1.2,
     cex.axis = 1.2)

lines(variogramLine(fit.CNT, maxdist = max(emp.CNT$dist)),
      lwd = 2,
      col = "black")

# ---- ADD THIS ----
usr <- par("usr")
text(x = usr[2] - 0.05*(usr[2]-usr[1]),
     y = usr[3] + 0.05*(usr[4]-usr[3]),
     labels = bquote(phi == .(round(phi.CNT, 4))),
     adj = c(1,0),
     cex = 1)

dev.off()