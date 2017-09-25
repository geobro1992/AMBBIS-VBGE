rm(list = ls())
library(R2WinBUGS)
library(ggplot2)
library(plyr)
library(reshape)

#------------
# Model 1
#-----------

df = read.csv("larvae_SVL.csv", sep = ",")

A_ob = df$Age
L = df$SVL
dat = na.omit(cbind(A_ob, L, df$Group))

bugs.data = list(
  L = dat[, 2],
  A_ob = dat[, 1] / 365,
  N = length(dat[, 1])
)

dat1 = data.frame(dat)
ggplot(dat1, aes(x = A_ob / 365, y = L)) +
  geom_point(aes(colour = as.factor(dat1[, 3])))

# Initial values
#inits <- function(){list(Linf = 90, k = 0.5, tau.eps = 10, m.eps = 10)}  
inits <- function(){list(Linf = runif(411, 39, 41), k = runif(411, 1, 1.2), sig = 5, tau = 5, t0 = 0)}  

# Parameters monitored
parameters <- c("Linf_mu", "Linf_std", "k_mu", "k_std", "tau_std", "t0")

# MCMC settings
ni <- 100000
nt <- 20
nb <- 20000
nc <- 3

# Call WinBUGS from R (BRT 40 min)
js.super <- bugs(bugs.data, inits, parameters, "PopDyProjLarvae.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 debug = TRUE, bugs.directory = "C:/Program Files (x86)/WinBUGS14", working.directory = getwd())

#-----------------------------------------------------------------------
# Posterior Plots
pl = data.frame(js.super$sims.list)

pl = pl[-7]
colnames(pl) = c("Linf", "Linf Sigma", "k", "k Sigma", "Tau", "t0")

dat.melt <- melt(pl)

ggplot() +
  geom_density(data = dat.melt, aes(value, fill = variable, colour = variable), 
               alpha = I(.1), size = 1) +
  facet_wrap(~ variable, scales = "free", nrow = 3)


k.l = density(pl[,1])
k.a = density(pl[,3])

plot(k.l, ylim = c(0,10), xlim = c(0, 4), lwd = 2, lty = 1, 
     ylab = "", xlab = "", main = "Posterior Distributions k", axes = FALSE)

par(new = TRUE)

plot(k.a, ylim = c(0, 10), xlim = c(0, 4), lwd = 2, lty = 2, ylab = "", xlab = "", main = "", axes = FALSE)
axis(1)
axis(2)
par(las = 0)
mtext("k", side = 1, line = 2.5, cex = 1.5)
mtext("Density", side = 2, line = 3, cex = 1.8)

# VBG curve plot
k.pred.l <- js.super$mean$k_mu
Linf.pred.l <- js.super$mean$Linf_mu

plot(curve(Linf.pred.l * (1 - exp(-k.pred.l * x)), 0, 10), pch = 16, cex = 0.5, ylim = c(0, 70),
     ylab = "", xlab = "", main = "VB Growth", axes = FALSE)
axis(1)
axis(2)
par(las = 0)
mtext("Age", side = 1, line = 2.5, cex = 1.5)
mtext("Length", side = 2, line = 3, cex = 1.8)

#-------------
# Model 2
#-----------
df = read.csv("FW_main.csv", sep = ",")
df = df[, c("MASTER.ID", "Date", "Sex", "Weight", "TL", "SVL_1stMeasurement")]

dl = split(df, df$MASTER.ID, drop = T)
names(dl) = NULL

# create empty vectors for all the VBGE data
TLt = matrix(NA, 353, 50)
SVLt = matrix(NA, 353, 50)
Wt = matrix(NA, 353, 50)
dt = matrix(NA, 353, 50)
dSVL = matrix(NA, 353, 50)
dW = matrix(NA, 353, 50)
dTL = matrix(NA, 353, 50)
ID = vector()
n = vector()
count = 0
# number of caps for each ID  
n.caps = sapply(dl, NROW)

# for each unique ID, get dt, dL, and L
for (i in 1:length(dl)) {
  if (n.caps[i] == 1) {
  } else {
    count = count + 1
      for (j in 2:length(dl[[i]][[1]])) {
      
        Wt[count, j-1] = dl[[i]][j-1, 4]
        Wt[count, j] = dl[[i]][j, 4]
      
        TLt[count, j-1] = dl[[i]][j-1, 5]
        TLt[count, j] = dl[[i]][j, 5]
      
        SVLt[count, j-1] = dl[[i]][j-1, 6]
        SVLt[count, j] = dl[[i]][j, 6]
      
        dW[count, j-1] = dl[[i]][j, 4] - dl[[i]][j-1, 4]
        dTL[count, j-1] = dl[[i]][j, 5] - dl[[i]][j-1, 5]
        dSVL[count, j-1] = dl[[i]][j, 6] - dl[[i]][j-1, 6]
      
        dt[count, j-1] = as.Date(dl[[i]][j, 2], format = "%m/%d/%Y") - as.Date(dl[[i]][j-1, 2], format = "%m/%d/%Y")
      } # j
      ID = append(ID, dl[[i]][1, 1])
      n = append(n, n.caps[i])
  } 
}   # i


# Put all together in a list
bugs.data = list(
  n.caps = n,
  L_ob = SVLt,
  dL = dSVL,
  dt = dt / 365,
  N = length(n)
)


# Initial values
inits <- function(){list(Linf = runif(353, 59, 61), k = runif(353, 0.4, 0.6), 
     tau = 5, sig = 5)}
     
parameters <- c("Linf_mu", "Linf_std", "k_mu", "k_std", "tau_std")

# MCMC settings
ni <- 100000
nt <- 20
nb <- 20000
nc <- 3

# Call WinBUGS from R (without heterogeneity  : PopDyProjAdult.txt)
js.super <- bugs(bugs.data, inits, parameters, "PopDyProjAdults.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 debug = TRUE, bugs.directory = "C:/Program Files (x86)/WinBUGS14", working.directory = getwd())

# Posterior Plots 

pl = data.frame(js.super$sims.list)
pl = pl[-6]
colnames(pl) = c("Linf", "Linf Sigma", "k", "k Sigma", "Tau")

dat.melt <- melt(pl)

ggplot() +
  geom_density(data = dat.melt, aes(value, fill = variable, colour = variable), 
               alpha = I(.1), size = 1) +
  facet_wrap(~ variable, scales = "free", nrow = 3)

# VBG curve plot
k.pred.a <- js.super$mean$k_mu
Linf.pred.a <- js.super$mean$Linf_mu

plot(curve(Linf.pred.a * (1 - exp(-k.pred.a * x)), 0, 10), pch = 16, cex = 0.5, ylim = c(0, 70),
     ylab = "", xlab = "", main = "VB Growth", axes = FALSE)

f = function(x){Linf.pred.l * (1 - exp(-k.pred.l * x))}
points(seq(from = 0, to = 10, length.out = 100), f(seq(from = 0, to = 10, length.out = 100)), pch = 17, cex = 0.5)
axis(1)
axis(2)
par(las = 0)
mtext("Age", side = 1, line = 2.5, cex = 1.5)
mtext("Length", side = 2, line = 3, cex = 1.8)
legend(6, 30, legend = c("Adult", "Larvae"), pch = c(16, 17))

#-----------------------------------------------------------------------
# Model 5 -Full Model
#-----------------------------------------------------------------------

#------------
# Model 1
#-----------

df = read.csv("larvae_SVL.csv", sep = ",")

A_ob = df$Age
S = df$SVL
dat = na.omit(cbind(A_ob, S, df$Group))
dat[which(dat[, 3] == 1), 1] = dat[which(dat[, 3] == 1), 1] - 25  

df = read.csv("FW_main.csv", sep = ",")
df = df[, c("MASTER.ID", "Date", "Sex", "Weight", "TL", "SVL_1stMeasurement")]

dl = split(df, df$MASTER.ID, drop = T)
names(dl) = NULL

# create empty vectors for all the VBGE data
TLt = matrix(NA, 353, 50)
SVLt = matrix(NA, 353, 50)
Wt = matrix(NA, 353, 50)
dt = matrix(NA, 353, 50)
dSVL = matrix(NA, 353, 50)
dW = matrix(NA, 353, 50)
dTL = matrix(NA, 353, 50)
ID = vector()
n = vector()
count = 0
# number of caps for each ID  
n.caps = sapply(dl, NROW)

# for each unique ID, get dt, dL, and L
for (i in 1:length(dl)) {
  if (n.caps[i] == 1) {
  } else {
      count = count + 1
        for (j in 2:length(dl[[i]][[1]])) {
      
          Wt[count, j-1] = dl[[i]][j-1, 4]
          Wt[count, j] = dl[[i]][j, 4]
      
          TLt[count, j-1] = dl[[i]][j-1, 5]
          TLt[count, j] = dl[[i]][j, 5]
      
          SVLt[count, j-1] = dl[[i]][j-1, 6]
          SVLt[count, j] = dl[[i]][j, 6]
      
          dW[count, j-1] = dl[[i]][j, 4] - dl[[i]][j-1, 4]
          dTL[count, j-1] = dl[[i]][j, 5] - dl[[i]][j-1, 5]
          dSVL[count, j-1] = dl[[i]][j, 6] - dl[[i]][j-1, 6]
      
          dt[count, j-1] = as.Date(dl[[i]][j, 2], format = "%m/%d/%Y") - as.Date(dl[[i]][j-1, 2], format = "%m/%d/%Y")
        } # j
    ID = append(ID, dl[[i]][1, 1])
    n = append(n, n.caps[i])
  } 
}   # i

# using first and last captures like this makes it a lot easier with the model subscripts (vectors rather than matrices)
fco = vector()
lco = vector()

for (i in 1:length(n)) {
  fco[i] = sum(n[1:i-1]) + 1
  lco[i] = sum(n[1:i])
} 

# Put all together in a list
bugs.data = list(
# Larval data
  LL.ob = dat[, 2],
  t.ob = dat[, 1] / 365,
  L.n = length(dat[, 1]),
# Adult data  
  n.caps = n,
  L.ob = SVLt,
  dL.ob = dSVL,
  dt = dt / 365,
  N = length(n)
)



#---------------------
# Model setup
#---------------------

# Initial values
inits <- function(){list(Linf.mu1 = 60, k1 = runif(411, 3.2, 3.5), k.tau1 = 4, t0 = 0,
                         Linf2 = runif(353, 50, 70), k2 = runif(353, 0.4, 0.6), Linf.tau2 = .1, k.tau2 = 2,
                         t = runif(411, 0.19, 0.21), t.a = runif(353, 2, 8), tau.L = 4, tau.a = 4, tau.t = 4,
                         t.m = 0.4)}  

# Parameters monitored
parameters <- c( "Linf.mu1", "Linf.mu2", "Linf.std2", "k.mu1", "k.std1", "k.mu2", "k.std2", 
                 "t.m", "L.m", "L.std", "a.std", "t.std", "t0", "t1.mu",
                 "t", "L.pred", "q.value", "L", "dL.pred", "p.value", "t.a", "t1")
 
# MCMC settings
ni <- 1000000
nt <- 500
nb <- 10000
nc <- 1

# Call WinBUGS from R (BRT 40 min)
VGBE <- bugs(bugs.data, inits, parameters, "PopDyProjFinal.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 debug = TRUE, bugs.directory = "C:/Program Files (x86)/WinBUGS14", working.directory = getwd())

save(VGBE, file = "VGBE_output.RData")




