library(ggplot2)
library(plyr)
library(reshape)

VGBE = load("VGBE_output.RData")

# extract relevant posteriors
pl = data.frame(VGBE$sims.list)
pl = pl[(1:14)]
colnames(pl) = c("Linf.mu1", "Linf.mu2", "Linf.std2", "k.mu1", "k.std1", "k.mu2", "k.std2",
                 "t.m", "L.m", "L.std", "a.std", "t.std", "t0", "t1.mu")

#-----------------------------------------
# Sample posterior plots for growth rate K
#-----------------------------------------
pp = pl[c(4, 6)]

k.l = density(pp[, 1])
k.a = density(pp[, 2])

plot(k.l, ylim = c(0, 10), xlim = c(0, 4), lwd = 2, lty = 1, 
     ylab = "", xlab = "", main = "Posterior Distributions k", axes = FALSE)
par(new = TRUE)
plot(k.a, ylim = c(0, 10), xlim = c(0, 4), lwd = 2, lty = 2, ylab = "", xlab = "", main = "", axes = FALSE)
axis(1)
axis(2)
par(las = 0)
mtext("k", side = 1, line = 2.5, cex = 1.5)
mtext("Density", side = 2, line = 3, cex = 1.8)


dat.melt <- melt(pp)

ggplot()+
  geom_density(data = dat.melt, aes(value, fill = variable, colour = variable), 
               alpha=I(.1), size =1)  +
  facet_wrap(~ variable, scales = "free", nrow = 3)


#---------------
# VBG curve plot
#---------------
# extract relevant data
dp = data.frame(VGBE$summary)
dd = dp[c(1:14), ]     # parameters of interest
la = dp[15:425, ]      # predicted age larvae
ll = dp[426:836, ]     # predicted length larvae            
al = dp[1248:2329, ]   # predicted lengths adults
aa = dp[3788:4140, ]   # predicted age adults

# Means of model parameters
Linf.pred.l = dd$mean[1]
Linf.pred.a = dd$mean[2]
k.pred.l = dd$mean[4]
k.pred.a = dd$mean[6]
t.m = dd$mean[8]
L.m = dd$mean[9]
t0.pred.l <-dd$mean[13]
t0.pred.a = dd$mean[14]  

# Posterior 95% intervals
k.a.u = k.pred.a + (2 * dd$mean[7]) 
k.a.l = k.pred.a - (2 * dd$mean[7])
L.a.u = Linf.pred.a + (2 * dd$mean[3])
L.a.l = Linf.pred.a - (2 * dd$mean[3])
k.l.u = k.pred.l + (2 * dd$mean[5])
k.l.l = k.pred.l - (2 * dd$mean[5])


# larval curve with 95% intervals
plot(t.m, L.m, cex = 1, ylim = c(0, 70), pch = 16, col = "red",
     ylab = "", xlab = "", main = "VB Growth", axes = FALSE, xlim = c(-.5, 8))
curve(Linf.pred.l * (1 - exp(-k.pred.l * (x - t0.pred.l))), 0, t.m, add = T, lwd = 2)
axis(1)
axis(2)
par(las = 0)
mtext("Age", side = 1, line = 2.5, cex = 1.5)
mtext("Length", side = 2, line = 3, cex = 1.8)

curve(Linf.pred.l * (1 - exp(-k.l.u * (x - t0.pred.l))), 0, t.m, add = T, lty = "dashed", lwd = 1.5)
curve(Linf.pred.l * (1 - exp(-k.l.l * (x - t0.pred.l))), 0, t.m, add = T, lty = "dashed", lwd = 1.5)

# adult curve with 95% intervals
f = function(x){Linf.pred.a * (1 - exp(-k.pred.a * (x - t0.pred.a)))}
curve(f(x), t.m, 8, lwd = 2, add = T)
axis(1)
axis(2)
par(las = 0)
mtext("Age", side = 1, line = 2.5, cex = 1.5)
mtext("Length", side = 2, line = 3, cex = 1.8)

curve(L.a.u * (1 - exp(-k.a.u * (x - t0.pred.a))), t.m, 8, add = T, lty = "dashed", lwd = 1.5)
curve(L.a.l * (1 - exp(-k.a.l * (x - t0.pred.a))), t.m, 8, add = T, lty = "dashed", lwd = 1.5)


#------------------------------
# age predictions distribution
#------------------------------

# function for just getting L[i,1]
Ls = al[1, ]
count = 1

for (i in 1:length(aa[, 1])) {
  count = count + n[i]
  Ls = rbind(Ls, al[count, ]) 
}
Ls = Ls[-354, ]

# add predicted ages from lengths on final curve
dci = cbind(aa, Ls)
dci = dci[order(dci[, 1]), ]
points(dci[, 5], dci[, 8], pch = 16, cex = 0.5, col = "grey")

la.pr = cbind(ll, la)
la.pr = la.pr[order(la.pr[, 8]), ]
points(la.pr[, 8], la.pr[, 1], pch = 16, cex = 0.2, col = "grey")

# add age and size at metamorphosis to plot
points(t.m, L.m, cex = 1, pch = 16, col = "red")

# histogram of predicted ages (mean and mode)
hist(da$X50., breaks = 50)
hist(da$mean, breaks = 50)

da[1, ]

