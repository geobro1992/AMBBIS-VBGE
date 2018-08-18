library(ggplot2)
library(plyr)
library(reshape)
library(coda)

load("VGBE_Output_6_5_18.RData")

# extract relevant posteriors
pl = data.frame(VGBE$sims.list)

# model fit bayesian p values
stps = pl[(426:5425)]

lstps = stps[(1:411)]
astps = stps[(412:5000)]

plot(1:411, colMeans(lstps))
plot(1:4589, colMeans(astps))

mean(colMeans(lstps))
mean(colMeans(astps))

# check convergence (need multiple chains)
pl = data.frame(VGBE$sims.list[1:14])
pl = as.mcmc.list(pl)
gelman.plot(VGBE, confidence = 0.95, transform=FALSE, autoburnin=TRUE, multivariate=TRUE)
traceplot(VGBE)
#-----------------------------------------
# Sample posterior plots for growth rate K
#-----------------------------------------
pp = pl[c(11:383)]

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
aa = dp[1566:1918, ]   # predicted age adults

# Means of model parameters
Linf.pred.l = dd$mean[1]
Linf.l.l = dd$X2.5.[1]
Linf.l.u = dd$X97.5.[1]

Linf.pred.a = dd$mean[7]
L.a.l = dd$X2.5.[7]
L.a.u = dd$X97.5.[7]

k.pred.l = dd$mean[2]
k.l.l = dd$X2.5.[2]
k.l.u = dd$X97.5.[2]

k.pred.a = dd$mean[8]
k.a.l = dd$X2.5.[8]
k.a.u = dd$X97.5.[8]

t.m = dd$mean[9]
t.m.l = dd$X2.5.[9]
t.m.u = dd$X97.5.[9]

L.m = dd$mean[6]
L.m.l = dd$mean[6] - 2*dd$mean[5]
L.m.u = dd$mean[6] + 2*dd$mean[5]

t0 = dd$mean[4]


# larval curve with 95% intervals
plot(0, 0, cex = 0.01, ylim = c(0, 80), pch = 16, col = "grey",
     ylab = "", xlab = "", main = "VB Growth", axes = FALSE, xlim = c(-.5, 5))
x1 = curve(Linf.pred.l * (1 - exp(-k.pred.l * (x - t0))), 0, t.m.l, add = T, lwd = 2)
axis(1)
axis(2)
par(las = 0)
mtext("Age", side = 1, line = 2.5, cex = 1.5)
mtext("Length", side = 2, line = 3, cex = 1.8)

x1u = curve(L.a.l * (1 - exp(-k.l.u * (x - t0))), 0, t.m.l, add = T, lty = "dashed", lwd = 1.5)
x1l = curve(L.m.l * (1 - exp(-k.l.l * (x - t0))), 0, t.m.l, add = T, lty = "dashed", lwd = 1.5)

# adult curve with 95% intervals
f = function(x){Linf.pred.a - (Linf.pred.a - L.m) * exp(-k.pred.a * (x - t0))}
x2 = curve(f(x), t.m.u, 15, lwd = 2, add = T)
axis(1)
axis(2)
par(las = 0)
mtext("Age", side = 1, line = 2.5, cex = 1.5)
mtext("Length", side = 2, line = 3, cex = 1.8)

x2u = curve(L.a.u - (L.a.u - L.m.u) * exp(-k.a.u * (x - t0)), t.m.u, 15, add = T, lty = "dashed", lwd = 1.5)
x2l = curve(L.a.l - (L.a.l - L.m.l) * exp(-k.a.l * (x - t0)), t.m.u, 15, add = T, lty = "dashed", lwd = 1.5)

abline(v=t.m.l)
abline(v=t.m.u)

dat = cbind(rbind(as.data.frame(x1), as.data.frame(x2)), 
            rbind(as.data.frame(x1l), as.data.frame(x2l)),
            rbind(as.data.frame(x1u), as.data.frame(x2u)))
names(dat) = c("x", "y", "xl", "yl", "xu", "yu")

p1 <- ggplot(dat, aes(x, y))+
  geom_line(size=1,linetype = 2)+
  geom_ribbon(data=dat,aes(ymin=yl,ymax=yu),alpha=0.3)+
  geom_vline(xintercept = t.m) +
  labs(x = "Age") +
  labs(y = "SVL / mm") +
  theme(axis.title.x = element_text(face="bold", size=20),
           axis.text.x  = element_text(size=15),
        axis.title.y = element_text(face="bold", size=20),
        axis.text.y  = element_text(size=15))

p1 + scale_x_log10() 

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
dci = cbind(aa, bugs.data$L.ob[,1])
dci = dci[order(dci[, 10]), ]
points(dci[, 5], dci[, 10], pch = 16, cex = 0.5, col = "grey")

la.pr = cbind(ll, la)
la.pr = la.pr[order(la.pr[, 8]), ]
points(la.pr[, 8], la.pr[, 1], pch = 16, cex = 0.2, col = "grey")

# add age and size at metamorphosis to plot
points(t.m, L.m, cex = 1, pch = 16, col = "red")

# histogram of predicted ages (mean and mode)
hist(dci[,1], breaks = 50)
hist(dci[,2], breaks = 50) 
hist(dci[,7], breaks = 50)

aa = aa[order(aa[,1]), ]
p <- ggplot(aa, aes(x = 1:353, y = aa$mean, 
                     ymin = aa$X25., ymax = aa$X75.)) + 
  geom_pointrange() + 
  theme_bw() + 
  labs(x = "Individual", y = "Age") +
  coord_flip() 

summary(dci)

p <- ggplot(pp, aes(x = 1:353, y = aa$mean, 
                    ymin = aa$X25., ymax = aa$X75.)) + 
  geom_pointrange() + 
  theme_bw() + 
  labs(x = "Individual", y = "Age") +
  coord_flip() 

yfc = vector()
for( i in 1:373){
  yfc[i] = format(ags[[i]][1,2],"%Y")
}

tas = apply(pp, 2, mean)

a10 = tas - as.numeric(yfc) + 2010
a11 = a10 + 1
a12 = a10 + 2
a13 = a10 + 3
a14 = a10 + 4
a15 = a10 + 5
a16 = a10 + 6

par(mfrow = c(1,1))
hist(a10[which(a10 > 0 & a10 < 16)], breaks = 18, 
     xlab = "Age", main = "Age Distribution in 2010")
hist(a11[which(a11 > 0 & a11 < 16)], breaks = 15)
hist(a12[which(a12 > 0 & a12 < 16)], breaks = 15)
hist(a13[which(a13 > 0 & a13 < 16)], breaks = 15)
hist(a14[which(a14 > 0 & a14 < 16)], breaks = 15)
hist(a15[which(a15 > 0 & a15 < 16)], breaks = 15)
hist(a16[which(a16 > 0 & a16 < 16)], breaks = 15)


VGBE$summary

#----------------------
# Size by year and age histograms
library(RODBC)
library(tidyverse)
library(lubridate)
library(forcats)

# Connect to Minnow server
channel <- odbcConnectAccess2007("//minnow.cc.vt.edu/cnre2/Eglin/projects/FlatwoodsSalamander/DriftFence/DriftFenceDatabase/Driftfence Database_GB_29Nov17_FRONT_END.accdb")

channel2 <- odbcConnectAccess2007("//minnow.cc.vt.edu/cnre2/Eglin/projects/FlatwoodsSalamander/LarvalSampling/LarvalDatabase/LarvalDatabase_GB_2Apr18_FRONT_END.accdb")

# Read in capture table
captures <- sqlFetch(channel,"Flatwoods Capture") %>%
  tbl_df() %>%
  glimpse()

size.freq <- captures %>%
  mutate(seas = (Date - 15552000)) %>%
  mutate(year = year(seas)) %>%
  filter(PondNumber != 53) %>%
  filter(year != 2018) %>% 
  filter(year != 2009) %>%
  mutate(Age = fct_recode(Age, 'Y' = 'J')) %>%
  mutate(Age = fct_recode(Age, 'A' = 'U')) %>%
  mutate(Age = fct_recode(Age, 'A' = 'Y')) %>%
  filter(Age == "M") 

size.freq$year = as.factor(size.freq$year)
levels(size.freq$year) = c("10-11", "11-12", "12-13", "13-14", "14-15", "15-16", "16-17", "17-18")

# read in larvae table
dips <- sqlFetch(channel2,"FlatwoodsLarvalCapture") %>%
  tbl_df() %>%
  rename(Survey = Survey_ID) %>%
  glimpse()

survs = sqlFetch(channel2,"Surveys") %>%
  tbl_df() %>%
  mutate(seas = (DATE - 15552000)) %>%
  rename(Survey = SURVEY_ID) %>%
  mutate(year = year(seas)) %>%
  mutate(month = month(DATE)) %>%
glimpse()

x = survs[, c(1,3,48,49)]
y = dips[, c(2,6,7)]

s.freq = merge(x, y, by = "Survey",all=F) %>% 
  filter(year > 2009) %>% 
  filter(month != 11) %>% 
  filter(month != 5) %>% 
  filter(SVL != 0)

s.freq$year = as.factor(s.freq$year)
s.freq$month = as.factor(s.freq$month)
s.freq$month = relevel(s.freq$month, "12")
levels(s.freq$month) = c("Dec", "Jan", "Feb", "Mar", "Apr")
#---------------------------------------------  
# adult size histogram by season

pdf(file = "size_freq_hist.pdf", 6, 15)
ggplot(size.freq, aes(x = SVL_1stMeasurement)) +
  geom_histogram(bins = 70, col = "grey", aes(fill = Age)) +
  theme_bw() +
  facet_wrap( ~ year, scales = 'fixed', nrow = 8) + 
  theme(legend.background = element_rect('transparent'),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'gray90'),
        strip.text = element_text(size = 10),
        panel.background = element_rect('gray80')) +
  scale_fill_brewer(palette = 'YlOrBr') +
  labs(x = 'Size (mm)',
       y = 'Frequency Density')
dev.off()

#-----------------------------------
# larvae sizes by date by season

pdf(file = "larvae_growth.pdf", 6, 15)
ggplot(s.freq, aes(x = SVL)) +
  geom_histogram(bins = 70, col = "grey") +
  theme_bw() +
  facet_wrap( ~ month, scales = 'fixed', nrow = 8) + 
  theme(legend.background = element_rect('transparent'),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'gray90'),
        strip.text = element_text(size = 10),
        panel.background = element_rect('gray80')) +
  labs(x = 'SVL',
       y = 'Frequency')
dev.off()




