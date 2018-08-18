#------------------------
# Model 1 - Larval Growth
#------------------------
sink("PopDyProjLarva.txt")
cat("
model{                                           
  for (i in 1:N) {
    L[i] ~ dnorm(L_Exp[i], sig)   
    L_Exp[i] <- Linf * (1.0 - exp(-k * (A[i] - t0)))
    A_ob[i] ~ dnorm(A[i], tau)
    A[i] ~ dunif(0.0001, 0.5)    

# posterior prediction
    L.pred[i] ~ dnorm(L_Exp[i], sig)
    p.value[i] <- step(L.pred[i] - L[i])
  }

# Priors  
  Linf ~ dunif(30, 70)
  k ~ dunif(0.01, 5)
  t0 ~ dunif(-2, 2)
  sig ~ dgamma(0.001, 0.0001)
  tau ~ dgamma(0.001, 0.0001)
}
  ", fill = TRUE)
sink()


#--------------------------------------------
# Model 2 - Larval Growth with heterogeneity
#--------------------------------------------
sink("PopDyProjLarvae.txt")
cat("
model{                                           
  for (i in 1:N) {
    L[i] ~ dnorm(L_Exp[i], sig)   
    L_Exp[i] <-  Linf[i] * (1.0 - exp(-k[i] * (A[i] - t0)))
    A_ob[i] ~ dnorm(A[i], tau)
    A[i] ~ dunif(0.0001, 0.5)    
    
# posterior prediction
    L.pred[i] ~ dnorm(L_Exp[i], sig)
    p.value[i] <- step(L.pred[i] - L[i])
    
    Linf[i] ~ dnorm(Linf_mu,  Linf_tau) I(20, 60)    
    k[i] ~ dnorm(k_mu, k_tau) I(0.1, 5)
  }
    
  Linf_std <- sqrt(1 / Linf_tau)
  k_std <- sqrt(1 / k_tau)
  tau_std <- sqrt(1 / tau)
    
  Linf_mu ~ dnorm(40, 0.001) 
  Linf_tau ~ dgamma(0.01, 0.001)
  k_mu ~ dunif(0.01, 5)
  k_tau ~ dgamma(0.01, 0.001)
  tau ~ dgamma(0.01, 0.001)
  sig ~ dgamma(0.01, 0.001)
  t0 ~ dunif(-2, 2)
}
  ", fill = TRUE)
sink()

#-----------------------
# Model 3 - Adult Growth
#-----------------------
sink("PopDyProjAdult.txt")
cat("
model{
  for (i in 1:N) {              # i indexes individual animal
    for (j in 1:n.caps[i]-1) {  # j indexes recap occasion
      dL[i, j] ~ dnorm(L_Exp[i, j], sig)   
      L_Exp[i, j] <-  (Linf - L[i, j]) * (1.0 - exp(-k * (dt[i, j])))
      L_ob[i, j] ~ dnorm(L[i, j], tau)
      L[i, j+1] ~ dunif(30, 100)
# posterior prediction
      L.pred[i, j] ~ dnorm(L_Exp[i, j], sig)
      p.value[i, j] <- step(L[i, j] + L.pred[i, j] - L[i, j+1])
    }
    L[i, 1] ~ dunif(30, 100)
  }

# Priors
  Linf ~ dnorm(60, 0.0001)
  k ~ dunif(0.01, 1)
  sig ~ dgamma(0.001, 0.0001)
  tau ~ dgamma(0.001, 0.0001)
}
  ", fill = TRUE)
sink()

#--------------------------------------------
# Model 4 - Adult Growth with Heterogeneity
#--------------------------------------------
sink("PopDyProjAdults.txt")
cat("
model{                                           
  for (i in 1:N) {                   # i indexes individual animal
    for (j in 1:n.caps[i]-1) {       # j indexes recap occasion
      dL[i, j] ~ dnorm(L_Exp[i, j], sig)   
      L_Exp[i, j] <-  (Linf[i] - L[i, j]) * (1.0 - exp(-k[i] * (dt[i, j])))
      L_ob[i, j] ~ dnorm(L[i, j], tau)
      L[i,j+1] ~ dunif(30, 100)
    
# posterior prediction
      L.pred[i, j] ~ dnorm(L_Exp[i, j], sig)
      p.value[i, j] <- step(L[i, j] + L.pred[i, j] - L[i, j+1])
    }
    L[i, 1] ~ dunif(30, 100)
    Linf[i] ~ dnorm(Linf_mu,  Linf_tau)     
    k[i] ~ dnorm(k_mu, k_tau) I(0, 1)
  }

# Priors
  Linf_std <- sqrt(1 / Linf_tau)
  k_std <- sqrt(1 / k_tau)
  tau_std <- sqrt(1 / tau)
    
  Linf_mu ~ dnorm(60, 0.001)
  Linf_tau ~ dgamma(0.01, 0.001)
  k_mu ~ dunif(0.01, 1)
  k_tau ~ dgamma(0.01, 0.001)
  tau ~ dgamma(0.01, 0.001)
  sig ~ dgamma(0.01, 0.001)
}
  ", fill = TRUE)
sink()

#--------------------------------------------
# Model 5 - Both with heterogeneity
#--------------------------------------------
sink("PopDyProjFinal.txt")
cat("
model{                                           
  for (i in 1:M) {
    S[i] ~ dnorm(S_Exp[i], sig1)   
    S_Exp[i] <-  a1[i] * (1.0 - exp(-k1[i] * (A[i] - t0)))
    A_ob[i] ~ dnorm(A[i], tau.e)
    A[i] ~ dunif(0.0001, 0.3)

# posterior prediction
    S.pred[i] ~ dnorm(S_Exp[i], sig1)
    q.value[i] <- step(S.pred[i] - S[i])
    
    a1[i] ~ dnorm(a.mu1,  tau.a)
    k1[i] ~ dnorm(k.mu1, tau.k)
  }
    
  for (i in 1:N) {                # i indexes individual animal
    for (j in 1:n.caps[i]-1) {    # j indexes recap occasion
      dL[i, j] ~ dnorm(L_Exp[i, j], sig2)   
      L_Exp[i, j] <-  (a2[i] - L[i, j]) * (1.0 - exp(-k2[i] * (dt[i, j])))
      L_ob[i, j] ~ dnorm(L[i, j], tau.e)
      L[i, j+1] ~ dunif(30, 100)
    
# posterior prediction
      L.pred[i, j] ~ dnorm(L_Exp[i, j], sig2)
      p.value[i, j] <- step(L[i, j] + L.pred[i, j] - L[i, j+1])
    }
        
    a2[i] ~ dnorm(a.mu2,  tau.a)   
    k2[i] ~ dnorm(k.mu2, tau.k)
    L[i, 1] ~ dunif(30, 100)
  }
    
# Priors
  a.mu1 ~ dnorm(40, 0.2)
  a.mu2 ~ dnorm(60, 0.001) 
  s.a ~ dunif(0, 100) 
  k.mu1 ~ dnorm(0, 0.000001)
  k.mu2 ~ dnorm(0, 0.000001)
  s.k ~ dunif(0, 100) 
  s.e ~ dunif(0, 100) 
  sig1 ~ dgamma(0.1, 0.1) 
  sig2 ~ dgamma(0.1, 0.1) 
  
  tau.a <- 1/pow(s.a,2)
	tau.k <- 1/pow(s.k,2)
	tau.e <- 1/pow(s.e,2)

  t0 ~ dunif(-2, 2)
  tmax ~ dunif(0.25, 0.5)
}
  ", fill = TRUE)
sink()

#----------------------------------------------
# Model 6 - Both with heterogeneity and linked
#----------------------------------------------
sink("PopDyProjFinal.txt")
cat("
model{                                           
  for (i in 1:L.n) {  # larval growth (known ages)
    LL[i] <-  Linf.mu1 * (1 - exp(-k1[i] * (t[i] - t0)))    # objective function
    LL.ob[i] ~ dnorm(LL[i], tau.L) I(0, )                   # likelihood function 
    
    t.ob[i] ~ dnorm(t[i], tau.t) I(0.001, t.ob[i])          # likelihood function
    t[i] ~ dunif(0.01, t.m)                                 # prior (with hyper parameters)
        
    k1[i] ~ dnorm(k.mu1, k.tau1) I(2, 4)                    # prior (with hyper parameters)

  L.m[i] <- Linf.mu1 * (1 - exp(-k1[i] * (t.m - t0)))

    # posterior prediction
    L.pred[i] ~ dnorm(LL[i], tau.L)I(0,)         
    q.value[i] <- step(L.pred[i] - LL[i])

    
  
  }
  

  for (i in 1:N) {                                                           # adult growth (unknown ages)
      for (j in 1:n.caps[i]-1) {                                               # j indexes recap occasion
        dL[i, j] <-  (Linf2[i] - L[i, j]) * (1.0 - exp(-k2[i] * (dt[i, j])))   # objective function   
        dL.ob[i, j] ~ dnorm(dL[i, j], tau.a)                                   # likelihood function
        L.ob[i, j] ~ dnorm(L[i, j], tau.L) I(0, )                              # likelihood function

        L[i, j] <- (Linf2[i] - L.meta) * (1.0 - exp(-k2[i] * (t.a[i] - (am[i, j] + t.m)))) + L.meta

# posterior prediction
        dL.pred[i, j] ~ dnorm(dL[i, j], tau.a)
        p.value[i, j] <- step(L[i, j] + dL.pred[i, j] - L[i, j+1])
     } 



    for (j in n.caps[i]:n.caps[i]) {
        L[i,j] <- (Linf2[i] - L.meta) * (1.0 - exp(-k2[i] * (t.a[i] - t.m))) + L.meta
    }    

    t1[i] <- t.m - log((Linf2[i] - L.meta) / Linf2[i]) / -k2[i]

    Linf2[i] ~ dnorm(Linf.mu2,  Linf.tau2) I(0, )               # prior (with hyperparameters)
    k2[i] ~ dnorm(k.mu2, k.tau2) I(0, 2)                        # prior (with hyperparameters)
    min.a[i] <- sum(dt[i,1:n.caps[i]-1])
    t.a[i] ~ dunif(min.a[i], 15)    
}      


for (i in 1:N) {
      L.m[i+L.n] <- Linf2[i] * (1 - exp(-k2[i] * (t.m - t0)))
}

L.meta <- mean(L.m[])

for(i in 1:761){
  s.meta[i] ~ dnorm(L.meta, tau.m)
}

# Priors
  
  Linf.std2 <- sqrt(1 / Linf.tau2)
  k.std1 <- sqrt(1 / k.tau1)
  k.std2 <- sqrt(1 / k.tau2)
    
  L.std <- sqrt(1 / tau.L)   
  a.std <- sqrt(1 / tau.a)   
  t.std <- sqrt(1 / tau.t)   
  m.std <- sqrt(1 / tau.m)   

  Linf.mu1 ~ dunif(30, 50)
  Linf.mu2 ~ dunif(50, 80)

  Linf.tau2 ~ dgamma(0.01, 0.01) 
    
  k.mu1 ~ dunif(1, 4)
  k.mu2 ~ dunif(0, 2) 

  k.tau1 ~ dgamma(0.01, 0.01) 
  k.tau2 ~ dgamma(0.01, 0.01)
    
  tau.L ~ dgamma(0.01, 0.01) 
  tau.a ~ dgamma(0.01, 0.01) 
  tau.t ~ dgamma(0.01, 0.01) 
  tau.m ~ dgamma(0.01, 0.01) 
  
  t0 ~ dunif(-2, 2)
  t.m ~ dunif(0.2, 0.5)
  t1.mu <- mean(t1[])
}
  ", fill = TRUE)
sink()


#----------------------------
# New model parameterization
#---------------------------

sink("VBGE2018.txt")
cat("
    model{

#LARVAL PROCESS
for( i in 1 : N.L) {
    pred.LL[i]<-Linf1*(1-exp(-k1*(age.p[i]-t0)))   #predicted length
    LL[i] ~ dnorm(pred.LL[i],prec.LL)
}  
    pred.m<-Linf1*(1-exp(-k1*(t.m-t0)))
    
    # Priors
    Linf1~dnorm(pred.m,prec.m)  
    k1~dunif(3,8)  
    t0 ~ dunif(-0.5,0.5)
    prec.LL<-1/(sd.LL*sd.LL)
    sd.LL~dunif(0,1)
    
    
    #METAMORPH PROCESS
for( i in 1 : N.M) {
    LM[i] ~ dnorm(pred.m, prec.m)
}
    
#Priors
prec.m<-1/(sd.m*sd.m)
sd.m~dunif(0,5)
    
    
    #ADULT PROCESS
for( i in 1 : N.A) {
    pred.L0[i]<-Linf2 - (Linf2 - pred.m)*(exp(-k2*(t.a[i]-t.m)))   #predicted length at release
    L0[i]~dnorm(pred.L0[i],prec.L)
    
    for(j in 1:J[i]){
      pred.L[i,j]<-Linf2 - (Linf2 - pred.m)*(exp(-k2*(t.a[i]+dt[i,j]-t.m))) #length at jth recapture
      L[i,j]~dnorm(pred.L[i,j],prec.L)
}}
    
    #PRIORS
for( i in 1 : N.A) {
    t.a[i]~dlnorm(log.centre,prec.t) I(t.m,30)
} 
    
Linf2~dnorm(60,.002)  
k2~dunif(0,5)  
    
t.m ~ dunif(0.2, 0.4)
    
prec.L<-1/(sd.L*sd.L)
sd.L~dunif(0,2)
    
centre~dunif(0,200)
log.centre<-log(centre)
prec.t~dgamma(.01,.01)
}
  ", fill = TRUE)
sink()