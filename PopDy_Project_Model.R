setwd("C:/Users/boa10gb/Documents/R/VB Growth")

#------------------------
# Model 1 - Larval Growth
#------------------------
sink("PopDyProjLarva.txt")
cat("
    model{                                           
    for (i in 1:N)     {
    
    L[i] ~ dnorm(L_Exp[i], sig)   
    L_Exp[i] <-  Linf*(1.0 - exp(-k*(A[i]-t0)))
    A_ob[i] ~ dnorm(A[i], tau)
    A[i] ~ dunif(0.0001, 0.5)    


    # posterior prediction
    L.pred[i] ~ dnorm(L_Exp[i], sig)
    p.value[i] <- step(L.pred[i] - L[i])

    
    
  }
    
    Linf ~ dnorm(60, 0.0001)
    k ~ dunif(0.01, 1)
    t0 ~ dunif(-2, 2)
    sig ~ dgamma(0.001, 0.0001)
    tau ~ dgamma(0.001, 0.0001)
}
    
    
    ",fill=TRUE)
sink()


#--------------------------------------------
# Model 2 - Larval Growth with heterogeneity
#--------------------------------------------
sink("PopDyProjLarvae.txt")
cat("
    model{                                           
    for (i in 1:N)     {
    
    L[i] ~ dnorm(L_Exp[i], sig)   
    L_Exp[i] <-  Linf[i]*(1.0 - exp(-k[i]*(A[i]-t0)))
    A_ob[i] ~ dnorm(A[i], tau)
    A[i] ~ dunif(0.0001, 0.5)    
    
    
    # posterior prediction
    L.pred[i] ~ dnorm(L_Exp[i], sig)
    p.value[i] <- step(L.pred[i] - L[i])
    
    
     Linf[i] ~ dnorm(Linf_mu,  Linf_tau) I(20,60)    
    k[i] ~ dnorm(k_mu, k_tau) I(0.1,5)
    
    }
    
    Linf_std <- sqrt(1/Linf_tau)
    k_std <- sqrt(1/k_tau)
    tau_std <- sqrt(1/tau)
    
    Linf_mu ~ dnorm(40, 0.001) 
    Linf_tau ~ dgamma(0.01, 0.001)
    k_mu ~ dunif(0.01, 5)
    k_tau ~ dgamma(0.01, 0.001)
    tau ~ dgamma(0.01, 0.001)
    sig ~ dgamma(0.01, 0.001)
    t0 ~ dunif(-2, 2)
    }
    
    
    ",fill=TRUE)
sink()

#-----------------------
# Model 3 - Adult Growth
#-----------------------
sink("PopDyProjAdult.txt")
cat("

model{
    for (i in 1:N)     {
        for (j in 1:n.caps[i]-1)   {    # j indexes recap occasion
           
            dL[i,j] ~ dnorm(L_Exp[i,j], sig)   
            L_Exp[i,j] <-  (Linf - L[i,j]) *(1.0 - exp(-k*(dt[i,j])))
            L_ob[i,j] ~ dnorm(L[i,j], tau)
            L[i,j+1] ~ dunif(30, 100)

            # posterior prediction
            L.pred[i,j] ~ dnorm(L_Exp[i,j], sig)
            p.value[i,j] <- step(L[i,j] + L.pred[i,j] - L[i,j+1])
        }
            L[i,1] ~ dunif(30, 100)

  }
  Linf ~ dnorm(60, 0.0001)
    k ~ dunif(0.01, 1)
    sig ~ dgamma(0.001, 0.0001)
    tau ~ dgamma(0.001, 0.0001)
}
    
    
    ",fill=TRUE)
sink()

#--------------------------------------------
# Model 4 - Adult Growth with Heterogeneity
#--------------------------------------------
sink("PopDyProjAdults.txt")
cat("
model{                                           
    for (i in 1:N)     {                   # i indexes individual animal
        for (j in 1:n.caps[i]-1)   {    # j indexes recap occasion
           
    dL[i,j] ~ dnorm(L_Exp[i,j], sig)   
    L_Exp[i,j] <-  (Linf[i] - L[i,j]) *(1.0 - exp(-k[i]*(dt[i,j])))
    L_ob[i,j] ~ dnorm(L[i,j], tau)
    L[i,j+1] ~ dunif(30, 100)
    
    # posterior prediction
    L.pred[i,j] ~ dnorm(L_Exp[i,j], sig)
    p.value[i,j] <- step(L[i,j] + L.pred[i,j] - L[i,j+1])

        }

    L[i,1] ~ dunif(30, 100)
    Linf[i] ~ dnorm(Linf_mu,  Linf_tau)     
    k[i] ~ dnorm(k_mu, k_tau) I(0,1)

  }

  Linf_std <- sqrt(1/Linf_tau)
  k_std <- sqrt(1/k_tau)
  tau_std <- sqrt(1/tau)
    
  Linf_mu ~ dnorm(60, 0.001)
  Linf_tau ~ dgamma(0.01, 0.001)
  k_mu ~ dunif(0.01, 1)
  k_tau ~ dgamma(0.01, 0.001)
  tau ~ dgamma(0.01, 0.001)
  sig ~ dgamma(0.01, 0.001)
}
    
    
    ",fill=TRUE)
sink()


#--------------------------------------------
# Model 5 - Both with heterogeneity
#--------------------------------------------
sink("PopDyProjFinal.txt")
cat("
    model{                                           
    for (i in 1:M)     {
          S[i] ~ dnorm(S_Exp[i], sig1)   
          S_Exp[i] <-  a1[i]*(1.0 - exp(-k1[i]*(A[i]-t0)))
          A_ob[i] ~ dnorm(A[i], tau.e)
          A[i] ~ dunif(0.0001, 0.3)

          # posterior prediction
          S.pred[i] ~ dnorm(S_Exp[i], sig1)
          q.value[i] <- step(S.pred[i] - S[i])
    
          a1[i] ~ dnorm(a.mu1,  tau.a)
          k1[i] ~ dnorm(k.mu1, tau.k)
    
        }
    
    for (i in 1:N)     {            # i indexes individual animal
    for (j in 1:n.caps[i]-1)   {    # j indexes recap occasion
    
    dL[i,j] ~ dnorm(L_Exp[i,j], sig2)   
    L_Exp[i,j] <-  (a2[i] - L[i,j]) *(1.0 - exp(-k2[i]*(dt[i,j])))
    L_ob[i,j] ~ dnorm(L[i,j], tau.e)
    L[i,j+1] ~ dunif(30, 100)
    
    # posterior prediction
    L.pred[i,j] ~ dnorm(L_Exp[i,j], sig2)
    p.value[i,j] <- step(L[i,j] + L.pred[i,j] - L[i,j+1])
    
    }
        
      a2[i] ~ dnorm(a.mu2,  tau.a)   
      
      k2[i] ~ dnorm(k.mu2, tau.k)
      
      L[i,1] ~ dunif(30, 100)
      
    
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
    
    tau.a <- 1/pow(s.a,2)		# next 3 lines convert sd to precision
	  tau.k <- 1/pow(s.k,2)
	  tau.e <- 1/pow(s.e,2)
    
    t0 ~ dunif(-2, 2)
    tmax ~ dunif(0.25, 0.5)
    }
    
    ",fill=TRUE)
sink()



#----------------------------------------------
# Model 6 - Both with heterogeneity and linked
#----------------------------------------------
sink("PopDyProjFinal.txt")
cat("
    model{                                           
    for (i in 1:L.n)     {                              # larval growth (known ages)

    LL[i] <-  Linf.mu1 * (1-exp(-k1[i]*(t[i]-t0)))    # objective function
    LL.ob[i] ~ dnorm(LL[i], tau.L)I(0,)                    # likelihood function 
    t.ob[i] ~ dnorm(t[i], tau.t)I(0,1)                      # likelihood function
    
        t[i] ~ dunif(0.01, t.m)                      # prior (with hyper parameters)
        
        k1[i] ~ dnorm(k.mu1, k.tau1)I(2,4)                  # prior (with hyper parameters)

    # posterior prediction
    L.pred[i] ~ dnorm(LL[i], tau.L)I(0,)         
    q.value[i] <- step(L.pred[i] - LL[i])

    
  
  }
    L.m <- Linf.mu1 * (1-exp(-k.mu1*(t.m-t0)))
  
    for (i in 1:N)     {                              # adult growth (unknown ages)
      for (j in 1:n.caps[i]-1)   {                    # j indexes recap occasion
    
      dL[i,j] <-  (Linf2[i] - L[i,j]) *(1.0 - exp(-k2[i]*(dt[i,j])))   # objective function   
      dL.ob[i,j] ~ dnorm(dL[i,j], tau.a)                               # likelihood function
      L.ob[i,j] ~ dnorm(L[i,j], tau.L)I(0,)                                 # likelihood function

        L[i,j+1] ~ dunif(40, 70)                                       # prior


    # posterior prediction
    dL.pred[i,j] ~ dnorm(dL[i,j], tau.a)
    p.value[i,j] <- step(L[i,j] + dL.pred[i,j] - L[i,j+1])

    
      }
           
            L[i,1] <- (Linf2[i] - L.m) *(1.0 - exp(-k2[i]*(t.a[i]-t.m))) + L.m
            
            t.a[i] ~ dlnorm(1,0.5)I(t.m,12)

            t1[i] <- t.m - log((Linf2[i]-L.m)/Linf2[i])/-k2[i]

            Linf2[i] ~ dnorm(Linf.mu2,  Linf.tau2)I(0,)              # prior (with hyperparameters)
            k2[i] ~ dnorm(k.mu2, k.tau2)I(0,2)                        # prior (with hyperparameters)
  }    


  # Priors
  
    Linf.std2 <- sqrt(1/Linf.tau2)
    
    k.std1 <- sqrt(1/k.tau1)
    k.std2 <- sqrt(1/k.tau2)
    
    L.std <- sqrt(1/tau.L)   
    a.std <- sqrt(1/tau.a)   
    t.std <- sqrt(1/tau.t)   

    Linf.mu1 ~ dunif(30, 70)
    Linf.mu2 ~ dunif(50, 70)


    Linf.tau2 ~ dgamma(0.01, 0.01) 
    
    k.mu1 ~ dunif(2, 4)
    k.mu2 ~ dlnorm(0, 10)I(0,1.5) 


    k.tau1 ~ dgamma(0.01, 0.01) 
    k.tau2 ~ dgamma(0.01, 0.01)
    
    tau.L ~ dgamma(0.01, 0.01) 
    tau.a ~ dgamma(0.01, 0.01) 
    tau.t ~ dgamma(0.01, 0.01) 
    
    t0 ~ dunif(-2, 2)
    t.m ~ dunif(0.3, 0.5)
    t1.mu <- mean(t1[])
  }
    
    ",fill=TRUE)
sink()



#  L.m <- mean(Linf1[])*(1.0 - exp(mean(-k1[])*(t.m-t0)))
# t.m[i] <- log((L.m[i]-Linf1[i])/Linf1[i])/k1[i] + t0
#   t.m ~ dunif(0.2, 0.4)


# L.a <- Linf.mu2*(1.0 - exp(-k.mu2*(t.m-t0)))

#  L.meta <- mean(L.m[]) - mean(L.a[])
# L.meta ~ dnorm(L.m, tau.L)
