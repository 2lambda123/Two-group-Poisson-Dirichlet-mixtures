# Example
# credo che soffra dove ci sono molte obs.
# rimedio: partire con piu' full sweeps e poi lasciar correre lo Spl-Mer
# piccolo cane che si morde la codea: uso SM per tante obs, ma con tante obs SM puo' faticare
# Warmstart?

prior_par11 <- list(m0=0,   s0=1, 
                    m1=3,   V1=3, 
                    a1= 1, b1=1, 
                    a_rho=1,  b_rho=9,
                    theta0=1, theta1=1, 
                    sigma0=0.75, sigma1=0.1, # parameters of the process
                    a_0 = 5, b_0 = .2,
                    kappaNLP=3, s1=2)

nsim=1000
BI=2
datadata = c(rnorm(20),rnorm(20,-3),rnorm(20,5))
zz = GIBBS_SAMPLER_BNPtesting_v5(NSIM = nsim , 
                                 burn_in = BI, 
                            thinning = 1, 
                            y= datadata, 
                            prior_par = prior_par11,
                            StartingSWEEP = 10,
                            SWEEPevery = 100,
                            densities = T, 
                            verbose = T,
                            verbose_step = 10, 
                            sed = i*100,
                            SM = .5,warmstart = T,
                            optThresh = .44,
                            batch = 100) 


plot(datadata,zz$Y)

hist(zz$Y,freq=F)
lines(zz$PostPred[,1],zz$PostPred[,2])
lines(zz$PostPred[,1],zz$PostPred[,3])
lines(zz$PostPred[,1],zz$PostPred[,4])
plot(zz$m1Con,type = 'l')
plot(zz$RhoCon,type = 'l')
plot(zz$mppi~zz$Y)




# Example

prior_par11 <- list(m0=0,   s0=1, 
                    m1=3,   V1=3, 
                    a1= 1, b1=1, 
                    a_rho=1,  b_rho=9,
                    theta0=1, theta1=1, 
                    sigma0=0.75, sigma1=0.1, # parameters of the process
                    a_0 = 5, b_0 = .2,
                    kappaNLP=3, s1=2)

nsim=1000
BI=2
datadata = c(rnorm(500),rnorm(50,-3),rnorm(50,5))
zz = GIBBS_SAMPLER_BNPtesting_v5(NSIM = nsim , 
                                 burn_in = BI, 
                            thinning = 1, 
                            y= datadata, 
                            prior_par = prior_par11,
                            StartingSWEEP = 100,
                            SWEEPevery = 100,
                            densities = T, 
                            verbose = T,
                            verbose_step = 10, 
                            sed = i*100,
                            SM = .5,
                            ##############################
                            warmstart = T,
                            #############################
                            optThresh = .44,
                            batch = 100) 


plot(datadata,zz$Y)

hist(zz$Y,freq=F)
lines(zz$PostPred[,1],zz$PostPred[,2])
lines(zz$PostPred[,1],zz$PostPred[,3])
lines(zz$PostPred[,1],zz$PostPred[,4])
plot(zz$m1Con,type = 'l')
plot(zz$RhoCon,type = 'l')
plot(zz$mppi~zz$Y)
