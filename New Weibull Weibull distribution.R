#_________________________________________________________________________________________________________
#                                  Upload the following packages
#_________________________________________________________________________________________________________
require(gamlss)
require(numDeriv)

#_________________________________________________________________________________________________________
#                                  Next Weibull implementation in gamlss
#_________________________________________________________________________________________________________
Weibull_<- function (mu.link = "log", sigma.link="log")
{
  mstats <- checklink("mu.link", "Weibull", substitute(mu.link), 
                       c("1/mu^2", "log", "identity", "own"))
  dstats <- checklink("sigma.link", "Weibull", substitute(sigma.link), 
                      c("inverse", "log", "identity", "own"))
  structure(
    list(family = c("Weibull_", "Weibull"),
         parameters = list(mu=TRUE, sigma=TRUE), 
         nopar = 2, 
         type = "Continuous",
         mu.link = as.character(substitute(mu.link)),  
         sigma.link = as.character(substitute(sigma.link)), 
         
         mu.linkfun = mstats$linkfun, 
         sigma.linkfun = dstats$linkfun, 
         
         mu.linkinv = mstats$linkinv, 
         sigma.linkinv = dstats$linkinv,
         
         mu.dr = mstats$mu.eta, 
         sigma.dr = dstats$mu.eta, 
#_________________________________________________________________________________________________________
#                                  First, Second and Cross Numerical Derivatives 
#_________________________________________________________________________________________________________
dldm = function(y,mu,sigma){
           lpdf<-function(t,x,sigma){log(dauxiWeibull_(t,x,sigma))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,method='simple')
           dldm
         },
         
         d2ldm2 = function(y,mu,sigma){
           lpdf<-function(t,x,sigma){log(dauxiWeibull_(t,x,sigma))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,method='simple')
           d2ldm2 <- -dldm * dldm
           d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15) 
           d2ldm2 
         },     
         
         dldd = function(y,mu,sigma){
           lpdf<-function(t,mu,x){log(dauxiWeibull_(t,mu,x))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,method='simple')
           dldd
         } ,
         
         d2ldd2 = function(y,mu,sigma){
           lpdf<-function(t,mu,x){log(dauxiWeibull_(t,mu,x))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,method='simple')
           d2ldd2 <- -dldd*dldd
           d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
           d2ldd2
         },
         
         d2ldmdd = function(y,mu,sigma){
           lpdf<-function(t,x,sigma){log(dauxiWeibull_(t,x,sigma))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,method='simple')
           lpdf<-function(t,mu,x){log(dauxiWeibull_(t,mu,x))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma)
           d2ldmdd = -(dldm * dldd)
           d2ldmdd<-ifelse(is.na(d2ldmdd)==TRUE,0,d2ldmdd)
           d2ldmdd                 
         },

   G.dev.incr  = function(y,mu,sigma) 
         { 
           -2*dWeibull_(y,mu,sigma,log=TRUE)
         } ,                     
         rqres = expression(   
           rqres(pfun="pWeibull_", type="Continuous", y=y, mu=mu, sigma=sigma)) ,
         
         mu.initial = expression(mu <- rep(1, length(y))), 
         sigma.initial = expression(sigma <- rep(1, length(y))), 
         
         mu.valid = function(mu) all(mu> 0) , 
         sigma.valid = function(sigma)all(sigma>0),
         
         y.valid = function(y) all(y>0)),
    class = c("gamlss.family","family"))
}
#_________________________________________________________________________________________________________
#                                  Weibull probability density function
#_________________________________________________________________________________________________________
dWeibull_<-function(x, mu =1, sigma =1, log = FALSE){
  
  if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
  
  fy1 <-sigma*mu^(sigma)*x^(sigma-1)*exp(-(mu*x)^(sigma))
  
  if(log==FALSE) fy<-fy1 else fy<-log(fy1)
  fy
}    
#_________________________________________________________________________________________________________
#                                  Weibull cumulative density function 
#_________________________________________________________________________________________________________
pWeibull_<-function(q, mu =1, sigma = 1,lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
  
  if (any(q < 0))  stop(paste("q must be positive", "\n", ""))          
  
  cdf1 <-(1-exp(-(mu*q)^(sigma)))
  
  if(lower.tail==TRUE) cdf<-cdf1 else cdf<- 1-cdf1
  if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
  cdf
}
#_________________________________________________________________________________________________________
#                                  Weibull quantile function 
#_________________________________________________________________________________________________________
qWeibull_<-  function(p, mu=1, sigma=1,lower.tail = TRUE, log.p = FALSE){   
  if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
  
  if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
  if (log.p==TRUE) p <- exp(p) else p <- p
  if (lower.tail==TRUE) p <- p else p <- 1-p
  
  q<- ((-log(-(p-1)))^(1/sigma))/mu
  q
}
#_________________________________________________________________________________________________________
#                                  Weibull random generating function
#_________________________________________________________________________________________________________
rWeibull_<- function(n, mu=0.1, sigma=0.1){
  if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
  
  if (any(n < 0))  stop(paste("n must be positive", "\n", "")) 
  uni<- runif(n = n,0,1)
  
  r<-qWeibull_(uni,mu=mu,sigma=sigma)
  r
}
#_________________________________________________________________________________________________________
#                                  Weibull auxiliar function for derivates 
#_________________________________________________________________________________________________________
dauxiWeibull_<-function(t,mu,sigma){ 
  pdfWei<- dWeibull_(t, mu=mu, sigma=sigma)
  fy1 <-pdfWei
  fy1}
#_________________________________________________________________________________________________________
#                                  End of Weibull implemented in gamlss
#_________________________________________________________________________________________________________

#_________________________________________________________________________________________________________
#                                  Next New Weibull Weibull (NWW) implementation in gamlss
#_________________________________________________________________________________________________________
NWW<- function (mu.link = "log", sigma.link="log", nu.link = "identity", tau.link = "identity")
{
  mstats <- checklink("mu.link", "New Weibull Weibull", substitute(mu.link), 
                         c("inverse", "log", "identity", "logit"))
  dstats <- checklink("sigma.link", "New Weibull Weibull", substitute(sigma.link), 
                      c("inverse", "log", "identity", "own"))
  vstats <- checklink(   "nu.link", "New Weibull Weibull", substitute(nu.link),    
                         c("1/nu^2", "log", "identity", "own"))
  tstats <- checklink(  "tau.link", "New Weibull Weibull", substitute(tau.link),   
                        c("1/tau^2", "log", "identity", "own")) 
  structure(
    list(family = c("NWW", "New Weibull Weibull"),
         parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE, tau=TRUE), 
         nopar = 4, 
         type = "Continuous",
         mu.link = as.character(substitute(mu.link)),  
         sigma.link = as.character(substitute(sigma.link)), 
         nu.link = as.character(substitute(nu.link)), 
         tau.link = as.character(substitute(tau.link)), 
         mu.linkfun = mstats$linkfun, 
         sigma.linkfun = dstats$linkfun, 
         nu.linkfun = vstats$linkfun,
         tau.linkfun = tstats$linkfun,  
         mu.linkinv = mstats$linkinv, 
         sigma.linkinv = dstats$linkinv,
         nu.linkinv = vstats$linkinv,
         tau.linkinv = tstats$linkinv, 
         mu.dr = mstats$mu.eta, 
         sigma.dr = dstats$mu.eta, 
         nu.dr = vstats$mu.eta,
         tau.dr = tstats$mu.eta, 
#_________________________________________________________________________________________________________
#                                  First, Second and Cross Numerical Derivatives 
#_________________________________________________________________________________________________________

         dldm = function(y,mu,sigma,nu,tau){ 
           lpdf<-function(t,x,sigma,nu,tau){log(dauxiNWW(t,x,sigma,nu,tau))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,tau=tau,method='simple')
           dldm
         },
         d2ldm2 = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,x,sigma,nu,tau){log(dauxiNWW(t,x,sigma,nu,tau))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,tau=tau,method='simple')
           d2ldm2 <- -dldm * dldm
           d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15) 
           d2ldm2 
         },     
         dldd = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,mu,x,nu,tau){log(dauxiNWW(t,mu,x,nu,tau))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,tau=tau,method='simple')
           dldd
         } ,
         d2ldd2 = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,mu,x,nu,tau){log(dauxiNWW(t,mu,x,nu,tau))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,tau=tau,method='simple')
           d2ldd2 <- -dldd*dldd
           d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
           d2ldd2
         },   
         dldv = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,mu,sigma,x,tau){log(dauxiNWW(t,mu,sigma,x,tau))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,tau=tau,method='simple')
           dldv
         },
         d2ldv2 = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,mu,sigma,x,tau){log(dauxiNWW(t,mu,sigma,x,tau))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,tau=tau,method='simple')
           d2ldv2<- -dldv * dldv
           d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2,-1e-15)                   
           d2ldv2
         },
         dldt = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,mu,sigma,nu,x){log(dauxiNWW(t,mu,sigma,nu,x))}
           dldt<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,nu=nu,x=tau)
           dldt
         } ,
         d2ldt2 = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,mu,sigma,nu,x){log(dauxiNWW(t,mu,sigma,nu,x))}
           dldt<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,nu=nu,x=tau,method='simple')
           d2ldt2<- -dldt * dldt
           d2ldt2 <- ifelse(d2ldt2 < -1e-15, d2ldt2,-1e-15)                                    
           d2ldt2
         },
         d2ldmdd = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,x,sigma,nu,tau){log(dauxiNWW(t,x,sigma,nu,tau))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,tau=tau,method='simple')
           lpdf<-function(t,mu,x,nu,tau){log(dauxiNWW(t,mu,x,nu,tau))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,tau=tau)
           d2ldmdd = -(dldm * dldd)
           d2ldmdd<-ifelse(is.na(d2ldmdd)==TRUE,0,d2ldmdd)
           d2ldmdd                 
         },
         d2ldmdv = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,x,sigma,nu,tau){log(dauxiNWW(t,x,sigma,nu,tau))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,tau=tau,method='simple')
           lpdf<-function(t,mu,sigma,x,tau){log(dauxiNWW(t,mu,sigma,x,tau))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,tau=tau)
           d2ldmdv = -(dldm * dldv)
           d2ldmdv				
         },
         
         d2ldmdt = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,x,sigma,nu,tau){log(dauxiNWW(t,x,sigma,nu,tau))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,tau=tau,method='simple')
           lpdf<-function(t,mu,sigma,nu,x){log(dauxiNWW(t,mu,sigma,nu,x))}
           dldt<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,nu=nu,x=tau)
           d2ldmdt <- -(dldm*dldt)
           d2ldmdt
         },
         
         d2ldddv = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,mu,x,nu,tau){log(dauxiNWW(t,mu,x,nu,tau))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,tau=tau,method='simple')
           lpdf<-function(t,mu,sigma,x,tau){log(dauxiNWW(t,mu,sigma,x,tau))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,tau=tau)
           d2ldddv = -(dldd * dldv)
           d2ldddv	
         },
         d2ldddt = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,mu,x,nu,tau){log(dauxiNWW(t,mu,x,nu,tau))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,tau=tau)
           lpdf<-function(t,mu,sigma,nu,x){log(dauxiNWW(t,mu,sigma,nu,x))}
           dldt<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,nu=nu,x=tau)
           d2ldddt <- -(dldd*dldt) 
           d2ldddt 
         },
         d2ldvdt = function(y,mu,sigma,nu,tau){
           lpdf<-function(t,mu,sigma,x,tau){log(dauxiNWW(t,mu,sigma,x,tau))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,tau=tau)
           lpdf<-function(t,mu,sigma,nu,x){log(dauxiNWW(t,mu,sigma,nu,x))}
           dldt<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,nu=nu,x=tau)
           d2ldvdt <- -(dldv*dldt) 
           d2ldvdt 
         },
        
 G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
         { 
           -2*dNWW(y,mu,sigma,nu,tau,log=TRUE)
         } ,                     
         rqres = expression(   
           rqres(pfun="pNWW", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)) ,
         
         mu.initial = expression( mu <- rep(1, length(y))),
         sigma.initial = expression( sigma <-rep(0.5, length(y))),
         nu.initial = expression(nu <- rep(0.8, length(y))), 
         tau.initial = expression(tau <-rep(1, length(y))), 
         
         mu.valid = function(mu) all(mu> 0), 
         sigma.valid = function(sigma) all(sigma > 0),
         nu.valid = function(nu) all(nu > 0), 
         tau.valid = function(tau) all(tau > 0), 
         y.valid = function(y) all(y>0)),
    class = c("gamlss.family","family"))
}
#_________________________________________________________________________________________________________
#                                  NWW probability density function
#_________________________________________________________________________________________________________
dNWW <- function(x, mu=1, sigma=1, nu=1, tau=1, log = FALSE){
  if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))  
  if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
  if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))
  
  pdfsn <- dWeibull_(x,mu=mu,sigma=sigma)
  cdfsn <- pWeibull_(x,mu=mu,sigma=sigma)
  
  fy1 <- ((nu*tau*pdfsn)/(cdfsn))*((-log(cdfsn))^(tau-1))*(exp(-nu*(-log(cdfsn))^(tau)))
  
  if(log==FALSE) fy<-fy1 else fy<-log(fy1)
  fy
}    
#_________________________________________________________________________________________________________
#                                  NWW cumulative density function 
#_________________________________________________________________________________________________________
pNWW <- function(q, mu=0.5, sigma=1, nu=1, tau=2, lower.tail = TRUE, log.p = FALSE){  
  if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))  
  if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
  if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))
  if (any(q < 0))  stop(paste("q must be positive", "\n", ""))          
  
  cdfsn<-pWeibull_(q,mu=mu,sigma=sigma)
  
  cdf1 <- exp(-nu*(-log(cdfsn))^(tau))
  
  if(lower.tail==TRUE) cdf<-cdf1 else cdf<- 1-cdf1
  if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
  cdf
}
#_________________________________________________________________________________________________________
#                                  NWW quantile function 
#_________________________________________________________________________________________________________
qNWW <-  function(p, mu=0.5, sigma=1, nu=1, tau=2, lower.tail = TRUE, log.p = FALSE){   
  if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))  
  if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
  if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))
  
  if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
  if (log.p==TRUE) p <- exp(p) else p <- p
  if (lower.tail==TRUE) p <- p else p <- 1-p
  
  u <- exp(-(log(p)/(-nu))^(1/tau))
  
  q <- qWeibull_(u, mu=mu , sigma=sigma)
  q
}
#_________________________________________________________________________________________________________
#                                  NWW random generating function
#_________________________________________________________________________________________________________
rNWW <- function(n, mu=0.5, sigma=1, nu=1, tau=2){
  if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
  if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))  
  if (any(n < 0))  stop(paste("n must be positive", "\n", "")) 
  
  uni<- runif(n = n,0,1)
  
  r <- qNWW(uni,mu =mu, sigma =sigma, nu=nu, tau=tau)
  r
}
#_________________________________________________________________________________________________________
#                                  NWW auxiliar function for derivates 
#_________________________________________________________________________________________________________
dauxiNWW <- function(t,mu,sigma,nu,tau){ 
  pdfsn <- dWeibull_(t,mu=mu,sigma=sigma)
  cdfsn <- pWeibull_(t,mu=mu,sigma=sigma)
  
  fy1 <- ((nu*tau*pdfsn)/(cdfsn))*((-log(cdfsn))^(tau-1))*(exp(-nu*(-log(cdfsn))^(tau)))
  fy1}

#_________________________________________________________________________________________________________
#                                  End of NWW implemented in gamlss
#_________________________________________________________________________________________________________
