######################################################################################
#########################Data and Parameter Generation(M2PL model)####################
######################################################################################
library(MASS)
set.seed(1818)
N<-10000              ##Number of examinees
J<-20                 ##Number of items
r<-2                  ##Number of dimensions of latent traits
Sigma<-diag(1,r)
theta<-mvrnorm(N,rep(0,r),Sigma)
a<-matrix(rlnorm(J*r,0.3,0.2),nrow=J,ncol=r)
a[1:r,]<-diag(1,r)
b<-matrix(rnorm(J),nrow=1)
b[1:r]<-0
M_b<-matrix(1,nrow=N,ncol=1)%*%b
p<-1/(1+exp(M_b-theta%*%t(a)))
nrep<-25              ##Number of replications
y<-as.list(rep(0,nrep))
set.seed(8888)
for(i in 1:nrep){
  u<-matrix(runif(N*J),nrow=N,ncol=J)
  y[[i]]<-ifelse(p>u,1,0)
}



######################################################################################
####################Polya-gamma algorithm for full data (M2PL model)##################
######################################################################################
library(truncnorm)
library(BayesLogit)
library(MCMCpack)
PGM2PL<-function(u,ntime,nburn){
  t1<-Sys.time()
  N<-nrow(u)
  J<-ncol(u)
  
  Rtheta<-as.list(rep(0,ntime))
  Ra<-as.list(rep(0,ntime))
  Rb<-matrix(0,ntime,J)
  
  theta0<-matrix(0,nrow=N,ncol=r)
  a0<-matrix(1,nrow=J,ncol=r)
  a0[1:r,]<-diag(1,r)
  b0<-matrix(0,nrow=J,ncol=1)
  Sigma0<-diag(1,r)
  M_b0<-matrix(1,nrow=N,ncol=1)%*%t(b0)
  for(i in 1:ntime){
    temp<-theta0%*%t(a0)-M_b0
    w0<-matrix(rpg(N*J,1,temp),N,J,byrow = F)
    
    ##theta, prior N(0,1)
    miu_theta<-rep(0,r)
    sig_theta<-Sigma0
    temp1<-(kronecker(w0,rep(1,r))*kronecker(rep(1,N),t(a0)))%*%a0+kronecker(rep(1,N),solve(sig_theta))
    V_theta<-lapply(c(1:N),function(x) solve(temp1[((x-1)*r+1):(x*r),]))
    temp2<-t(a0)%*%t((u-1/2)+M_b0*w0)+solve(sig_theta)%*%miu_theta%*%matrix(1,nrow=1,ncol=N)
    m_theta<-lapply(c(1:N),function(x) V_theta[[x]]%*%temp2[,x])
    theta0<-matrix(0,N,r)
    for(j in 1:N){
      theta0[j,]<-mvrnorm(1,m_theta[[j]],V_theta[[j]])
    }
    Rtheta[[i]]<-theta0
    
    ##b,prior N(0,10)
    miu_b<-0
    sig_b<-10
    V_b<-1/(colSums(w0)+1/sig_b)
    m_b<-V_b*(colSums((theta0%*%t(a0))*w0-(u-1/2))+miu_b/sig_b)
    b0<-rnorm(J,m_b,sqrt(V_b))
    b0[1:r]<-0
    Rb[i,]<-b0
    M_b0<-matrix(1,nrow=N,ncol=1)%*%b0
    
    ##a,prior N(0,10)
    miu_a<-0
    sig_a<-10
    V_a<-1/(t(w0)%*%theta0^2+1/sig_a)
    a_temp<-lapply(c(1:r),function(d) t(theta0[,d])%*%((u-1/2)+w0*M_b0-w0*(theta0[,-d]%*%t(a0[,-d])))+miu_a/sig_a)
    m_a<-V_a*t(Reduce("rbind",a_temp))
    a0<-matrix(rtruncnorm(1,a=0,b=Inf,mean=t(m_a),sd=t(sqrt(V_a))),J,r,byrow = T)
    a0[1:r,]<-diag(1,r)
    Ra[[i]]<-a0
  }
  t0<- Sys.time()
  time<-difftime(t0,t1,units="hours")
  theta_hat<-Reduce("+",Rtheta[-(1:nburn)])/(ntime-nburn)
  a_hat<-Reduce("+",Ra[-(1:nburn)])/(ntime-nburn)
  b_hat<-apply(Rb[-(1:nburn),],2,mean)
  return(list(theta_hat=theta_hat,a_hat=a_hat,b_hat=b_hat,time=time))
} 

######################################################################################
##################################Simulation (M2PL model)#############################
######################################################################################
ntime<-10000
nburn<-5000
nrep<-length(y)
start <- Sys.time()
results_full<-as.list(rep(0,nrep))
for(x in 1:nrep){
  results_full[[x]]<-PGM2PL(y[[x]],ntime,nburn)
}
end <- Sys.time()
difftime(end,start,units="hours")

#####bias, RMSE and time
Bias_full<-Reduce("+",lapply(results_full,function(x) c(mean(x$theta_hat-theta),mean(x$a_hat-a),mean(x$b_hat-b[1,]))))/nrep
RMSE_full<-sqrt(Reduce("+",lapply(results_full,function(x) c(mean((x$theta_hat-theta)^2),mean((x$a_hat-a)^2),mean((x$b_hat-b[1,])^2))))/nrep)
names(Bias_full)<-c("theta_Bias","a_Bias","b_Bias")
names(RMSE_full)<-c("theta_RMSE","a_RMSE","b_RMSE")

max_time<-max(unlist(lapply(results_full,function(x) x$time)))
max_time
