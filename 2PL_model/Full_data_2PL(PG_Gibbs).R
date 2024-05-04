######################################################################################
#########################Data and Parameter Generation(2PL model)#####################
######################################################################################
N<-10000                    ##Number of examinees
set.seed(1818)              
J<-20                       ##Number of items
theta<-matrix(rnorm(N),ncol=1)
a<-matrix(rlnorm(J,0.3,0.2),nrow=1)
b<-matrix(rnorm(J),nrow=1)
#the matrix form of a, b and theta
M_a<-matrix(1,nrow=N,ncol=1)%*%a
M_b<-matrix(1,nrow=N,ncol=1)%*%b
M_theta<-theta%*%matrix(1,nrow=1,ncol=J)
p<-1/(1+exp(-(M_theta-M_b)*M_a))
nrep<-25                   ##Number of replications
y<-as.list(rep(0,nrep))
set.seed(666)
for(i in 1:nrep){
  u<-matrix(runif(N*J),nrow=N,ncol=J)
  y[[i]]<-ifelse(p>u,1,0)
}


######################################################################################
####################Polya-Gamma algorithm for full data (2PL model)###################
######################################################################################
library(truncnorm)
library(BayesLogit)
PG2PL<-function(u,ntime,nburn){
  t1<-Sys.time()
  N<-nrow(u)
  J<-ncol(u)
  Rtheta<-matrix(0,N,ntime)
  Ra<-matrix(0,ntime,J)
  Rb<-matrix(0,ntime,J)
  theta0<-matrix(0,nrow=N,ncol=1)
  a0<-matrix(1,nrow=1,ncol=J)
  b0<-matrix(0,nrow=1,ncol=J)
  
  M_theta0<-theta0%*%matrix(1,nrow=1,ncol=J)
  M_a0<-matrix(1,nrow=N,ncol=1)%*%a0
  M_b0<-matrix(1,nrow=N,ncol=1)%*%b0
  for(i in 1:ntime){
    temp<-M_a0*abs(M_theta0-M_b0)
    w0<-matrix(rpg(N*J,1,temp),N,J,byrow = F)
    
    ##theta, prior N(0,1)
    miu_theta<-0
    sig_theta<-1
    V_theta<-1/(rowSums(M_a0^2*w0)+1/sig_theta)
    m_theta<-V_theta*(rowSums(M_a0^2*w0*M_b0+M_a0*(u-1/2))+miu_theta/sig_theta)
    theta0<-rnorm(N,m_theta,sqrt(V_theta))
    Rtheta[,i]<-theta0
    M_theta0<-theta0%*%matrix(1,nrow=1,ncol=J)
    
    ##a, prior N(0,10)
    miu_a<-0
    sig_a<-10
    V_a<-1/(colSums((M_theta0-M_b0)^2*w0)+1/sig_a)
    m_a<-V_a*(colSums((M_theta0-M_b0)*(u-1/2))+miu_a/sig_a)
    a0<-rtruncnorm(J,a=0,b=Inf,mean=m_a,sd=sqrt(V_a))
    Ra[i,]<-a0
    M_a0<-matrix(1,nrow=N,ncol=1)%*%a0
    
    ##b, prior N(0,10)
    miu_b<-0
    sig_b<-10
    V_b<-1/(colSums(M_a0^2*w0)+1/sig_b)
    m_b<-V_b*(colSums(M_a0^2*w0*M_theta0-M_a0*(u-1/2))+miu_b/sig_b)
    b0<-rnorm(J,m_b,sqrt(V_b))
    Rb[i,]<-b0
    M_b0<-matrix(1,nrow=N,ncol=1)%*%b0
  }
  t0<- Sys.time()
  time<-difftime(t0,t1,units="hours")
  theta_hat<-apply(Rtheta[,-(1:nburn)],1,mean)
  a_hat<-apply(Ra[-(1:nburn),],2,mean)
  b_hat<-apply(Rb[-(1:nburn),],2,mean)
  return(list(theta_hat=theta_hat,a_hat=a_hat,b_hat=b_hat,time=time))
} 


######################################################################################
##################################Simulation(2PL model)###############################
######################################################################################
ntime<-10000
nburn<-5000
nrep<-length(y)
start <- Sys.time()
results_full<-as.list(rep(0,nrep))
for(x in 1:nrep){
  results_full[[x]]<-PG2PL(y[[x]],ntime,nburn)
}
end <- Sys.time()
difftime(end,start,units="hours")

#####bias, RMSE and time
Bias_full<-Reduce("+",lapply(results_full,function(x) c(mean(x$theta_hat-theta),mean(x$a_hat-a[1,]),mean(x$b_hat-b[1,]))))/nrep
RMSE_full<-sqrt(Reduce("+",lapply(results_full,function(x) c(mean((x$theta_hat-theta)^2),mean((x$a_hat-a[1,])^2),mean((x$b_hat-b[1,])^2))))/nrep)
names(Bias_full)<-c("theta_Bias","a_Bias","b_Bias")
names(RMSE_full)<-c("theta_RMSE","a_RMSE","b_RMSE")

max_time<-max(unlist(lapply(results_full,function(x) x$time)))
max_time

