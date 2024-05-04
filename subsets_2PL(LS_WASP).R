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

#################################################################################################
###############################LS-WASP algorithm for 2PL model###################################
#################################################################################################

###Polya-Gamma algorithm (parallel version for each subset) using Modified Likelihood method
PG2PL_group_Par<-function(x){
  t1<-Sys.time()
  u<-y_subset[[x]]
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
    V_a<-1/(k*colSums((M_theta0-M_b0)^2*w0)+1/sig_a)
    m_a<-V_a*(k*colSums((M_theta0-M_b0)*(u-1/2))+miu_a/sig_a)
    a0<-rtruncnorm(J,a=0,b=Inf,mean=m_a,sd=sqrt(V_a))
    Ra[i,]<-a0
    M_a0<-matrix(1,nrow=N,ncol=1)%*%a0
    
    ##b, prior N(0,10)
    miu_b<-0
    sig_b<-10
    V_b<-1/(k*colSums(M_a0^2*w0)+1/sig_b)
    m_b<-V_b*(k*colSums(M_a0^2*w0*M_theta0-M_a0*(u-1/2))+miu_b/sig_b)
    b0<-rnorm(J,m_b,sqrt(V_b))
    Rb[i,]<-b0
    M_b0<-matrix(1,nrow=N,ncol=1)%*%b0
  }
  t0<- Sys.time()
  time<-difftime(t0,t1,units="hours")
  return(list(Rtheta=Rtheta[,-(1:nburn)],Ra=Ra[-(1:nburn),],Rb=Rb[-(1:nburn),],time=time))
}


#########Preparation before parallelization
library(parallel)                       #Loading Parallel Package
cores <- detectCores()                  #Calculating the number of cores of a computer
cl <- makeCluster(cores-2,type="PSOCK") #Initializing the number of cores to be used  

#########Preparation before Partitioning
set.seed(1818)
r_num<-sample(c(1:N),N)
k<-2          ##Number of subsets (k=2,5,10,20)
s_num<-N/k    ##Sample size of subsets 
ntime<-10000
nburn<-5000
nrep<-length(y)
N<-nrow(y[[1]])
J<-ncol(y[[1]])

wasp<-as.list(rep(0,nrep))
start1 <- Sys.time()
for(i in 1:nrep){
  ################################LS-WASP algorithm#############################
  
  ####Stage 1: Partitioning of the Full Data
  y_subset<-lapply(c(1:k),function(x) y[[i]][r_num[((x-1)*s_num+1):(x*s_num)],])
  
  clusterExport(cl,c('y_subset','ntime','nburn','k'))
  clusterEvalQ(cl, library(truncnorm))
  clusterEvalQ(cl, library(BayesLogit))
  
  ####Stage 2: Parallel Sampling of Each Subset using Modified Likelihood method
  start <- Sys.time()
  results_group <- parLapply(cl=cl,1:k,PG2PL_group_Par) 
  
  ####Stage 3: Assembling and Integrating Sampled Parameters from Each Subset
  group_temp<-lapply(results_group,function(x){
    Rtheta<-apply(x$Rtheta,1,mean)
    Ra<-t(x$Ra)
    m_a<-apply(Ra,1,mean)
    V_a<-apply((Ra-m_a%*%t(rep(1,(ntime-nburn))))^2,1,mean)
    U_a<-(Ra-m_a%*%t(rep(1,(ntime-nburn))))/sqrt(V_a)
    
    Rb<-t(x$Rb)
    m_b<-apply(Rb,1,mean)
    V_b<-apply((Rb-m_b%*%t(rep(1,(ntime-nburn))))^2,1,mean)
    U_b<-(Rb-m_b%*%t(rep(1,(ntime-nburn))))/sqrt(V_b)
    return(list(Rtheta=Rtheta,m_a=m_a,m_b=m_b,V_a=V_a,U_a=U_a,V_b=V_b,U_b=U_b)) 
  })
  theta_subset<-unlist(lapply(group_temp,function(x) x$Rtheta))
  m_a_bar<-apply(Reduce('cbind',lapply(group_temp,function(x) x$m_a)),1,mean)
  V_a_bar<-(apply(sqrt(Reduce('cbind',lapply(group_temp,function(x) x$V_a))),1,mean))^2
  U_a<-lapply(group_temp,function(x) x$U_a)
  a_subset<-apply(Reduce("cbind",lapply(U_a,function(x) m_a_bar%*%t(rep(1,(ntime-nburn)))+sqrt(V_a_bar)%*%t(rep(1,(ntime-nburn)))*x)),1,mean)
  
  m_b_bar<-apply(Reduce('cbind',lapply(group_temp,function(x) x$m_b)),1,mean)
  V_b_bar<-(apply(sqrt(Reduce('cbind',lapply(group_temp,function(x) x$V_b))),1,mean))^2
  U_b<-lapply(group_temp,function(x) x$U_b)
  b_subset<-apply(Reduce("cbind",lapply(U_b,function(x) m_b_bar%*%t(rep(1,(ntime-nburn)))+sqrt(V_b_bar)%*%t(rep(1,(ntime-nburn)))*x)),1,mean)
  end <- Sys.time()
  
  time_total<-difftime(end,start,units="hours")
  time_subset<-max(unlist(lapply(results_group,function(x) x$time)))
  wasp[[i]]<-list(theta_hat=theta_subset,a_hat=a_subset,b_hat=b_subset,time=time_subset,time_total=time_total)
}
end1 <- Sys.time()
difftime(end1,start1,units="hours")
stopCluster(cl) ###Turn off parallel mode

#####bias, RMSE and time
Bias_wasp<-Reduce("+",lapply(wasp,function(x) c(mean(x$theta_hat-theta[r_num,]),mean(x$a_hat-a[1,]),mean(x$b_hat-b[1,]))))/nrep
RMSE_wasp<-sqrt(Reduce("+",lapply(wasp,function(x) c(mean((x$theta_hat-theta[r_num,])^2),mean((x$a_hat-a[1,])^2),mean((x$b_hat-b[1,])^2))))/nrep)
names(Bias_wasp)<-c("theta_Bias","a_Bias","b_Bias")
names(RMSE_wasp)<-c("theta_RMSE","a_RMSE","b_RMSE")

max_time<-max(unlist(lapply(wasp,function(x) x$time_total)))
max_time


