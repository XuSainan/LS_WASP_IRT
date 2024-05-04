######################################################################################
######################Data and Parameter Generation(M2PL model)#######################
######################################################################################
library(MASS)
set.seed(1818)
N<-10000             ##Number of examinees
J<-20                ##Number of items
r<-2                 ##Number of dimensions of latent traits
Sigma<-diag(1,r)
theta<-mvrnorm(N,rep(0,r),Sigma)
a<-matrix(rlnorm(J*r,0.3,0.2),nrow=J,ncol=r)
a[1:r,]<-diag(1,r)
b<-matrix(rnorm(J),nrow=1)
b[1:r]<-0
M_b<-matrix(1,nrow=N,ncol=1)%*%b
p<-1/(1+exp(M_b-theta%*%t(a)))
nrep<-25             ##Number of replications
y<-as.list(rep(0,nrep))
set.seed(8888)
for(i in 1:nrep){
  u<-matrix(runif(N*J),nrow=N,ncol=J)
  y[[i]]<-ifelse(p>u,1,0)
}



#################################################################################################
###############################LS-WASP algorithm for M2PL model##################################
#################################################################################################

###Polya-Gamma algorithm (parallel version for each subset) using Modified Likelihood method
PGM2PL_group_Par<-function(x){
  t1<-Sys.time()
  u<-y_subset[[x]]
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
    
    ##b, prior N(0,10)
    miu_b<-0
    sig_b<-10
    V_b<-1/(k*colSums(w0)+1/sig_b)
    m_b<-V_b*(k*colSums((theta0%*%t(a0))*w0-(u-1/2))+miu_b/sig_b)
    b0<-rnorm(J,m_b,sqrt(V_b))
    b0[1:r]<-0
    Rb[i,]<-b0
    M_b0<-matrix(1,nrow=N,ncol=1)%*%b0
    
    ##a, prior N(0,10)
    miu_a<-0
    sig_a<-10
    V_a<-1/(k*(t(w0)%*%theta0^2)+1/sig_a)
    a_temp<-lapply(c(1:r),function(d) k*(t(theta0[,d])%*%((u-1/2)+w0*M_b0-w0*(theta0[,-d]%*%t(a0[,-d]))))+miu_a/sig_a)
    m_a<-V_a*t(Reduce("rbind",a_temp))
    a0<-matrix(rtruncnorm(1,a=0,b=Inf,mean=t(m_a),sd=t(sqrt(V_a))),J,r,byrow = T)
    a0[1:r,]<-diag(1,r)
    Ra[[i]]<-a0
  }
  t0<- Sys.time()
  time<-difftime(t0,t1,units="hours")
  return(list(Rtheta=Rtheta[-(1:nburn)],Ra=Ra[-(1:nburn)],Rb=Rb[-(1:nburn),],time=time))
}


#########Preparation before parallelization
library(parallel)                       #Loading Parallel Package
cores <- detectCores()                  #Calculating the number of cores of a computer
cl <- makeCluster(cores-2,type="PSOCK") #Initializing the number of cores to be used


#########Preparation before Partitioning
set.seed(1818)
r_num<-sample(c(1:N),N)
k<-5                ##Number of subsets (k=2,5,10)
s_num<-N/k          ##Sample size of subsets
ntime<-10000
nburn<-5000
nrep<-length(y)
N<-nrow(y[[1]])
J<-ncol(y[[1]])
r<-ncol(theta)

wasp<-as.list(rep(0,nrep))
start1 <- Sys.time()
for(i in 1:nrep){
  ################################LS-WASP algorithm#############################
  
  ####Stage 1: Partitioning of the Full Data
  y_subset<-lapply(c(1:k),function(x) y[[i]][r_num[((x-1)*s_num+1):(x*s_num)],])
  
  clusterExport(cl,c('y_subset','ntime','nburn','k','r'))
  clusterEvalQ(cl, library(truncnorm))
  clusterEvalQ(cl, library(BayesLogit))
  clusterEvalQ(cl, library(MASS))
  
  ####Stage 2: Parallel Sampling of Each Subset using Modified Likelihood method
  start <- Sys.time()
  results_group <- parLapply(cl=cl,1:k,PGM2PL_group_Par) 
  
  ####Stage 3: Assembling and Integrating Sampled Parameters from Each Subset
  group_temp<-lapply(results_group,function(x){
    Rtheta<-Reduce("+",x$Rtheta)/(ntime-nburn)
    Ra1<-x$Ra
    m_a<-Reduce("+",Ra1)/(ntime-nburn)
    V_a<-Reduce("+",lapply(Ra1,function(x) (x-m_a)^2))/(ntime-nburn)
    U_a<-lapply(Ra1,function(x) rbind(diag(1,r),((x[-(1:r),]-m_a[-(1:r),])/sqrt(V_a[-(1:r),]))))
    
    Rb<-t(x$Rb)
    m_b<-apply(Rb,1,mean)
    V_b<-apply((Rb-m_b%*%t(rep(1,(ntime-nburn))))^2,1,mean)
    U_b<-rbind(matrix(0,r,(ntime-nburn)),((Rb[-(1:r),]-m_b[-(1:r)]%*%t(rep(1,(ntime-nburn))))/sqrt(V_b[-(1:r)])))
    return(list(Rtheta=Rtheta,m_a=m_a,m_b=m_b,V_a=V_a,U_a=U_a,V_b=V_b,U_b=U_b)) 
  })
  theta_subset<-Reduce(rbind,lapply(group_temp,function(x) x$Rtheta))
  m_a_bar<-Reduce("+",lapply(group_temp,function(x) x$m_a))/k
  V_a_bar<-(Reduce("+",lapply(group_temp,function(x) sqrt(x$V_a)))/k)^2
  U_a_temp<-lapply(group_temp,function(x) x$U_a)
  U_a<-NULL
  for(j in 1:k){
    U_a<-c(U_a,U_a_temp[[j]])
  }
  a_draw<-lapply(U_a,function(x) m_a_bar+sqrt(V_a_bar)*x)
  a_subset<-Reduce("+",a_draw)/(k*(ntime-nburn))
  
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
Bias_wasp<-Reduce("+",lapply(wasp,function(x) c(mean(x$theta_hat-theta[r_num,]),mean(x$a_hat-a),mean(x$b_hat-b[1,]))))/nrep
RMSE_wasp<-sqrt(Reduce("+",lapply(wasp,function(x) c(mean((x$theta_hat-theta[r_num,])^2),mean((x$a_hat-a)^2),mean((x$b_hat-b[1,])^2))))/nrep)
names(Bias_wasp)<-c("theta_Bias","a_Bias","b_Bias")
names(RMSE_wasp)<-c("theta_RMSE","a_RMSE","b_RMSE")

max_time<-max(unlist(lapply(wasp,function(x) x$time_total)))
max_time
