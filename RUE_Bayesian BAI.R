library(nnet)
######################################
#rm(list = ls())
######################################
##### updating variance estimates ######
update_variance<-function(mean_r,mean_r2,size,M_r){
  ## element in mean_r: reward mean of each arm
  ## element in mean_r2: mean of reward^2 of each arm
  ## element in size: size of each arm
    N<-max(sum(size),1);
    N1<-max(sum(size[size>=1]-1),1);#max(sum(size[size>=1]-1)-1,1); ## remove the arms of size not more than 1. 
    sigmar2<-sum(size*mean_r2-size*mean_r^2)/N1;
    ######
    u<-mean_r-M_r;
    Nstar<-max(N-sum(size^2)/N,1);
    sigmamu2<-sum(size*u^2)/Nstar;
    if(N1==0){sigmar2<-0}
    if(N==0){sigmamu2<-0;sigmar2<-0}
  ###
  output=list("sigmamu2"=sigmamu2,"sigmar2"=sigmar2)
  return(output)
}


#### the estimation of the random-effect model ####
blup<-function(sigmamu2,sigmar2,mean_r,size,M_r){
  #blup estimator
  K<-length(size);
    if((sigmamu2==0)&&(sigmar2==0)){
      w<-rep(0,K);
    } else {
      w<-size*sigmamu2/(size*sigmamu2+sigmar2);
    }
    w[size==0]<-0;
    AA<-(1-w)*size;
    sum_1ws<-max(sum((1-w)*size),1);
    bar<-sum(mean_r*AA)/sum_1ws;
    muhat<-w*mean_r+(1-w)*bar;
    tau2<-sigmar2*(w/size+(1-w)^2/sum_1ws);
  ###
  output=list("muhat"=muhat,"tau2"=tau2)
  return(output)
}


TS<-function(timek,alpha,beta){
  thetahat<-c();
  K<-length(timek)
  for(k in 1:K){thetahat[k]<-rbeta(1,shape1=alpha[k],shape2=beta[k]);}
  K2<-c(1:K)[thetahat==max(thetahat)]
  if(length(K2)>1){
    xidx<-sample(K2,1);
  } else {xidx<-which.max(thetahat)}
}

TTTS<-function(timek,alpha,beta){
  BernChoice<-rbinom(1,1,prob=0.5)
  I<-TS(timek,alpha,beta);
  if(BernChoice==1){
    xidx<-I
  } else {
    J<-I;
    kk<-0;
    while((J==I)&(kk<=length(timek))){
      J<-TS(timek,alpha,beta);
      kk<-kk+1;
    }
    xidx<-J
  }
}


######### main code #######
BAI_RE<-function(theta,B=50,method="BLUP",N=10^5,DEG=2,MODEL="Binary",recommend="mean"){
  K<-length(theta);
  H1<-ceiling(sum(1/(max(theta)-theta[-which.max(theta)])^2));
  RightNum<-0;
  ##
  for(b in 1:B){
    alpha<-rep(1,K);beta<-rep(1,K);
    timek<-rep(0,K);mean_r<-rep(0,K);mean_r2<-rep(0,K);
    sigmamu2<-0;sigmar2<-0;M_r<-0;
    for(t in 1:N){
      thetahat<-c();
      if(method=="BLUP"){ ### BLUP estimates (random-effect model)
      size<-timek;
      ### estimate the variances 
      update.variance<-update_variance(mean_r,mean_r2,size,M_r);
      sigmamu2<-update.variance$sigmamu2;
      sigmar2<-update.variance$sigmar2;
      ### update 
      estimation<-blup(sigmamu2,sigmar2,mean_r,size,M_r);
      muhat<-estimation$muhat;
      tau2<-estimation$tau2;
      bound<-(tau2)^0.5;
      }
      if(method=="TTTS"){
        if(t<=(K)){xidx<-t
        } else if((t<=2*K)&(t>K)){xidx<-t-K;
        } else {
          xidx<-TTTS(timek,alpha,beta);}
      } else {
       for(k in 1:K){
       #################### begin 
        if(method=="BLUP"){
          if(timek[k]<=1){
            thetahat[k]<-1000;#n^0.5
          } else {
            thetahat[k]<-muhat[k]+bound[k]*(DEG*log(N))^0.5;
          }
        ####
        } else if(method=="UCBE"){
          if(timek[k]<1){
            thetahat[k]<-1000;
          } else {
            ucb_a<-2*N/H1;
            thetahat[k]<-mean_r[k]+(ucb_a/timek[k])^0.5;
          }
        }
       }
      ################## end
       K2<-c(1:K)[thetahat==max(thetahat)]
       if(length(K2)>1){
        xidx<-sample(K2,1);
       } else {xidx<-which.max(thetahat)}
      ## end 
      }
      ################
      if(MODEL=="Bernoulli"){
        r<-rbinom(n=1,1,prob=theta[xidx]);
      } else if (MODEL=="Unif") {
        rx<-runif(1,0,1);
        thetax<-theta[xidx];
        r<-sum(rx<thetax)
      } else if (MODEL=="Gaussian"){
        r<-rnorm(n=1,theta[xidx],sd=1);
      } else if (MODEL=="Beta4"){
        thetax<-rbeta(n=1,shape1=4*theta[xidx],shape2=4*(1-theta[xidx]));
        r<-thetax;
      }
      alpha[xidx]<-alpha[xidx]+r;
      beta[xidx]<-beta[xidx]+1-r;
      timek[xidx]<-timek[xidx]+1;
      M_r<-M_r+1/t*(r-M_r);
      mean_r[xidx]<-mean_r[xidx]+1/timek[xidx]*(r-mean_r[xidx]);
      mean_r2[xidx]<-mean_r2[xidx]+1/timek[xidx]*(r^2-mean_r2[xidx]);
      ## recommendation step
        if(t==N){
          if(recommend=="times"){
            istar<-which.max(timek);
          } else if (recommend=="mean"){
            istar<-which.max(mean_r);
          } else if (recommend=="postmean"){
            size<-timek;
            update.variance<-update_variance(mean_r,mean_r2,size,M_r);
            sigmamu2<-update.variance$sigmamu2;
            sigmar2<-update.variance$sigmar2;
            estimation<-blup(sigmamu2,sigmar2,mean_r,size,M_r);
            muhat<-estimation$muhat;
            istar<-which.max(muhat);
          }
          if(istar==1){RightNum<-RightNum+1;}
      }
    ##
    }
  }
  ###
  return(RightNum/B)
}



####################################################################################
### phase-based uniform exploration methods, including SR, SH, and Two-stage  ######
####################################################################################
halving<-function(St,pt){
  ## St the index of arms at round t, and pt the corresponding values
  nt<-length(St);
  pt_sort<-sort.int(pt,decreasing=T,index.return=T)
  St<-St[pt_sort$ix[1:ceiling(nt/2)]];
  results<-list("St"=St,"pt"=pt)
  return(results)
}


# Sequential Halving Algorithm
SH = function(n, p, TT,MODEL="Bernoulli"){
  r_max = floor(log(n, 2))
  S = 1:n
  psum_hat<-rep(0,n)
  nsum<-rep(0,n)
  while(length(S)>1){
    times = max(floor(TT/(length(S)*ceiling(log(n, 2)))),1);#
    for(i in S){
      if(MODEL=="Bernoulli"){
        psum_hat[i]<-psum_hat[i]+sum(rbinom(times, 1, p[i]))
      } else if(MODEL=="Beta4"){
        psum_hat[i]<-psum_hat[i]+sum(rbeta(times,shape1=4*p[i],shape2=4*(1-p[i])))
      } else if(MODEL=="Gaussian"){
        psum_hat[i]<-psum_hat[i]+sum(rnorm(times,p[i],sd=1))
      }
      nsum[i]<-nsum[i]+times;
    }
    p_hat<-psum_hat/nsum
    HR<-halving(S,p_hat[S]);
    S<-HR$St;
  }
  return(S)
}


# Successive Rejects Algorithm
delete.one<-function(St,pt){
  ## St the index of arms at round t, and pt the corresponding values
  St<-St[-which.is.max(-pt)];
  return(St)
}

###
SR = function(n, p, TT,MODEL="Bernoulli"){
  log_n<-0.5+sum(1/(2:n));
  S = 1:n
  psum_hat<-rep(0,n)
  nsum<-rep(0,n)
  k<-1;
  while(length(S)>1){
    nk<-1/log_n*(TT-n)/(n+1-k);
    nk_1<-1/log_n*(TT-n)/(n+2-k)
    if(k==1){
      times = max(ceiling(nk),1);
    } else {times = max(ceiling(nk)-ceiling(nk_1),1);}
    ###
    for(i in S){
      if(MODEL=="Bernoulli"){
        psum_hat[i]<-psum_hat[i]+sum(rbinom(times, 1, p[i]))
      } else if(MODEL=="Beta4"){
        psum_hat[i]<-psum_hat[i]+sum(rbeta(times,shape1=4*p[i],shape2=4*(1-p[i])))
      } else if(MODEL=="Gaussian"){
        psum_hat[i]<-psum_hat[i]+sum(rnorm(times,p[i],sd=1))
      }
      nsum[i]<-nsum[i]+times;
    }
    p_hat<-psum_hat/nsum
    S<-delete.one(S,p_hat[S]);
    k<-k+1
  }
  return(S)
}


#############################
# The Two-Stage algorithm
#############################
J_star<-function(St,pt){
  ## St the index of arms at round t, and pt the corresponding values
  B_conf<-(n*log(TT)/(q*TT))^0.5;
  L<-pt-B_conf;U<-pt+B_conf;
  St<-St[U>=max(L)];
  return(St)
}

###
TwoStage = function(n, p, q=0.5, TT,MODEL="Bernoulli"){
  S = 1:n
  psum_hat<-rep(0,n)
  nsum<-rep(0,n)
  k<-1;
  while((length(S)>1)&(k==1)){
    if((k==1)|(q==1)){times<-floor(q*TT/length(S));
    } else {times<-floor((1-q)*TT/length(S));}
    ###
    for(i in S){
      if(MODEL=="Bernoulli"){
        psum_hat[i]<-psum_hat[i]+sum(rbinom(times, 1, p[i]))
      } else if(MODEL=="Beta4"){
        psum_hat[i]<-psum_hat[i]+sum(rbeta(times,shape1=4*p[i],shape2=4*(1-p[i])))
      } else if(MODEL=="Gaussian"){
        psum_hat[i]<-psum_hat[i]+sum(rnorm(times,p[i],sd=1))
      } 
      nsum[i]<-nsum[i]+times;
    }
    p_hat<-psum_hat/nsum
    S<-J_star(S,p_hat[S]);
    k<-k+1
  }
  if(length(S)==1){
    Best<-S;
  } else {Best<-S[which.max(p_hat[S])];}
  return(Best)
}


#########################################################
###### run the SH, SR, or TwoStage algorithms
#########################################################
accuracy_phase = function(n, p, q=0.5, TT, B, MODEL="Bernoulli",method="SH"){
  accura <- rep(0,4)
  ReportIdx<-round(c(TT/8,TT/4,TT/2,TT))
  for(j in 1:4){
    TTj<-ReportIdx[j];
    for(k in 1:B){
      if(method=="SH"){
        i_star = SH(n, p, TTj, MODEL)
      } else if(method=="SR"){i_star = SR(n, p, TTj, MODEL)
      } else if(method=="TwoStage"){i_star = TwoStage(n, p, q,TTj, MODEL)}
      if(1==i_star){accura[j] = accura[j] + 1}  
    }
  }
  return(accura/B)
}


#### calculate H ####
HSize<-function(p){
  H<-ceiling(sum(1/(max(p)-p[-which.max(p)])^2))
  return(H)
}

#####################################################################
##### one example for comparing our method (RUE) to previous ####
#####################################################################
B<-500;
######
L<-4;
R<-c(1/8,1/4,1/2,1);
#ns = c(length(p4),length(p5),length(p6))
ns<-c(40,80)
SR_accuracy = rep(0, length(ns)*L)
SH_accuracy = rep(0, length(ns)*L)
TTTS_accuracy<-rep(0, length(ns)*L); 
REU_accuracy<-rep(0, length(ns)*L); ## it's UCBE
RUE_accuracy2<-rep(0, length(ns)*L); ## our approach
for(j in 1:length(ns)){
  n = ns[j]
  idx<-L*(j-1)+c(1:L);
  ######
  #####
  p<-rep(0.5,n);p[2]<-0.5-1/(5*n);p[n]<-0.25
  ma<-(p[n]-p[2])/(n-2);
  for(j in 3:n){
    p[j]<-p[j-1]+ma;
  }
  TT0 = HSize(p);
  TT0<-min(TT0,200000);TT0<-max(TT0,n^2+1);TT<-2*TT0;
  SHBAI = accuracy_phase(n, p, q=0.5, TT, B, MODEL="Bernoulli",method="SR");
  SR_accuracy[idx] <- SHBAI;
  SHBAI = accuracy_phase(n, p, q=0.5, TT, B, MODEL="Bernoulli",method="SH");
  SH_accuracy[idx] <- SHBAI;
  for(i in 1:L){
    TTi<-round(TT*R[i]);
    ReBAI = BAI_RE(theta=p,method="TTTS",N=TTi,B=B,DEG=1,MODEL="Bernoulli",recommend="mean");
    TTTS_accuracy[idx[i]] <- ReBAI;
    ReBAIN = BAI_RE(theta=p,method="UCBE",N=TTi,B=B,DEG=2,MODEL="Bernoulli",recommend="mean");
    REU_accuracy2[idx[i]] <- ReBAIN;
    ReBAI = BAI_RE(theta=p,method="BLUP",N=TTi,B=B,DEG=2,MODEL="Bernoulli",recommend="postmean");
    RUE_accuracy2[idx[i]] <- ReBAI;
  }
  print(cbind(SR_accuracy,SH_accuracy,TTTS_accuracy,REU_accuracy2,RUE_accuracy2))
}
results<-cbind(SR_accuracy,SH_accuracy,TTTS_accuracy,REU_accuracy2,RUE_accuracy2)
print(results)

