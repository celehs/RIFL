### console inputs: 
###  S = sample size 500, 1000, 2000
###  M = separation level 1,2,3,4,5
###  P = number of majority sites 6 or 8
args <- commandArgs(TRUE)
print(args)
## Parse the arguments
if(length(args) == 0) {
  print("No arguments supplied.")
} else {
  for (i in 1:length(args)) {
    eval(parse(text = args[i]))
  }
}


library(MASS)
library(igraph)
library(nleqslv)

# import RIFL helper functions and functions to estimate TATE
source("tate.R")
source("helper/helper_fns_uni.R")


######################
### main functions ###
######################
RIFL.TATE = function(source.data.list, target.data, or.formula, ps.formula, gamma.formula,link, seed, median_B=500,M=500){
  L = length(source.data.list)
  source_est = rep(NA,L)
  source_se = rep(NA,L)
  for (l in 1:L){
    source.data = source.data.list[[l]]
    tate.site = tate(source.data,target.data, or.formula, ps.formula, gamma.formula,link)
    source_est[l] = tate.site[1]
    source_se[l] = tate.site[2]
  }
  
  Lhat = L_hat(source_est)
  se_Lhat = se_L_hat(source_se)
  
  for (myrho in seq(0.25,1.2,0.025)){  
    set.seed(seed)
    V.mat = replicate(M,RIFL.one(Lhat = Lhat, se_Lhat = se_Lhat, rho = myrho, alpha = 0.05/20))
    Vhat.mat = V.mat[1:(dim(V.mat)[1]/2),]
    size = apply(Vhat.mat,2,sum)
    percent = sum(size>(dim(Lhat)[2]*thres)) / M
    if (percent >= 0.1) {
      # print(myrho)
      break
    }
  }
  print(myrho)
  print(percent)
  
  set.seed(seed)
  V.mat = replicate(M,RIFL.one(Lhat = Lhat,se_Lhat = se_Lhat, rho = myrho, alpha = 0.05/20))
  Vhat.mat = V.mat[1:(dim(V.mat)[1]/2),]
  cali_M = which(apply(Vhat.mat,2,sum) > (dim(Lhat)[2]*thres))
  Vhat.CI = sapply(cali_M,function(m){IVW(source_est, source_se, which(Vhat.mat[,m]==1), 0.05*19/20)})
  Vtilde.mat = V.mat[-c(1:(dim(V.mat)[1]/2)),]
  Vtilde.CI = sapply(cali_M,function(m){IVW(source_est, source_se, which(Vtilde.mat[,m]==1), 0.05*19/20)})
  
  CI.low.RIFL = min(Vhat.CI[2,])
  CI.up.RIFL = max(Vhat.CI[3,])
  est.RIFL = (CI.low.RIFL + CI.up.RIFL)/2
  res.RIFL = c(est.RIFL,CI.low.RIFL,CI.up.RIFL)
  
  CI.low.RIFL.alt = min(Vtilde.CI[2,])
  CI.up.RIFL.alt = max(Vtilde.CI[3,])
  est.RIFL.alt = (CI.low.RIFL.alt + CI.up.RIFL.alt)/2
  res.RIFL.alt = c(est.RIFL.alt,CI.low.RIFL.alt,CI.up.RIFL.alt)
  
  ### true oracle, Xiudi implemented this 11/19
  res.oracle = IVW(source_est,source_se,true.prevail)
  
  ### MC
  MC = gen_V_m(Lhat, se_Lhat, rho = 1, alpha = 0.05)[1,]
  res.MC = IVW(source_est,source_se,which(MC==1))
  
  ### median estimator
  res.median = median_CI(source_est,source_se,median_B)
  
  return(list(c(res.oracle,res.median,res.MC,res.RIFL,res.RIFL.alt),source_est,source_se,V.mat))
}

do.one.RIFL = function(or.formula,ps.formula,gamma.formula,link,seed){
  set.seed(seed)
  target.data = mvrnorm(N, mu = mean.target, Sigma = Cov)
  target.data = data.frame(target.data)
  colnames(target.data) = sapply(c(1:p), function(x){paste("X",x,sep="")})
  
  source.data.list = vector("list",K)
  target.data.list = vector("list",K)
  if (link == "logit"){
    for (k in 1:K){
      X = mvrnorm(n, mu = means[,k], Sigma = Cov)
      prob0 = exp(X %*% beta.ctrl+ intercepts[k]) / (1 + exp(X %*% beta.ctrl+ intercepts[k]))
      prob1 = exp(trt.effects[k]+ intercepts[k] + X %*% beta.ctrl) / (1 + exp(trt.effects[k] + intercepts[k]+ X %*% beta.ctrl))
      Y0 = rbinom(n,1,prob0)
      Y1 = rbinom(n,1,prob1)
      probA = exp(X %*% betaA + betaA.interaction*X[,1]*X[,2])/(1+exp(X %*% betaA + betaA.interaction*X[,1]*X[,2]))
      A = rbinom(n,1,probA)
      Y = A*Y1 + (1-A)*Y0
      
      source.data = data.frame(cbind(Y,A,X))
      colnames(source.data) = c("Y","A",sapply(c(1:p), function(x){paste("X",x,sep="")}))
      source.data.list[[k]] = source.data
    }
  } else if (link == "identity"){
    for (k in 1:K){
      X = mvrnorm(n, mu = means[,k], Sigma = Cov)
      mean0 = X %*% beta.ctrl + intercepts[k]
      mean1 = trt.effects[k] + X %*% beta.ctrl + intercepts[k]
      Y0 = rnorm(n,mean0,1)
      Y1 = rnorm(n,mean1,1)
      probA = exp(X %*% betaA+ betaA.interaction*X[,1]*X[,2])/(1+exp(X %*% betaA+ betaA.interaction*X[,1]*X[,2]))
      A = rbinom(n,1,probA)
      Y = A*Y1 + (1-A)*Y0
      
      source.data = data.frame(cbind(Y,A,X))
      colnames(source.data) = c("Y","A",sapply(c(1:p), function(x){paste("X",x,sep="")}))
      source.data.list[[k]] = source.data
    }
  }
  RIFL.TATE(source.data.list,target.data,or.formula,ps.formula,gamma.formula,link,seed)
}


###########################
#### simulation set up ####
###########################
A1gen <- function(rho,p){
  A1=matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      A1[i,j]<-rho^(abs(i-j))
    }
  }
  A1
}
rho = 0.6
p=10
mean.target = c(rep(0.5,2),rep(0,p-2))
mean.source = rep(0,p)

means = cbind(mean.target,mean.target,mean.target,
              mean.source,mean.source,mean.source,
              mean.target,mean.source,
              mean.target,mean.source)

Cov <- A1gen(rho,p)
betaA = c(0.5,-0.5,rep(0,p-2))
betaA.interaction = 0.1
beta.ctrl = c(rep(0.5,5),rep(0.1,3),rep(0,2))
intercepts = c(0.05,-0.05,0.1,-0.1,0.05,-0.05,0.1,-0.1,0,0)
n = S
N = 10000
a = M
K = 10
Lmaj = P
thres = 0.5  # thres = 0.5 corresponds to the majority rule
#  change to 0.7 to get RIFL with 80% rule

if (Lmaj == 6){
  trt.effects = c(rep(-1,6),rep(-1-0.2*a,2),rep(-1-0.1*a,2))
  true.prevail = c(1:6)
}

if (Lmaj == 8){
  trt.effects = c(rep(-1,8),rep(-1-0.2*a,1),rep(-1-0.1*a,1)) 
  true.prevail = c(1:8)
}

nsim = 200
allres = matrix(0,nsim,15)
site.est = matrix(0,nsim,10)
site.se = matrix(0,nsim,10)
V.mat.list = vector("list",nsim)
for (i in 1:nsim){
  print(i)
  res.one = do.one.RIFL("Y~.","A~X1+X2","Y~.","identity",i*3)
  allres[i,] = res.one[[1]]
  site.est[i,] = res.one[[2]]
  site.se[i,] = res.one[[3]]
  V.mat.list[[i]] = res.one[[4]]
}

filename = paste("n",n,"diff",a,"TATEmaj",P,".RData",sep="") 
save(allres,site.est,site.se,V.mat.list,file=filename)


###########################
### process the results ###
###########################
res_summary = NULL
for (n in c(500,1000,2000)){
  for (m in c(1:5)){
    for (s in c(4)){
      filename = paste("n",n,"diff",m,"nu",s,"tateMNBmaj6.RData",sep="")
      try({load(filename)
        cov = mean(allres[,2] < (-1) & allres[,3]>(-1))
        length = mean(allres[,3] - allres[,2])
        res_summary = rbind(res_summary,c(n,m,s,cov,length))})
    }
  }
}

res_summary = data.frame(res_summary)
colnames(res_summary) = c("n","sep","nu","cov","length")
res_summary$nu = (res_summary$nu-1)*0.1+0.5

res_boot = res_summary
res_mnb = data.frame(est = rep(0,dim(res_boot)[1]),
                     cov = res_boot$cov,
                     width=res_boot$length,
                     n = res_boot$n,
                     separation = res_boot$sep,
                     method = rep("MNB",dim(res_boot)[1]))

result_summary = NULL
for (n in c(500,1000,2000)){
  for (m in c(1:5)){
    filename = paste("n",n,"diff",m,"TATEmaj6.RData",sep="")
    load(filename)
    est = c()
    cov = c()
    width=c()
    for (i in 1:5){
      est = c(est,mean(allres[,(3*i-2)]))
      cov = c(cov,mean(allres[,(3*i-1)] < (-1) & allres[,(3*i)]>(-1)))
      width = c(width,mean(allres[,(3*i)] - allres[,(3*i-1)] ))
      
    }
    points = allres[,7]
    bias = points- (-1)
    se = (allres[,9] - allres[,8])/2/qnorm(0.975)
    sd = sqrt(var(bias))
    se.new = sd/mean(se)*se
    temp = mean(bias)/se.new
    cv_alpha = sqrt(qchisq(1-0.05,df=1, ncp=temp^2))
    CI.new = cbind(points - se.new*cv_alpha, points + se.new*cv_alpha)
    cov.new = mean((CI.new[,1] < (-1))*(CI.new[,2] > (-1)))
    len.new = mean(2*se.new*cv_alpha)
    
    est = c(est,est[3])
    cov = c(cov,cov.new)
    width = c(width,len.new)
    res_summary = data.frame(est=est,
                             cov=cov,
                             width=width,
                             n=rep(n,6),
                             separation = rep(m,6),
                             method=c("oracle","median","VMC","RIFL.hat","RIFL","OBA"))
    result_summary = rbind(result_summary,res_summary)
  }
}

result_summary = rbind(result_summary,res_mnb)


res_summary = NULL
for (n in c(500,1000,2000)){
  for (m in c(1:5)){
    for (s in c(4)){
      filename = paste("n",n,"diff",m,"nu",s,"tateMNBmaj8.RData",sep="")
      try({load(filename)
        cov = mean(allres[,2] < (-1) & allres[,3]>(-1))
        length = mean(allres[,3] - allres[,2])
        res_summary = rbind(res_summary,c(n,m,s,cov,length))})
    }
  }
}

res_summary = data.frame(res_summary)
colnames(res_summary) = c("n","sep","nu","cov","length")
res_summary$nu = (res_summary$nu-1)*0.1+0.5

res_boot = res_summary
res_mnb = data.frame(est = rep(0,dim(res_boot)[1]),
                     cov = res_boot$cov,
                     width=res_boot$length,
                     n = res_boot$n,
                     separation = res_boot$sep,
                     method = rep("MNB",dim(res_boot)[1]))

result_summary = NULL
for (n in c(500,1000,2000)){
  for (m in c(1:5)){
    filename = paste("n",n,"diff",m,"TATEmaj8.RData",sep="")
    load(filename)
    est = c()
    cov = c()
    width=c()
    for (i in 1:5){
      est = c(est,mean(allres[,(3*i-2)]))
      cov = c(cov,mean(allres[,(3*i-1)] < (-1) & allres[,(3*i)]>(-1)))
      width = c(width,mean(allres[,(3*i)] - allres[,(3*i-1)] ))
      
    }
    points = allres[,7]
    bias = points- (-1)
    se = (allres[,9] - allres[,8])/2/qnorm(0.975)
    sd = sqrt(var(bias))
    se.new = sd/mean(se)*se
    temp = mean(bias)/se.new
    cv_alpha = sqrt(qchisq(1-0.05,df=1, ncp=temp^2))
    CI.new = cbind(points - se.new*cv_alpha, points + se.new*cv_alpha)
    cov.new = mean((CI.new[,1] < (-1))*(CI.new[,2] > (-1)))
    len.new = mean(2*se.new*cv_alpha)
    
    est = c(est,est[3])
    cov = c(cov,cov.new)
    width = c(width,len.new)
    res_summary = data.frame(est=est,
                             cov=cov,
                             width=width,
                             n=rep(n,6),
                             separation = rep(m,6),
                             method=c("oracle","median","VMC","RIFL.hat","RIFL","OBA"))
    result_summary = rbind(result_summary,res_summary)
  }
}



result_summary = rbind(result_summary,res_mnb)


result_summary = result_summary %>% filter(method != "RIFL.hat")
get_legend = function(myggplot){
  tmp = ggplot_gtable(ggplot_build(myggplot))
  leg = which(sapply(tmp$grobs,function(x) x$name) == "guide-box")
  legend = tmp$grobs[[leg]]
  return(legend)
}

library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggpubr)
datplot = result_summary %>% filter(n == 500 & method!="oracle" & method!="RIFL.hat")
p1 = ggplot(datplot,aes(x=separation,y=cov, group=as.factor(method), col=as.factor(method)))+
  #ggtitle("n=500, coverage")+
  xlab("separation")+ylab("coverage")+geom_point()+geom_line()+ 
  geom_hline(yintercept=0.95, linetype="dashed")+
  theme(legend.direction="horizontal",legend.title=element_blank())
mylegend = get_legend(p1)
p1 = p1+theme(legend.position="none")
p2 = ggplot(datplot,aes(x=separation,y=width, group=as.factor(method), col=as.factor(method)))+
  #ggtitle("n=500, length")+
  xlab("separation")+ylab("length")+geom_point()+geom_line()+theme(legend.position="none")+
  geom_hline(yintercept=result_summary$width[result_summary$method=="oracle" & result_summary$n == 500][1], linetype="dashed")
p12 = ggarrange(p1,p2, ncol=2, nrow=1)
p12 = annotate_figure(p12, top = text_grob("n=500", 
                                           color = "black", face = "bold", size = 10))
#grid.arrange(arrangeGrob(p1,p2,ncol=2,nrow=1),mylegend,nrow=2,heights=c(9,1))


datplot = result_summary %>% filter(n == 1000 & method!="oracle" & method!="RIFL.hat")
p3 = ggplot(datplot,aes(x=separation,y=cov, group=as.factor(method), col=as.factor(method)))+
  #ggtitle("n=1000, coverage")+
  xlab("separation")+ylab("coverage")+geom_point()+geom_line()+ 
  geom_hline(yintercept=0.95, linetype="dashed")+
  theme(legend.direction="horizontal",legend.title=element_blank())
#mylegend = get_legend(p1)
p3 = p3+theme(legend.position="none")
p4 = ggplot(datplot,aes(x=separation,y=width, group=as.factor(method), col=as.factor(method)))+
  #ggtitle("n=1000, length")+
  xlab("separation")+ylab("length")+geom_point()+geom_line()+theme(legend.position="none")+
  geom_hline(yintercept=result_summary$width[result_summary$method=="oracle" & result_summary$n == 1000][1], linetype="dashed")
p34 = ggarrange(p3,p4, ncol=2, nrow=1)
p34 = annotate_figure(p34, top = text_grob("n=1000", 
                                           color = "black", face = "bold", size = 10))

#grid.arrange(arrangeGrob(p1,p2,ncol=2,nrow=1),mylegend,nrow=2,heights=c(9,1))

datplot = result_summary %>% filter(n == 2000 & method!="oracle" & method!="RIFL.hat")
p5 = ggplot(datplot,aes(x=separation,y=cov, group=as.factor(method), col=as.factor(method)))+
  #ggtitle("n=2000, coverage")+
  xlab("separation")+ylab("coverage")+geom_point()+geom_line()+ 
  geom_hline(yintercept=0.95, linetype="dashed")+
  theme(legend.direction="horizontal",legend.title=element_blank())
#mylegend = get_legend(p1)
p5 = p5+theme(legend.position="none")
p6 = ggplot(datplot,aes(x=separation,y=width, group=as.factor(method), col=as.factor(method)))+
  #ggtitle("n=2000, length")+
  xlab("separation")+ylab("length")+geom_point()+geom_line()+theme(legend.position="none")+
  geom_hline(yintercept=result_summary$width[result_summary$method=="oracle" & result_summary$n == 2000][1], linetype="dashed")

p56 = ggarrange(p5,p6, ncol=2, nrow=1)
p56 = annotate_figure(p56, top = text_grob("n=2000", 
                                           color = "black", face = "bold", size = 10))

#grid.arrange(arrangeGrob(p1,p2,p3,p4,p5,p6,ncol=2,nrow=3),mylegend,nrow=2,heights=c(9,1))
grid.arrange(arrangeGrob(p12,p34,p56,ncol=1,nrow=3),mylegend,nrow=2,heights=c(9,1))
