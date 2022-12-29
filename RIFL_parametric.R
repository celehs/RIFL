### console inputs:
###  N sample size 500, 1000 or 2000
###  M separation level 1,2,3,4,5
###  P number of majority sites 6 or 8

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

# import RIFL helper functions
source('helper/helper_fns.R')

# simulation set up
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
L = 10
a = 0.1*M
nsim = 500
L_maj = P

if (L_maj == 6){
  beta1 = c(rep(0.5,5),rep(0.1,3),rep(0,2))
  beta2 = beta3 = beta4 = beta5 = beta6 = beta1
  beta7 = c(rep(0.5 - 3*a,5),rep(0.1,3),rep(0,2))
  beta8 = c(rep(0.5 - 2*a,5),rep(0.1,3),rep(0,2))
  beta9 = c(rep(0.5 - a,5),rep(0.1,3),rep(0,2))
  beta10 = c(rep(0.5 + a,5),rep(0.1,3),rep(0,2))
  true.prevail = c(1:6)
}

if (L_maj == 8){
  beta1 = c(rep(0.5,5),rep(0.1,3),rep(0,2))
  beta2 = beta3 = beta4 = beta5 = beta6 = beta1
  beta7 = c(rep(0.5 - 3*a,5),rep(0.1,3),rep(0,2))
  beta8 = beta1
  beta9 = c(rep(0.5 - a,5),rep(0.1,3),rep(0,2))
  beta10 = beta1
  true.prevail = c(1:6,8,10)
}


betas = cbind(beta1,beta2,beta3,beta4,beta5,beta6,beta7,beta8,beta9,beta10)
intercepts = c(0.05,-0.05,0.1,-0.1,0.05,-0.05,0.1,-0.1,0,0)
Cov <- A1gen(rho,dim(betas)[1])
n_list = rep(N,10)
thres = 0.5

### function to generate data for one simulation
gen_dat = function(seed, n_list,intercepts,betas,Cov){
  set.seed(seed)
  p = dim(betas)[1]
  L = dim(betas)[2]
  
  theta_list = vector("list",L)
  V_list = vector("list",L)
  source_est = rep(NA,L)
  source_se = rep(NA,L)
  for (l in 1:L){
    X = mvrnorm(n=n_list[l],mu=rep(0,p),Sigma=Cov)
    exp_val = intercepts[l] + X %*% betas[,l]
    prob = exp(exp_val)/(1+exp(exp_val)) 
    y = rbinom(n_list[l], 1, prob)
    
    mod = glm(y~X,family=binomial())
    theta_list[[l]] = mod$coefficients[-1]
    V_list[[l]] = vcov(mod)[-1,-1]
    source_est[l] = mod$coefficients[2]
    source_se[l] = sqrt(V_list[[l]][1,1])
  }
  return(list(theta_list,V_list,source_est,source_se))
}

### take as input site-specific: source_est, source_se, theta_list, V_list
do.one = function(seed, source_est, source_se, theta_list, V_list, n_list, median_B=500,M=1000){
  Lhat = L_hat(source_est)
  se_Lhat = se_L_hat(source_se)
  Dhat = D_hat(theta_list)
  se_Dhat = se_D_hat(theta_list,V_list,n_list)
  
  res.oracle = IVW(source_est,source_se,true.prevail)
  
  ### median estimator
  res.median = median_CI(source_est,source_se,median_B)
  
  ### maximum clique estimator
  MC = gen_V_m(Dhat, se_Dhat, Lhat, se_Lhat, rho = 1, alpha = 0.05, univariate = FALSE)[1,]
  res.MC = IVW(source_est,source_se,which(MC==1))
  
    ### RIFL
    for (myrho in seq(0.4,1.2,0.025)){
      set.seed(seed)
      V.mat = replicate(M,RIFL.one(Dhat=Dhat, se_Dhat=se_Dhat, Lhat = Lhat, 
                                   se_Lhat = se_Lhat, rho = myrho, alpha = 0.05/20, univariate = FALSE))
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
  V.mat = replicate(M,RIFL.one(Dhat=Dhat, se_Dhat=se_Dhat, Lhat = Lhat, 
                                 se_Lhat = se_Lhat, rho = myrho, alpha = 0.05/20, univariate = FALSE))
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
  
  return(list(c(res.oracle,res.median,res.MC,res.RIFL,res.RIFL.alt),V.mat))
}

one.sim.parametric = function(seed,n_list,intercepts,betas,Cov,median_B = 500,M=1000){
  estimates = gen_dat(seed, n_list,intercepts,betas,Cov)
  source_est = estimates[[3]]
  source_se = estimates[[4]]
  theta_list = estimates[[1]]
  V_list = estimates[[2]]
  do.one(seed, source_est, source_se, theta_list, V_list, n_list, median_B, M)
}


allres = matrix(NA,nsim,15)
allVmat = vector("list",nsim)
for (i.sim in 1:nsim){
  print(i.sim)
  one.res = one.sim.parametric(i.sim*9,n_list,intercepts,betas,Cov)
  allres[i.sim,] = one.res[[1]]
  allVmat[[i.sim]] = one.res[[2]]
}
filename = paste("n",N,"diff",10*a,"maj",L_maj,"param.RData",sep="")
save(allres,allVmat,file=filename)

##########################
### process the result ###
##########################
res_summary = NULL
for (n in c(500,1000,2000)){
  for (m in c(1:5)){
    for (s in c(4)){
      filename = paste("n",n,"diff",m,"nu",s,"paramMNBmaj8.RData",sep="")
      try({load(filename)
        cov = mean(allres[,2] < 0.5 & allres[,3]>0.5)
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
    filename = paste("n",n,"diff",m,"maj8param.RData",sep="")
    load(filename)
    est = c()
    cov = c()
    width=c()
    for (i in 1:5){
      est = c(est,mean(allres[,(3*i-2)]))
      cov = c(cov,mean(allres[,(3*i-1)] < 0.5 & allres[,(3*i)]>0.5))
      width = c(width,mean(allres[,(3*i)] - allres[,(3*i-1)] ))
    }
    points = allres[,7]
    bias = points- 0.5
    se = (allres[,9] - allres[,8])/2/qnorm(0.975)
    sd = sqrt(var(bias))
    se.new = sd/mean(se)*se
    temp = mean(bias)/se.new
    cv_alpha = sqrt(qchisq(1-0.05,df=1, ncp=temp^2))
    CI.new = cbind(points - se.new*cv_alpha, points + se.new*cv_alpha)
    cov.new = mean((CI.new[,1] < 0.5)*(CI.new[,2] > 0.5))
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
#grid.arrange(arrangeGrob(p1,p2,p3,p4,p5,p6,ncol=2,nrow=3),mylegend,nrow=2,heights=c(8.5,1.5))
