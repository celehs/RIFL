# args <- commandArgs(TRUE)
# print(args)
# ## Parse the arguments
# if(length(args) == 0) {
#   print("No arguments supplied.")
# } else {
#   for (i in 1:length(args)) {
#     eval(parse(text = args[i]))
#   }
# }

N = 500 #sample size {500,1000,2000}
S = 5   #separation level {2,3,4,5,6,7}
P = 6   #number of majority site {6,8}
Q = 1   #job indexing, also used for seed {1,2,3,4,5,..25}

library(glmnet)
library(CVXR)
library(MASS)
library(igraph)

# import code for high-dimensional debiased inference
source('R/check.R')
source('R/helpers.R')
source('R/LF.R')
source('R/Methods.R')

# import RIFL helper functions
source('helper/helper_fns.R')

n = N
p = 500
L = 10
L_maj = P
nsim = 20

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
Cov <- (A1gen(rho,p))/2
mu <- rep(0,p)

beta1 = beta2 = beta3 = beta4 = beta5 = beta6 = beta7 =beta8 = beta9 = beta10 = rep(0,p)
if (L_maj == 6){
  beta1[1:11] = 0.1*seq(-5,5)
  beta2 = beta3 = beta4 = beta5 = beta6 = beta1
  beta7[1:5] = beta1[1:5]
  beta7[6:11] = beta1[6:11] + rep(0.15 + 0.05*S, 6)
  beta8 = beta7
  beta9[1:5] = beta1[1:5]
  beta9[6:11] = beta1[6:11] + rep(0.1 + 0.05*S, 6)
  beta10 = beta9
  true.prevail = c(1:6)
}

if (L_maj == 8){
  beta1[1:11] = 0.1*seq(-5,5)
  beta2 = beta3 = beta4 = beta5 = beta6 = beta7 = beta8 = beta1
  beta9[1:5] = beta1[1:5]
  beta9[6:11] = beta1[6:11] + rep(0.15 + 0.05*S, 6)
  beta10[1:5] = beta1[1:5]
  beta10[6:11] = beta1[6:11] + rep(0.1 + 0.05*S, 6)
  true.prevail = c(1:8)
}


beta = cbind(beta1,beta2,beta3,beta4,beta5,beta6,beta7,beta8,beta9,beta10)
intercepts = c(0.05,0.05,0.05,0,0,0,0.05,0,0.05,0)

median_B = 500
M = 500
thres = 0.5  # thres = 0.5 corresponds to the majority rule
             #  change to 0.7 to get RIFL with 80% rule
# the code here focus on theta_11 = 0.5
#   change uni.loading to uni.loading = matrix(c(rep(0,7),1,rep(0,p-8)),nrow=p,ncol=1)
#   to get results for theta_8 = 0.2
pair.truth.mat.list = pair.est.debias.mat.list = pair.est.plugin.mat.list =
  pair.se.mat.list = pair.last.term.mat.list = pair.est.trun.mat.list = rep(list(NA), nsim)
allres = matrix(NA,nsim,15)
allVmat = vector("list",nsim)
allsource.est = matrix(0,nsim,L)
allsource.se = matrix(0,nsim,L)
for(i.sim in 1:nsim){
  print(paste0('sim=',i.sim))
  seed = i.sim + (Q-1)*nsim
  set.seed(seed)
  
  X.list = vector("list",L)
  y.list = vector("list",L)
  source_est = rep(NA,L)
  source_se = rep(NA,L)
  for (l in 1:L){
    X1 = mvrnorm(n, mu, Cov)
    exp_val1 = X1%*%beta[,l] + intercepts[l]
    y1 = exp_val1 + rnorm(n)
    X.list[[l]] = X1
    y.list[[l]] = y1
    
    uni.loading = matrix(c(rep(0,10),1,rep(0,p-11)),nrow=p,ncol=1)
    LF.result = LF(X1,y1,uni.loading,"linear")
    source_est[l] = LF.result$est.debias.vec
    source_se[l] = LF.result$se.vec
  }
  
  pair.est.plugin.mat = pair.est.debias.mat = pair.est.trun.mat = pair.se.mat = pair.last.term.mat = matrix(0, nrow=L, ncol=L)
  for (l in 1:(L-1)){
    for (k in (l+1):L){
      X1 = X.list[[l]]
      y1 = y.list[[l]]
      X2 = X.list[[k]]
      y2 = y.list[[k]]
      
      idx.copy1 = sample(1:n, size=round(n/2), replace=F)
      idx.copy2 = setdiff(1:n, idx.copy1)
      ## training initial estimator ##
      outLas1 = cv.glmnet(X1[idx.copy1,], y1[idx.copy1], family='gaussian',alpha=1,
                          intercept=T, standardize=T)
      outLas2 = cv.glmnet(X2[idx.copy1,], y2[idx.copy1], family='gaussian',alpha=1,
                          intercept=T, standardize=T)
      beta.init1 = as.vector(coef(outLas1, s = outLas1$lambda.min))
      beta.init2 = as.vector(coef(outLas2, s = outLas2$lambda.min))
      
      ## bias correction ##
      loading = beta.init1 - beta.init2
      loading = loading[-1]
      loading.mat = as.matrix(loading)
      Est1 = LF(X1[idx.copy2,], y1[idx.copy2], loading.mat, model="linear", intercept=T, 
                intercept.loading=F, beta.init=beta.init1)
      Est2 = LF(X2[idx.copy2,], y2[idx.copy2], loading.mat, model="linear", intercept=T, 
                intercept.loading=F, beta.init=beta.init2)
      delta1 = Est1$est.debias.vec - Est1$est.plugin.vec
      delta2 = Est2$est.debias.vec - Est2$est.plugin.vec
      
      ## infos
      pair.est.plugin.mat[l,k] = sum(loading^2)
      pair.est.debias.mat[l,k] = sum(loading^2) + 2*delta1 - 2*delta2
      pair.est.trun.mat[l,k] = max(pair.est.debias.mat[l,k], 0)
      pair.se.mat[l,k] = sqrt(4*(Est1$se.vec)^2 + 4*(Est2$se.vec)^2 + 1/length(idx.copy2))
      pair.last.term.mat[l,k] = sum((loading - (beta[,l]-beta[,k]))^2)
    }
  }
  for (l in 1:(L-1)){
    for (k in (l+1):L){
      pair.est.plugin.mat[k,l] = pair.est.plugin.mat[l,k]
      pair.est.debias.mat[k,l] = pair.est.debias.mat[l,k]
      pair.est.trun.mat[k,l] = pair.est.trun.mat[l,k]
      pair.se.mat[k,l] = pair.se.mat[l,k]
      pair.last.term.mat[k,l] = pair.last.term.mat[l,k] 
    }
  }
  
  pair.truth.mat = matrix(NA, nrow=L, ncol=L)
  for(l in 1:L){
    for(k in 1:L){
      pair.truth.mat[l,k] = sum((beta[,l] - beta[,k])^2)
    }
  }
  
  Dhat = pair.est.trun.mat
  se_Dhat = pair.se.mat
  diag(se_Dhat) = 1
  Lhat = L_hat(source_est)
  se_Lhat = se_L_hat(source_se)
  
  ### true oracle, Xiudi implemented this 11/19
  res.oracle = IVW(source_est,source_se,true.prevail)
  
  ### median estimator
  res.median = median_CI(source_est,source_se,median_B)
  
  ### maximum clique estimator
  MC = gen_V_m(Dhat, se_Dhat, Lhat, se_Lhat, rho = 1, alpha = 0.05, univariate = FALSE)[1,]
  res.MC = IVW(source_est,source_se,which(MC==1))
  
  ### RIFL
  for (myrho in seq(0.2,1.2,0.025)){
    set.seed(seed)
    V.mat = replicate(M,RIFL.one(Dhat=Dhat, se_Dhat=se_Dhat, Lhat = Lhat, 
                                 se_Lhat = se_Lhat, rho = myrho, alpha = 0.05/20, univariate = FALSE))
    Vhat.mat = V.mat[1:(dim(V.mat)[1]/2),]
    size = apply(Vhat.mat,2,sum)
    percent = sum(size>(dim(Lhat)[2]*thres)) / M
    #print(c(myrho,percent))
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
  
  res.five = c(res.oracle,res.median,res.MC,res.RIFL,res.RIFL.alt)
  
  ### storing infos for each simulation ###
  allres[i.sim,] = res.five
  allsource.est[i.sim,] = source_est
  allsource.se[i.sim,] = source_se
  allVmat[[i.sim]] = V.mat
  pair.est.plugin.mat.list[[i.sim]] = pair.est.plugin.mat
  pair.est.debias.mat.list[[i.sim]] = pair.est.debias.mat
  pair.est.trun.mat.list[[i.sim]] = pair.est.trun.mat
  pair.se.mat.list[[i.sim]] = pair.se.mat
  pair.last.term.mat.list[[i.sim]] = pair.last.term.mat
  pair.truth.mat.list[[i.sim]] = pair.truth.mat
}

filename = paste("n",N,"diff",S,"seed",Q,"maj",L_maj,"high.RData",sep="")
save(allres,allsource.est,allsource.se,allVmat,
     pair.est.plugin.mat.list,
     pair.est.debias.mat.list,
     pair.est.trun.mat.list,
     pair.se.mat.list,
     pair.last.term.mat.list,file=filename)

###########################
### process the results ###
###########################
for (n in c(500,1000,2000)){
  for (s in c(2:7)){
    for (p in c(6,8)){
      res = NULL
      source.est = NULL
      source.se = NULL
      Dhat.list = NULL
      se.Dhat.list = NULL
      for (q in c(1:25)){
        filename = paste("n",n,"diff",s,"seed",q,"maj",p,"high.RData",sep="")
        load(filename)
        res = rbind(res,allres)
        source.est = rbind(source.est,allsource.est)
        source.se = rbind(source.se,allsource.se)
        Dhat.list = append(Dhat.list,pair.est.trun.mat.list)
        se.Dhat.list = append(se.Dhat.list, pair.se.mat.list)
      }
      allres = res
      allsource.est = source.est
      allsource.se = source.se
      allDhat = Dhat.list
      allseDhat = se.Dhat.list
      filename = paste("n",n,"diff",s,"maj",p,"high.RData",sep="")
      save(allres,allsource.est,allsource.se,allDhat,allseDhat, file=filename)
    }
  }
}

truth = 0.5
result_summary = NULL
for (n in c(500,1000,2000)){
  for (m in c(2:7)){
    filename = paste("n",n,"diff",m,"maj8high.RData",sep="")
    load(filename)
    est = c()
    cov = c()
    width=c()
    for (i in 1:5){
      est = c(est,mean(allres[,(3*i-2)]))
      cov = c(cov,mean(allres[,(3*i-1)] < truth & allres[,(3*i)]>truth))
      width = c(width,mean(allres[,(3*i)] - allres[,(3*i-1)] ))
      
    }
    points = allres[,7]
    bias = points - truth
    se = (allres[,9] - allres[,8])/2/qnorm(0.975)
    sd = sqrt(var(bias))
    se.new = sd/mean(se)*se
    temp = mean(bias)/se.new
    cv_alpha = sqrt(qchisq(1-0.05,df=1, ncp=temp^2))
    CI.new = cbind(points - se.new*cv_alpha, points + se.new*cv_alpha)
    cov.new = mean((CI.new[,1] < truth)*(CI.new[,2] > truth))
    len.new = mean(2*se.new*cv_alpha)
    
    est = c(est,est[3])
    cov = c(cov,cov.new)
    width = c(width,len.new)
    res_summary = data.frame(est=est,
                             cov=cov,
                             width=width,
                             n=rep(n,6),
                             separation = rep(m-1,6),
                             method=c("oracle","median","VMC","RIFL.hat","RIFL","OBA"))
    result_summary = rbind(result_summary,res_summary)
  }
}

result_summary = result_summary %>% filter(separation != 6 & method !="RIFL (adaptive)")

get_legend = function(myggplot){
  tmp = ggplot_gtable(ggplot_build(myggplot))
  leg = which(sapply(tmp$grobs,function(x) x$name) == "guide-box")
  legend = tmp$grobs[[leg]]
  return(legend)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(5)

library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggpubr)

datplot = result_summary %>% filter(n == 500 & method != "RIFL.hat" & method != "oracle")
p1 = ggplot(datplot,aes(x=separation,y=cov, group=method, col=method))+
  #ggtitle("n=500, coverage")+
  xlab("separation")+ylab("coverage")+geom_point()+geom_line()+ 
  geom_hline(yintercept=0.95, linetype="dashed")+
  scale_color_manual(values=cols[-2])+
  theme(legend.direction="horizontal",legend.title=element_blank())
mylegend = get_legend(p1)
p1 = p1+theme(legend.position="none")
p2 = ggplot(datplot,aes(x=separation,y=width, group=as.factor(method), col=as.factor(method)))+
  #ggtitle("n=500, length")+
  xlab("separation")+ylab("length")+geom_point()+geom_line()+theme(legend.position="none")+
  geom_hline(yintercept=result_summary$width[result_summary$method=="oracle" & result_summary$n == 500][1], linetype="dashed")+
  scale_color_manual(values=cols[-2])

p12 = ggarrange(p1,p2, ncol=2, nrow=1)
p12 = annotate_figure(p12, top = text_grob("n=500", 
                                           color = "black", face = "bold", size = 10))

#grid.arrange(arrangeGrob(p1,p2,ncol=2,nrow=1),mylegend,nrow=2,heights=c(9,1))


datplot = result_summary %>% filter(n == 1000 & method != "RIFL.hat" & method != "oracle")
p3 = ggplot(datplot,aes(x=separation,y=cov, group=as.factor(method), col=as.factor(method)))+
  #ggtitle("n=1000, coverage")+
  xlab("separation")+ylab("coverage")+geom_point()+geom_line()+ 
  geom_hline(yintercept=0.95, linetype="dashed")+
  scale_color_manual(values=cols[-2])+
  theme(legend.direction="horizontal",legend.title=element_blank())
#mylegend = get_legend(p1)
p3 = p3+theme(legend.position="none")
p4 = ggplot(datplot,aes(x=separation,y=width, group=as.factor(method), col=as.factor(method)))+
  #ggtitle("n=1000, length")+
  xlab("separation")+ylab("length")+geom_point()+geom_line()+theme(legend.position="none")+
  geom_hline(yintercept=result_summary$width[result_summary$method=="oracle" & result_summary$n == 1000][1], linetype="dashed")+
  scale_color_manual(values=cols[-2])
p34 = ggarrange(p3,p4, ncol=2, nrow=1)
p34 = annotate_figure(p34, top = text_grob("n=1000", 
                                           color = "black", face = "bold", size = 10))

#grid.arrange(arrangeGrob(p1,p2,ncol=2,nrow=1),mylegend,nrow=2,heights=c(9,1))

datplot = result_summary %>% filter(n == 2000 & method != "RIFL.hat" & method != "oracle")
p5 = ggplot(datplot,aes(x=separation,y=cov, group=as.factor(method), col=as.factor(method)))+
  #ggtitle("n=2000, coverage")+
  xlab("separation")+ylab("coverage")+geom_point()+geom_line()+ 
  geom_hline(yintercept=0.95, linetype="dashed")+
  scale_color_manual(values=cols[-2])+
  theme(legend.direction="horizontal",legend.title=element_blank())
#mylegend = get_legend(p1)
p5 = p5+theme(legend.position="none")
p6 = ggplot(datplot,aes(x=separation,y=width, group=as.factor(method), col=as.factor(method)))+
  #ggtitle("n=2000, length")+
  xlab("separation")+ylab("length")+geom_point()+geom_line()+theme(legend.position="none")+
  geom_hline(yintercept=result_summary$width[result_summary$method=="oracle" & result_summary$n == 2000][1], linetype="dashed")+
  scale_color_manual(values=cols[-2])
p56 = ggarrange(p5,p6, ncol=2, nrow=1)
p56 = annotate_figure(p56, top = text_grob("n=2000", 
                                           color = "black", face = "bold", size = 10))

#grid.arrange(arrangeGrob(p1,p2,p3,p4,p5,p6,ncol=2,nrow=3),mylegend,nrow=2,heights=c(9,1))
grid.arrange(arrangeGrob(p12,p34,p56,ncol=1,nrow=3),mylegend,nrow=2,heights=c(9,1))

