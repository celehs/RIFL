### console inputs: 
###  S = sample size 500, 1000, 2000
###  M = separation level 1,2,3,4,5
###  W = number of majority sites 6 or 8
###  Q = seed 1 to 25
###  P = choice of nu, use value 4 to get nu = 0.8

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


########################
### helper functions ###
########################
source("tate.R")
source("helper/helper_fns_uni.R")


# do one bootstrap for one source dataset
one.boot.parametric = function(source_data_original,target_data,nu,method, or.formula, ps.formula, gamma.formula, link){
  source_est_boot = rep(NA,K)
  source_se_boot = rep(NA,K)
  m_list = rep(NA,K)
  for (k in 1:K){
    source_data_full = source_data_original[[k]]
    nn = dim(source_data_full)[1]
    if (method == "power"){
      m = floor(nn^nu)
    } else if (method == "multiple"){
      m = ceiling(0.5^nu*nn)
    }
    m_list[k] = m
    index = sample(1:nn,m,replace=TRUE)
    source_data_boot = source_data_full[index,]
    estimates = tate(source_data_boot, target_data, or.formula, ps.formula, gamma.formula, link)
    source_est_boot[k] = estimates[1]
    source_se_boot[k] = estimates[2]
  }
  
  Lhat_boot = L_hat(source_est_boot)
  se_Lhat_boot = se_L_hat(source_se_boot)
  
  MC_boot = gen_V_m(Lhat_boot, se_Lhat_boot, rho = 1, alpha = 0.05)[1,]
  IVW(source_est_boot,source_se_boot,which(MC_boot==1))[1]
}

do.one.sim = function(nu,seed,method,or.formula, ps.formula, gamma.formula, link, BB = 1000){
  data_original = vector("list",K)
  set.seed(seed)
  
  source_est = rep(NA,K)
  source_se = rep(NA,K)
  
  target.data = mvrnorm(N, mu = mean.target, Sigma = Cov)
  target.data = data.frame(target.data)
  colnames(target.data) = sapply(c(1:p), function(x){paste("X",x,sep="")})
  
  source.data.list = vector("list",K)
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
      
      tate.site = tate(source.data,target.data, or.formula, ps.formula, gamma.formula,link)
      source_est[k] = tate.site[1]
      source_se[k] = tate.site[2]
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
      tate.site = tate(source.data,target.data, or.formula, ps.formula, gamma.formula,link)
      source_est[k] = tate.site[1]
      source_se[k] = tate.site[2]
    }
  }
  
  Lhat = L_hat(source_est)
  se_Lhat = se_L_hat(source_se)
  MC = gen_V_m(Lhat, se_Lhat, rho = 1, alpha = 0.05)[1,]
  est.MC = IVW(source_est,source_se,which(MC==1))[1]
  
  res.MC.boot = replicate(BB,one.boot.parametric(source.data.list,target.data,nu,method, or.formula, ps.formula, gamma.formula, link))
  if (method == "power"){
    m = floor(n_list[1]^nu)
  } else if (method == "multiple"){
    m = ceiling(0.5^nu*n_list[1])
  }
  t.hat = quantile(sqrt(m)*(res.MC.boot - est.MC),c(0.975,0.025),na.rm=T)
  CI.MNB = est.MC -t.hat/sqrt(n_list[1])
  return(c(est.MC,CI.MNB))
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
Lmaj = P
if (Lmaj == 6){
  trt.effects = c(rep(-1,6),rep(-1-0.2*a,2),rep(-1-0.1*a,2))
  true.prevail = c(1:6)
}

if (Lmaj == 8){
  trt.effects = c(rep(-1,8),rep(-1-0.2*a,1),rep(-1-0.1*a,1))
  true.prevail = c(1:8)
}

K = 10
nsim = 20
nu_list = c(0.5,0.6,0.7,0.8,0.9,1)
n_list = rep(S,10)

allres = matrix(0,nsim,3)
for (i in 1:nsim){
  allres[i,] = do.one.sim(nu_list[P],((Q-1)*nsim+i)*3,"power",
                          "Y~.","A~X1+X2","Y~.","identity",500)
}

filename = paste("n",S,"diff",M,"nu",P,"seed",Q,"tateMNBmaj",W,".RData",sep="") #change
save(allres,file=filename)

for (n in c(500,1000,2000)){
  for (m in c(1:5)){
    for (pp in c(6,8)){
      res = NULL
      for (q in c(1:25)){
        filename = paste("n",n,"diff",m,"nu4seed",q,"tateMNBmaj",pp,".RData",sep="")
        load(filename)
        res = rbind(res,allres)
      }
      allres = res
      filename = paste("n",n,"diff",m,"nu4tateMNBmaj",pp,".RData",sep="")
      save(allres,file=filename)
    }
  }
}
