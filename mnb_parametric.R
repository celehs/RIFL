##############################################################
### implements M out of n bootstrap for parametric examples ##
##############################################################

### console inputs:
###   N sample size, 500, 1000 or 2000
###   M separation level, 1,2,3,4,5
###   Q random seed, 1,2,3,4,5
###   S value of nu, use 4 to get nu = 0.8
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
source('helper/helper_fns.R')

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
a = 0.1*M
nsim = 100 
L = 10
p = 10
L_maj = 8
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
  true.prevail = c(1:8)
}

betas = cbind(beta1,beta2,beta3,beta4,beta5,beta6,beta7,beta8,beta9,beta10)
intercepts = c(0.05,-0.05,0.1,-0.1,0.05,-0.05,0.1,-0.1,0,0)
Cov <- A1gen(rho,dim(betas)[1])
n_list = rep(N,10)
nu_list = c(0.5,0.6,0.7,0.8,0.9,1)


get_beta = function(yy,xx){
  mod = glm(yy~xx,family = binomial)
  list(mod$coefficients[-1],vcov(mod)[-1,-1])
}

# do one bootstrap for one source dataset
one.boot.parametric = function(data_original,nu,method){
  K = length(data_original)
  source_est_boot = rep(NA,K)
  source_se_boot = rep(NA,K)
  theta_list_boot = vector("list",K)
  V_list_boot = vector("list",K)
  m_list = rep(NA,K)
  for (k in 1:K){
    y = data_original[[k]][[1]]
    X = data_original[[k]][[2]]
    nn = length(y)
    if (method == "power"){
      m = floor(nn^nu)
    } else if (method == "multiple"){
      m = ceiling(0.5^nu*nn)
    }
    m_list[k] = m
    index = sample(1:nn,m,replace=TRUE)
    estimates = get_beta(y[index],X[index,])
    source_est_boot[k] = estimates[[1]][1]
    theta_list_boot[[k]] = estimates[[1]]
    source_se_boot[k] = sqrt(estimates[[2]][1,1])
    V_list_boot[[k]] = estimates[[2]]
  }
  
  Lhat_boot = L_hat(source_est_boot)
  se_Lhat_boot = se_L_hat(source_se_boot)
  Dhat_boot = D_hat(theta_list_boot)
  se_Dhat_boot = se_D_hat(theta_list_boot,V_list_boot,m_list)
  
  MC_boot = gen_V_m(Dhat_boot, se_Dhat_boot, Lhat_boot, se_Lhat_boot, rho = 1, alpha = 0.05, univariate = FALSE)[1,]
  IVW(source_est_boot,source_se_boot,which(MC_boot==1))[1]
}

do.one.sim = function(nu,seed,method,BB = 1000){
  data_original = vector("list",L)
  set.seed(seed)
  
  theta_list = vector("list",L)
  V_list = vector("list",L)
  source_est = rep(NA,L)
  source_se = rep(NA,L)
  
  for (l in 1:L){
    X = mvrnorm(n=n_list[l],mu=rep(0,p),Sigma=Cov)
    exp_val = intercepts[l] + X %*% betas[,l]
    prob = exp(exp_val)/(1+exp(exp_val)) 
    y = rbinom(n_list[l], 1, prob)
    data_original[[l]] = list(y,X)
    
    mod = glm(y~X,family=binomial())
    theta_list[[l]] = mod$coefficients[-1]
    V_list[[l]] = vcov(mod)[-1,-1]
    source_est[l] = mod$coefficients[2]
    source_se[l] = sqrt(V_list[[l]][1,1])
  }
  
  Lhat = L_hat(source_est)
  se_Lhat = se_L_hat(source_se)
  Dhat = D_hat(theta_list)
  se_Dhat = se_D_hat(theta_list,V_list,n_list)
  MC = gen_V_m(Dhat, se_Dhat, Lhat, se_Lhat, rho = 1, alpha = 0.05, univariate = FALSE)[1,]
  est.MC = IVW(source_est,source_se,which(MC==1))[1]
  
  res.MC.boot = replicate(BB,one.boot.parametric(data_original,nu,method))
  if (method == "power"){
    m = floor(n_list[1]^nu)
  } else if (method == "multiple"){
    m = ceiling(0.5^nu*n_list[1])
  }
  t.hat = quantile(sqrt(m)*(res.MC.boot - est.MC),c(0.975,0.025),na.rm=T)
  CI.MNB = est.MC -t.hat/sqrt(n_list[1])
  return(c(est.MC,CI.MNB))
}

allres = matrix(0,nsim,3)
for (i in 1:nsim){
  allres[i,] = do.one.sim(nu_list[S],((Q-1)*100+i)*9,"power") 
}
filename = paste("n",N,"diff",M,"nu",S,"seed",Q,"paramMNBmaj8.RData",sep="") 
save(allres,file=filename)

#######################
### process results ###
#######################
for (n in c(500,1000,2000)){
  for (m in c(1:5)){
    res = NULL
    for (q in c(1:5)){
      filename = paste("n",n,"diff",m,"nu4seed",q,"paramMNBmaj8.RData",sep="")
      load(filename)
      res = rbind(res,allres)
    }
    allres = res
    filename = paste("n",n,"diff",m,"nu4paramMNBmaj8.RData",sep="")
    save(allres,file=filename)
  }
}