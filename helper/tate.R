### expit
expit = function(x){
  exp(x)/(1+exp(x))
}

### logit
logit = function(x){
  log(x/(1-x))
}

### estimating equation for the density ratio parameter gamma
density.ratio.ee = function(gamma, source.basis, target.basis){
  as.numeric(t(source.basis) %*% (exp(source.basis %*% gamma))/dim(source.basis)[1]) - apply(target.basis,2,mean)
}

### compute TATE
tate = function(source.data, target.data, or.formula, ps.formula, gamma.formula, link){
  ORformula = as.formula(or.formula)
  PSformula = as.formula(ps.formula)
  GAformula = as.formula(gamma.formula)
  n = dim(source.data)[1]
  N = dim(target.data)[1]
  p = dim(target.data)[2]
  
  trt.index = which(colnames(source.data) == "A")
  outcome.index = which(colnames(source.data) == "Y")
  
  # create design matrix for outcome regression on both source and target
  source.X = model.matrix(ORformula, data=source.data[,-trt.index])
  target.data.withY = data.frame(cbind(rep(1,N),target.data))
  colnames(target.data.withY)[1] = "Y"
  target.X = model.matrix(ORformula, data=target.data.withY)
  
  # create design matrix for propensity score on source
  source.X.ps = model.matrix(PSformula,data=source.data[,-outcome.index])
  
  # create the matrix of basis to mean match to get density ratio
  # on both source and target
  source.X.gamma = model.matrix(GAformula,data=source.data[,-trt.index])
  target.X.gamma = model.matrix(GAformula,data=target.data.withY)
  rm(target.data.withY)
  
  #source.X = as.matrix(cbind(rep(1,n),source.data[,-c(outcome.index,trt.index)]))
  #target.X = as.matrix(cbind(rep(1,N),target.data))
  
  # estimate the density ratio by matching the mean
  gamma.hat = nleqslv(x=rep(0,dim(source.X.gamma)[2]),fn=density.ratio.ee, 
                      source.basis=source.X.gamma,
                      target.basis = target.X.gamma)$x
  density.ratio = exp(as.matrix(source.X.gamma) %*% gamma.hat)
  
  # estimate outcome regression by glm
  if (link == "logit"){
    mod1 = glm(ORformula,data=source.data[source.data$A==1,-trt.index],family=binomial())
    mod0 = glm(ORformula,data=source.data[source.data$A==0,-trt.index],family=binomial())
  } else if (link == "identity"){
    mod1 = glm(ORformula,data=source.data[source.data$A==1,-trt.index],family=gaussian)
    mod0 = glm(ORformula,data=source.data[source.data$A==0,-trt.index],family=gaussian)
  }
  mu1.source = predict(mod1,newdata = source.data[-trt.index],type="response")
  mu0.source = predict(mod0,newdata = source.data[-trt.index],type="response")
  mu1.target = predict(mod1,newdata=target.data,type="response")
  mu0.target = predict(mod0,newdata=target.data,type="response")
  
  # estimate propensity score by glm
  mod.ps = glm(PSformula,data=source.data[,-outcome.index],family=binomial())
  propensity = predict(mod.ps, type="response")
  
  # TATE estimator
  augment.TATE = mean(density.ratio*(source.data$A/propensity*(source.data$Y - mu1.source) - (1-source.data$A)/(1-propensity)*(source.data$Y - mu0.source)))
  augment = mean(source.data$A/propensity*(source.data$Y - mu1.source) - (1-source.data$A)/(1-propensity)*(source.data$Y - mu0.source))
  DR.source = mean(mu1.source - mu0.source) + augment
  DR.TATE = mean(mu1.target - mu0.target) + augment.TATE
  
  # influence function of the finite dimensional parameter
  if (link == "logit"){
    Gprime = function(xx){
      xx*(1-xx)
    }
  } else if (link == "identity"){
    Gprime = function(xx){
      rep(1,length(xx))
    }
  }
  info.beta1 = t(source.X) %*% diag(Gprime(mu1.source)*source.data$A) %*% (source.X)/n
  IF.beta1 = solve(info.beta1) %*% t(source.X) %*% diag((source.data$Y-mu1.source)*source.data$A)
  
  info.beta0 = t(source.X) %*% diag(Gprime(mu0.source)*(1-source.data$A)) %*% (source.X)/n
  IF.beta0 = solve(info.beta0) %*% t(source.X) %*% diag((source.data$Y-mu0.source)*(1-source.data$A))
  
  info.gamma = t(source.X.gamma) %*% diag(as.numeric(exp(source.X.gamma %*% gamma.hat))) %*% source.X.gamma/n
  U.gamma = t(source.X.gamma) %*% diag(as.numeric(exp(source.X.gamma %*% gamma.hat)))
  IF.gamma = solve((-1)*info.gamma) %*% (U.gamma - apply(U.gamma,1,mean))
  IF.gamma.ignore = solve(info.gamma) %*% (t(target.X.gamma) - apply(target.X.gamma,2,mean))
  
  info.alpha = t(source.X.ps) %*% diag(propensity*(1-propensity)) %*% source.X.ps/n
  IF.alpha = solve(info.alpha) %*% t(source.X.ps) %*% diag(source.data$A - propensity)
  
  # partial derivative with respect to finite dimensional parameter
  delta.d.gamma = t(source.X.gamma) %*% (density.ratio*(source.data$A/propensity*(source.data$Y - mu1.source) - 
                                                          (1-source.data$A)/(1-propensity)*(source.data$Y - mu0.source)))/n
  delta.d.alpha = t(source.X.ps) %*% (density.ratio*((-1)*source.data$A*(1-propensity)/propensity*(source.data$Y - mu1.source)-
                                                       (1-source.data$A)*propensity/(1-propensity)*(source.data$Y - mu0.source)))/n
  delta.d.beta1 = t(source.X) %*% ((-1)*density.ratio*source.data$A/propensity*Gprime(mu1.source))/n
  delta.d.beta0 = t(source.X) %*% (density.ratio * (1-source.data$A)/(1-propensity)*Gprime(mu0.source))/n
  
  # IF of tate
  IF.M = t(IF.beta1) %*% (t(target.X) %*% Gprime(mu1.target))/N - t(IF.beta0) %*% (t(target.X) %*% Gprime(mu0.target))/N
  IF.M.ignore = (mu1.target - mu0.target) - mean(mu1.target - mu0.target)
  IF.TATE = IF.M + as.numeric(t(IF.gamma) %*% delta.d.gamma) +
    as.numeric(t(IF.alpha) %*% delta.d.alpha) +
    as.numeric(t(IF.beta1) %*% delta.d.beta1) +
    as.numeric(t(IF.beta0) %*% delta.d.beta0) +
    density.ratio*(source.data$A/propensity*(source.data$Y - mu1.source) - (1-source.data$A)/(1-propensity)*(source.data$Y - mu0.source)) - augment.TATE
  IF.ignore = IF.M.ignore + as.numeric(t(IF.gamma.ignore) %*% delta.d.gamma)
  # ignored the influence function pieces that come from the target data
  
  # return
  se = sqrt(var(IF.TATE)/n)
  se.full = sqrt(var(IF.TATE)/n + var(IF.ignore)/N)
  return(c(DR.TATE,se,se.full))
}
