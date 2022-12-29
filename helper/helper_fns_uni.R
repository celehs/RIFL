### inverse variance weighting
IVW = function(tau, se.tau, est.prevail, alpha.level = 0.05){
  est.IVW = sum( tau[est.prevail] / (se.tau[est.prevail])^2 ) / sum( 1/(se.tau[est.prevail])^2 )
  se.IVW = 1 / sqrt( sum( 1/(se.tau[est.prevail])^2 ) )
  CI.IVW = est.IVW + c(-1,1)*qnorm(1-alpha.level/2)*se.IVW
  return(c(est.IVW,CI.IVW))
}

### calculate the matrix hat L, hat beta l - hat beta k
L_hat = function(source_est){
  K = length(source_est)
  dist_mat = matrix(0,K,K)
  for (l in 1:K){
    for (k in 1:K){
      dist_mat[l,k] = source_est[l] - source_est[k]
    }
  }
  return(dist_mat)
}


### calculate se of hat L
se_L_hat = function(source_se){
  K = length(source_se)
  se_mat = matrix(0,K,K)
  for (k in 1:K){
    for (l in 1:K){
      se_mat[l,k] = sqrt(source_se[l]^2 + source_se[k]^2)
    }
  }
  #diag(se_mat) = 0
  return(se_mat)
}

# Generate \hat{L}[m]
gen_disthat_m = function(dist, se_dist) {
  dist_hat = matrix( rnorm(length(dist), mean=as.vector(dist), sd=as.vector(se_dist)),
                     nrow(dist), nrow(dist) )
  diag(dist_hat) = 0
  dist_hat[lower.tri(dist_hat)] = t(dist_hat)[lower.tri(t(dist_hat))]
  #Dhat[upper.tri(Dhat)] = t(Dhat)[upper.tri(t(Dhat))]
  return(dist_hat)
}

### calculate voting matrix for one resample
gen_H_m = function(Lhat_m, se_Lhat, rho = 0.3, alpha = 0.05) {
  K = dim(Lhat_m)[1]
  z = qnorm(1-alpha/(K*(K-1)))
  Shat_m = abs(Lhat_m/se_Lhat)
  H_m = 1*(Shat_m <= rho*z)
  diag(H_m) = 1
  return(H_m)
}

### calculate estimates of prevailing set for one resample
### two methods: Vhat = maximum clique, Vtilde = sites receive more than half votes
gen_V_m = function(Lhat_m, se_Lhat, rho = 0.3, alpha = 0.05) {
  K = dim(Lhat_m)[1]
  H_m = gen_H_m(Lhat_m, se_Lhat, rho, alpha)
  voting.graph = igraph::as.undirected(igraph::graph_from_adjacency_matrix(H_m))
  max.clique.graph = igraph::largest_cliques(voting.graph)
  temp = max.clique.graph[1]
  # temp = unique(igraph::as_ids(Reduce(c,max.clique.graph)))
  max.clique = sort(igraph::as_ids(Reduce(c,temp)))
  max.clique.set = rep(0,K)
  max.clique.set[max.clique] = 1
  
  Vtilde_index = as.numeric(rowSums(H_m) > K/2)
  return(rbind(max.clique.set,Vtilde_index))
}

### wrapper to do one resample
RIFL.one = function(Lhat, se_Lhat, rho = 0.3, alpha = 0.05){
  Lhat_m = gen_disthat_m(Lhat,se_Lhat)
  V_m = gen_V_m(Lhat_m, se_Lhat, rho, alpha)
  Vhat_m = V_m[1,]
  Vtilde_m = V_m[2,]
  return(c(Vhat_m,Vtilde_m))
}

# helper function: parametric bootstrap for median estimator
median_CI = function(est,se, median_B){
  K = length(est)
  median_est = median(est)
  median_resample = replicate(median_B, median(rnorm(n = K, mean = est, sd = se)))
  median_se = sd(median_resample)
  CI = median_est + c(-1,1)*qnorm(0.975)*median_se
  #CI = quantile(median_resample,c(0.025,0.975))
  return(c(median_est,CI))
}