### compute the estimated distance matrix L
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

### compute the estimated distance matrix D, based on global dissimilarity
D_hat = function(theta_list){
  K = length(theta_list)
  D_mat = matrix(0,K,K)
  for (k in 1:K){
    for (l in 1:K){
      D_mat[k,l] = sum((theta_list[[k]] - theta_list[[l]])^2)
    }
  }
  return(D_mat)
}


### calculate se of hat D
se_D_hat = function(theta_list,V_list,n_list){
  K = length(theta_list)
  se_D_mat = matrix(0,K,K)
  for (k in 1:K){
    for (l in 1:K){
      gamma_hat = theta_list[[k]] - theta_list[[l]]
      se_D_mat[k,l] = sqrt(4 * t(gamma_hat) %*% V_list[[l]] %*% gamma_hat + 4 * t(gamma_hat) %*% V_list[[k]] %*% gamma_hat + 1/min(n_list[l],n_list[k]))
    }
  }
  return(se_D_mat)
}

### median estimator and CI based on parametric bootstrap
median_CI = function(est,se, median_B){
  K = length(est)
  median_est = median(est)
  median_se = sd(replicate(median_B, median(rnorm(n = K, mean = est, sd = se))))
  CI = median_est + c(-1,1)*qnorm(0.975)*median_se
  return(c(median_est,CI))
}

### inverse variance weighting
IVW = function(tau, se.tau, est.prevail, alpha.level = 0.05){
  est.IVW = sum( tau[est.prevail] / (se.tau[est.prevail])^2 ) / sum( 1/(se.tau[est.prevail])^2 )
  se.IVW = 1 / sqrt( sum( 1/(se.tau[est.prevail])^2 ) )
  CI.IVW = est.IVW + c(-1,1)*qnorm(1-alpha.level/2)*se.IVW
  return(c(est.IVW,CI.IVW))
}

### generate the voting matrix using one set of resampled dissimilarity
gen_H_m = function(Dhat_m, se_Dhat, Lhat_m, se_Lhat, rho = 0.3, alpha = 0.05, univariate = FALSE) {
  K = dim(Lhat_m)[1]
  if (univariate == FALSE){
    z = qnorm(1-alpha/(2*K*(K-1)))
  } else {
    z = qnorm(1-alpha/(K*(K-1)))
  }
  Shat_m = pmax(abs(Dhat_m/se_Dhat),abs(Lhat_m/se_Lhat))
  H_m = 1*(Shat_m <= rho*z)
  diag(H_m) = 1
  return(H_m)
}

# Generate \hat{V}[m] based on maximum clique of voting matrix
#   and \tilde{V}[m] based on rho sum of voting matrix
gen_V_m = function(Dhat_m, se_Dhat, Lhat_m, se_Lhat, rho = 0.3, alpha = 0.05, univariate = FALSE) {
  K = dim(Lhat_m)[1]
  H_m = gen_H_m(Dhat_m, se_Dhat, Lhat_m, se_Lhat, rho, alpha,univariate)
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


# Generate \hat{D}[m] or \hat{L}[m]
gen_disthat_m = function(dist, se_dist) {
  dist_hat = matrix( rnorm(length(dist), mean=as.vector(dist), sd=as.vector(se_dist)),
                     nrow(dist), nrow(dist) )
  diag(dist_hat) = 0
  dist_hat[lower.tri(dist_hat)] = t(dist_hat)[lower.tri(t(dist_hat))]
  #Dhat[upper.tri(Dhat)] = t(Dhat)[upper.tri(t(Dhat))]
  return(dist_hat)
}

### wrapper function to do one resampling and generate estimated prevailing set
RIFL.one = function(Dhat, se_Dhat, Lhat, se_Lhat, rho = 0.3, alpha = 0.05, univariate = FALSE){
  if (univariate == TRUE){
    Lhat_m = gen_disthat_m(Lhat,se_Lhat)
    Dhat_m = Lhat_m
  } else {
    Lhat_m = gen_disthat_m(Lhat,se_Lhat)
    Dhat_m = gen_disthat_m(Dhat,se_Dhat)
  }
  V_m = gen_V_m(Dhat_m, se_Dhat, Lhat_m, se_Lhat, rho, alpha, univariate)
  Vhat_m = V_m[1,]
  Vtilde_m = V_m[2,]
  return(c(Vhat_m,Vtilde_m))
}