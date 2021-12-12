create_adj_mat <- function(p, X_pi1 = 0.2, X_signal = quote(runif(1, 0.5, 1)),
                           XY_pi1 = 5/p, XY_signal = quote(runif(1, 0.5, 1)),
                           sig_sign = quote((2*rbinom(1, 1, 0.5)-1))){
  adj_mat <- matrix(0, nrow = p+1, ncol = p+1)
  upper_tri <- upper.tri(adj_mat)
  adj_mat[upper_tri] <- replicate(sum(upper_tri), rbinom(1, 1, X_pi1) * 
                                    eval(sig_sign) * eval(X_signal))
  
  adj_mat[, p+1] <- 0
  adj_mat[1:p, p+1] <- replicate(p, rbinom(1, 1, XY_pi1) * 
                                   eval(sig_sign) * eval(XY_signal))
  
  return(adj_mat)
}

get_causes_effect <- function(adj_mat, levels){
  p <- NROW(adj_mat) - 1
  causes_effect <- adj_mat[1:p, p+1]  # first count direct effects
  for(node_i in p:1){
    effects <- causes_effect[node_i] * adj_mat[1:(node_i-1), node_i]
    effects[levels[1:(node_i-1)] >= levels[node_i]] <- 0
    causes_effect[1:(node_i-1)] <- causes_effect[1:(node_i-1)] + c(effects)
  }
  causes_effect[abs(causes_effect) < 1e-6] <- 0
  
  return(causes_effect)
}

generate_data <- function(n, adj_mat, noise = quote(rnorm(1))){
  p <- NROW(adj_mat) - 1
  data <- matrix(NA, n, p+1)
  data[, 1] <- replicate(n, eval(noise))
  for(level_i in 2:(p+1)){
    data[, level_i] <- matrix(data[, 1:(level_i-1)], nrow = n) %*% adj_mat[1:(level_i-1), level_i] +
      replicate(n, eval(noise))
  }
  
  return(list(X = data[, 1:p], y = c(data[, p+1])))
}

generate_y <- function(X, adj_mat, noise = quote(rnorm(1))){
  p <- NROW(adj_mat) - 1
  n <- NROW(X)
  
  y <- X %*% adj_mat[1:p, p+1] + replicate(n, eval(noise))
  
  return(y)
}





