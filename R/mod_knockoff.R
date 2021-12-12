library(here)
library(knockoff)

source(here("R", "kn_utils.R"))


# modified from the knockoff R package
# create knockoff matrix
create.fixed.causal <- function(X, ignore, sigma=NULL, y=NULL, randomize=F) {
  method = 'sdp'
  
  # Validate dimensions, if using fixed-X knockoffs
  n = nrow(X); p = ncol(X)
  if (n <= p)
    stop('Input X must have dimensions n > p')
  else if (n < 2*p) {
    warning('Input X has dimensions p < n < 2p. ',
            'Augmenting the model with extra rows.',immediate.=T)
    X.svd = svd(X, nu=n, nv=0)
    u2 = X.svd$u[,(p+1):n]
    X = rbind(X, matrix(0, 2*p-n, p))
    if( is.null(sigma) ) {
      if( is.null(y) ) {
        stop('Either the noise level "sigma" or the response variables "y" must
             be provided in order to augment the data with extra rows.')
      }
      else{
        sigma = sqrt(mean((t(u2) %*% y)^2)) # = sqrt(RSS/(n-p))
      }
    }
    if (randomize)
      y.extra = rnorm(2*p-n, sd=sigma)
    else
      y.extra = with_seed(0, rnorm(2*p-n, sd=sigma))
    y = c(y, y.extra)
  }
  # Normalize X, if using fixed-X knockoffs
  X = normc(X, center=F)
  
  Xk = create_sdp_cause(X, ignore, randomize)
  
  structure(list(X=X, Xk=Xk, y=y), class='knockoff.variables')
}

# modified from the knockoff R package
create_sdp_cause <- function(X, ignore, randomize) {
  # Compute SVD and U_perp.
  X.svd = decompose(X, randomize)
  
  # Check for rank deficiency.
  tol = 1e-5
  d = X.svd$d
  d_inv = 1 / d
  d_zeros = d <= tol*max(d)
  if (any(d_zeros)) {
    warning(paste('Data matrix is rank deficient.',
                  'Model is not identifiable, but proceeding with SDP knockoffs'),immediate.=T)
    d_inv[d_zeros] = 0
  }
  
  # Compute the Gram matrix and its (pseudo)inverse.
  G = (X.svd$v %*diag% d^2) %*% t(X.svd$v)
  G_inv = (X.svd$v %*diag% d_inv^2) %*% t(X.svd$v)
  
  # Optimize the parameter s of Equation 1.3 using SDP.
  s = create.solve_sdp(G)
  s[s <= tol] = 0
  s[ignore] <- 0
  
  # Construct the knockoff according to Equation 1.4:
  C.svd = canonical_svd(2*diag(s) - (s %diag*% G_inv %*diag% s))
  X_ko = X - (X %*% G_inv %*diag% s) + 
    (X.svd$u_perp %*diag% sqrt(pmax(0, C.svd$d))) %*% t(C.svd$v)
}


# modified from the knockoff R package
stat.glmnet_lambdasmax.cause <- function(X, X_k, y, family='gaussian', ...) {
  # Randomly swap columns of X and Xk
  swap = rbinom(ncol(X),1,0.5)
  swap.M = matrix(swap,nrow=nrow(X),ncol=length(swap),byrow=TRUE)
  X.swap  = X * (1-swap.M) + X_k * swap.M
  Xk.swap = X * swap.M + X_k * (1-swap.M)
  
  # Compute statistics
  Z = lasso_max_lambda_glmnet(cbind(X.swap, Xk.swap), y, method='glmnet', family=family, ...)
  p = ncol(X)
  orig = 1:p
  W = pmax(Z[orig], Z[orig+p]) * sign(Z[orig] - Z[orig+p])
  
  upper_level <- which(colSums(abs(X - X_k)) < 1e-8)
  W[upper_level] <- 0
  
  # Correct for swapping of columns of X and Xk
  W = W * (1-2*swap)
}



causal_discover.kn <- function(X, y, levels, alpha = 0.2){
  level_values <- sort(unique(levels))
  n_level <- length(level_values)
  level_ends <- sapply(level_values, function(level_value){
    max(which(levels <= level_value))
  })
  
  causes <- NULL
  for(level_i in 1:n_level){
    X_level <- X[, 1:level_ends[level_i]]
    
    if(level_i == 1) ignore <- NULL
    else ignore <- 1:level_ends[level_i-1]
    
    XXk <- create.fixed.causal(X_level, ignore)
    X_level <- XXk$X
    Xk_level_gene <- function(X){return(XXk$Xk)}
    
    # kn_select <- cknockoff::cknockoff(X = XXk$X, y = y,
    #                                   knockoffs = XXk$Xk, alpha = alpha / n_level)$selected
    kn_select <- knockoff.filter(X = X_level, y = y,
                                 knockoffs = Xk_level_gene,
                                 statistic = stat.glmnet_lambdasmax.cause)$selected
    kn_select <- kn_select[!(kn_select %in% ignore)]
    
    causes <- c(causes, kn_select)
  }
  
  return(causes)
}

