# Imported from knockoff R package


# https://github.com/msesia/knockoff-filter/blob/master/R/knockoff/R/util.R

# Fast versions of diag(d) %*% X and X %*% diag(d).
`%diag*%` <- function(d, X) d * X
`%*diag%` <- function(X, d) t(t(X) * d)

# Efficient test for matrix positive-definiteness
# 
# Computes the smallest eigenvalue of a matrix A to verify whether
# A is positive-definite
#' @keywords internal
is_posdef = function(A, tol=1e-9) {
  p = nrow(matrix(A))
  
  if (p<500) {
    lambda_min = min(eigen(A)$values)
  }
  else {
    oldw <- getOption("warn")
    options(warn = -1)
    lambda_min = RSpectra::eigs(A, 1, which="SM", opts=list(retvec = FALSE, maxitr=100, tol))$values
    options(warn = oldw)
    if( length(lambda_min)==0 ) {
      # RSpectra::eigs did not converge. Using eigen instead."
      lambda_min = min(eigen(A)$values)
    }
  }
  return (lambda_min>tol*10)
}

# Reduced SVD with canonical sign choice.
# 
# Our convention is that the sign of each vector in U is chosen such that the
# coefficient with the largest absolute value is positive.
#' @keywords internal
canonical_svd = function(X) {
  X.svd = tryCatch({
    svd(X)
  }, warning = function(w){}, error = function(e) {
    stop("SVD failed in the creation of fixed-design knockoffs. Try upgrading R to version >= 3.3.0")
  }, finally = {})
  
  for (j in 1:min(dim(X))) {
    i = which.max(abs(X.svd$u[,j]))
    if (X.svd$u[i,j] < 0) {
      X.svd$u[,j] = -X.svd$u[,j]
      X.svd$v[,j] = -X.svd$v[,j]
    }
  }
  return(X.svd)
}

# Scale the columns of a matrix to have unit norm.
#' @keywords internal
normc = function(X,center=T) {
  X.centered = scale(X, center=center, scale=F)
  X.scaled = scale(X.centered, center=F, scale=sqrt(colSums(X.centered^2)))
  X.scaled[,] # No attributes
}

# Generate a random matrix with i.i.d. normal entries.
#' @keywords internal
rnorm_matrix = function(n, p, mean=0, sd=1) {
  matrix(rnorm(n*p, mean, sd), nrow=n, ncol=p)
}

# Generate a random, sparse regression problem.
#' @keywords internal
random_problem = function(n, p, k=NULL, amplitude=3) {
  if (is.null(k)) k = max(1, as.integer(p/5))
  X = normc(rnorm_matrix(n, p))
  nonzero = sample(p, k)
  beta = amplitude * (1:p %in% nonzero)
  y.sample <- function() X %*% beta + rnorm(n)
  list(X = X, beta = beta, y = y.sample(), y.sample = y.sample)
}

# Evaluate an expression with the given random seed, then restore the old seed.
#' @keywords internal
with_seed = function(seed, expr) {
  seed.old = if (exists('.Random.seed')) .Random.seed else NULL
  set.seed(seed)
  on.exit({
    if (is.null(seed.old)) {
      if (exists('.Random.seed')) 
        rm(.Random.seed, envir=.GlobalEnv)
    } else {
      .Random.seed <<- seed.old
    }
  })
  expr
}



# https://github.com/msesia/knockoff-filter/blob/master/R/knockoff/R/solve_sdp.R
create.solve_sdp <- function(Sigma, gaptol=1e-6, maxit=1000, verbose=FALSE) {
  # Check that covariance matrix is symmetric
  stopifnot(isSymmetric(Sigma))
  # Convert the covariance matrix to a correlation matrix
  G = cov2cor(Sigma)
  p = dim(G)[1]
  
  # Check that the input matrix is positive-definite
  if (!is_posdef(G)) {
    warning('The covariance matrix is not positive-definite: knockoffs may not have power.', immediate.=T)
  }
  
  # Convert problem for SCS
  
  # Linear constraints
  Cl1 = rep(0,p)
  Al1 = -Matrix::Diagonal(p)
  Cl2 = rep(1,p)
  Al2 = Matrix::Diagonal(p)
  
  # Positive-definite cone
  d_As = c(diag(p))
  As = Matrix::Diagonal(length(d_As), x=d_As)
  As = As[which(Matrix::rowSums(As) > 0),] 
  Cs = c(2*G)
  
  # Assemble constraints and cones
  A = cbind(Al1,Al2,As)
  C = matrix(c(Cl1,Cl2,Cs),1)
  K=NULL
  K$s=p
  K$l=2*p
  
  # Objective
  b = rep(1,p)
  
  # Solve SDP with Rdsdp
  OPTIONS=NULL
  OPTIONS$gaptol=gaptol
  OPTIONS$maxit=maxit
  OPTIONS$logsummary=0
  OPTIONS$outputstats=0
  OPTIONS$print=0
  if(verbose) cat("Solving SDP ... ")
  sol = Rdsdp::dsdp(A,b,C,K,OPTIONS)
  if(verbose) cat("done. \n")
  
  # Check whether the solution is feasible
  if( ! identical(sol$STATS$stype,"PDFeasible")) {
    warning('The SDP solver returned a non-feasible solution. Knockoffs may lose power.')
  }
  
  # Clip solution to correct numerical errors (domain)
  s = sol$y
  s[s<0]=0
  s[s>1]=1
  
  # Compensate for numerical errors (feasibility)
  if(verbose) cat("Verifying that the solution is correct ... ")
  psd = 0
  s_eps = 1e-8
  while ((psd==0) & (s_eps<=0.1)) {
    if (is_posdef(2*G-diag(s*(1-s_eps),length(s)),tol=1e-9)) {
      psd  = 1
    }
    else {
      s_eps = s_eps*10
    }
  }
  s = s*(1-s_eps)
  s[s<0]=0
  if(verbose) cat("done. \n")
  
  # Verify that the solution is correct
  if (all(s==0)) {
    warning('In creation of SDP knockoffs, procedure failed. Knockoffs will have no power.',immediate.=T)
  }
  
  # Scale back the results for a covariance matrix
  return(s*diag(Sigma))
}

#' Vectorize a matrix into the SCS format
#'  
#' @rdname vectorize_matrix
#' @keywords internal
create.vectorize_matrix = function(M) {
  # Scale the off-diagonal entries by sqrt(2)
  vectorized_matrix = M
  vectorized_matrix[lower.tri(M,diag=FALSE)] = M[lower.tri(M,diag=FALSE)] * sqrt(2)
  # Stack the lower triangular elements column-wise
  vectorized_matrix = vectorized_matrix[lower.tri(vectorized_matrix,diag=TRUE)]
}



# https://github.com/msesia/knockoff-filter/blob/master/R/knockoff/R/create_fixed.R
#' Compute the SVD of X and construct an orthogonal matrix U_perp such that U_perp * U = 0.
decompose <- function(X, randomize) {
  n = nrow(X); p = ncol(X)
  stopifnot(n >= 2*p)
  
  result = canonical_svd(X)
  Q = qr.Q(qr(cbind(result$u, matrix(0,n,p))))
  u_perp = Q[,(p+1):(2*p)]
  if (randomize) {
    Q = qr.Q(qr(rnorm_matrix(p,p)))
    u_perp = u_perp %*% Q
  }
  result$u_perp = u_perp
  result
}



# https://github.com/msesia/knockoff-filter/blob/master/R/knockoff/R/stats_glmnet.R
lasso_max_lambda_glmnet <- function(X, y, nlambda=500, intercept=T, standardize=T, ...) {
  if (!requireNamespace('glmnet', quietly=T))
    stop('glmnet is not installed', call.=F)
  
  # Standardize the variables
  if( standardize ){
    X = scale(X)
  }
  
  n = nrow(X); p = ncol(X)
  if (!methods::hasArg(family) ) family = "gaussian"
  else family = list(...)$family
  
  if (!methods::hasArg(lambda) ) {
    if( identical(family, "gaussian") ) {
      if(!is.numeric(y)) {
        stop('Input y must be numeric.')
      }
      # Unless a lambda sequence is provided by the user, generate it
      lambda_max = max(abs(t(X) %*% y)) / n
      lambda_min = lambda_max / 2e3
      k = (0:(nlambda-1)) / nlambda
      lambda = lambda_max * (lambda_min/lambda_max)^k
    }
    else {
      lambda = NULL
    }
  }
  
  fit <- glmnet::glmnet(X, y, lambda=lambda, intercept=intercept, standardize=F, standardize.response=F, ...)
  
  first_nonzero <- function(x) match(T, abs(x) > 0) # NA if all(x==0)
  if(family=="multinomial") {
    indices <- sapply(fit$beta, function(beta) apply(beta, 1, first_nonzero))
    indices <- apply(indices, 1, min)
  } else {
    indices <- apply(fit$beta, 1, first_nonzero)
  }
  names(indices) <- NULL
  ifelse(is.na(indices), 0, fit$lambda[indices] * n)
}



