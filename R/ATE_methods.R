library(Matching)

est_propensity <- function(treats, covars){
  # method <- match.arg(method)
  
  prop_scores_est <- glm(treats ~ covars, family = binomial(link = "logit"))$fitted.values
  
  prop_scores_est <- c(prop_scores_est)
  
  return(prop_scores_est)
}


check_covar_balance <- function(treats, covars, prop_scores){
  bin_num <- 10
  bin_breaks <- seq(from = 0, to = 1, length.out = bin_num+1)
  
  strata_belong <- sapply(prop_scores, function(prop_score){
    max(which(prop_score >= bin_breaks))
  })
  
  balance_difference <- sapply(1:bin_num, function(bin_ind){
    in_strata <- strata_belong == bin_ind
    treat_strata <- treats & in_strata
    control_strata <- !treats & in_strata
    
    difference <- colMeans(X[treat_strata, ]) - colMeans(X[control_strata, ])
    
    return(difference)
  })
  
  return(balance_difference)
}



ATE_by_Hajek <- function(treats, covars, outcomes, prop_scores){
  weights_treat <- treats / prop_scores
  weights_control <- (1-treats) / (1-prop_scores)
  
  ATE.Hajek <- sum(weights_treat * outcomes) / sum(weights_treat) - 
    sum(weights_control * outcomes) / sum(weights_control)
  
  return(ATE.Hajek)
}

ATE_by_DR <- function(treats, covars, outcomes, prop_scores){
  tryCatch({
    fitted_treat <- glm(outcomes ~ covars, weights = as.integer(treats),
                        family = gaussian)$fitted.values
    fitted_control <- glm(outcomes ~ covars, weights = as.integer(!treats),
                          family = gaussian)$fitted.values
  },
  error = function(msg){browser()})
  
  
  mean_treat <- mean(treats * (outcomes - fitted_treat) / prop_scores + fitted_treat)
  mean_control <- mean((!treats) * (outcomes - fitted_control) / (1-prop_scores) + fitted_control)
  
  ATE.DR <- mean_treat - mean_control
  
  return(ATE.DR)
}

est_sd_bootstrap <- function(estimator, treats, covars, outcomes, prop_scores,
                             boot_sample = 100){
  n <- length(treats)
  boot_ests <- replicate(boot_sample, {
    boot_indeces <- sample(1:n, n, replace = T)
    boot_est <- estimator(treats[boot_indeces], covars[boot_indeces, ],
                          outcomes[boot_indeces], prop_scores[boot_indeces])
    
    return(boot_est)
  })
  
  sd_est <- sd(boot_ests)
  
  return(sd_est)
}

causal_infer.Hajek <- function(treats, covars, outcomes, prop_scores){
  ATE <- ATE_by_Hajek(treats, covars, outcomes, prop_scores)
  sd <- est_sd_bootstrap(ATE_by_Hajek, treats, covars, outcomes, prop_scores)
  
  return(list(ATE = ATE, sd = sd))
}

causal_infer.DR <- function(treats, covars, outcomes, prop_scores){
  ATE <- ATE_by_DR(treats, covars, outcomes, prop_scores)
  sd <- est_sd_bootstrap(ATE_by_DR, treats, covars, outcomes, prop_scores)
  
  return(list(ATE = ATE, sd = sd))
}

causal_infer.Match <- function(treats, covars, outcomes, prop_scores){
  match.result <- Match(Y = outcomes, Tr = treats, X = prop_scores, estimand = "ATE",
                        M = 1, BiasAdjust = T)
  
  return(list(ATE = match.result$est, sd = match.result$se))
}





BH <- function(pvals, alpha, weights = rep(1, length(pvals)) / length(pvals)){
  n <- length(pvals)
  
  adjust_pvals <- sort(pvals / weights) / (1:n)
  nrejs <- max(0, which(adjust_pvals <= alpha))
  
  rejs <- which(pvals <= nrejs * alpha * weights)
  
  return(list(nrejs = nrejs, rejs = rejs))
}

causal_discover.ATE <- function(X, y, levels, causal_infer, alpha = 0.05, pscore_mat = NULL){
  p <- NCOL(X)
  n <- NROW(X)
  zvals <- sapply(1:p, function(treat_ind){
    treats <- c(X[, treat_ind] >= median(X[, treat_ind]))
    
    covar_indeces <- setdiff(which(levels <= levels[treat_ind]), treat_ind)
    covars <- matrix(X[, covar_indeces], nrow = n)
    
    outcomes <- y
    
    if(is.null(pscore_mat)){
      prop_scores <- est_propensity(treats, covars)
    } else{
      prop_scores <- c(pscore_mat[, treat_ind])
    }
    
    trunc_pscore <- 0.1
    kept <- which((prop_scores > trunc_pscore) & (prop_scores < 1-trunc_pscore))
    
    if(length(kept) > 10){
      est.result <- causal_infer(treats[kept], covars[kept, ], outcomes[kept], prop_scores[kept])
      zval <- est.result$ATE / est.result$sd
    } else{
      zval <- 0
    }
    
    
    
    return(zval)
  })
  
  pvals <- 2 * pnorm(abs(zvals), lower.tail = FALSE)
  causes <- BH(pvals, alpha)$rejs
  
  return(causes)
}








