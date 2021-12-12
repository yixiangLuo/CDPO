library(abind)    # abind

library(foreach)
library(doParallel)


# compute fdr and power for methods on a linear problem
get_fdp_power <- function(n, p, levels, pi1, sig_sign, noise, alpha,
                          fig_x_var, sample_size, n_cores,
                          expr_name, X_title, X_seed){

    registerDoParallel(n_cores)

    org_data <- list()
    FDR_Power <- NULL
    
    # set.seed(X_seed)
    
    

    for(var_i in 1:length(fig_x_var$value)){
        fig_x_value <- fig_x_var$value[var_i]
        sign_strength <- fig_x_value

        update_progress(expr_name, X_title, sign_strength, sample_size, action = "start")

        results <- foreach(iter = 1:sample_size, .options.multicore = list(preschedule = F)) %dopar% {
        # results <- lapply(1:sample_size, function(iter){
            # print(iter)

            set.seed(iter)
            
            adj_mat <- create_adj_mat(p, X_pi1 = pi1, X_signal = sign_strength,
                                      XY_pi1 = pi1, XY_signal = sign_strength,
                                      sig_sign = sig_sign)
            
            true_effects <- get_causes_effect(adj_mat, levels)
            H0 <- true_effects == 0
            
            # y <- generate_y(X, adj_mat, noise)
            
            data <- generate_data(n, adj_mat, noise)
            X <- data$X
            y <- data$y
            # save(X, y, alpha, file = "debug.RData")

            fdp_power <- NULL
            method_names <- NULL
            
            pscore_mat <- sapply(1:p, function(treat_ind){
                treats <- c(X[, treat_ind] >= median(X[, treat_ind]))

                covar_indeces <- setdiff(which(levels <= levels[treat_ind]), treat_ind)
                covars <- X[, covar_indeces]

                prop_scores <- est_propensity(treats, covars)
                prop_scores <- prop_scores + 1e-10 * runif(n, min = -1, max = 1)

                return(prop_scores)
            })
            
            # Hajek
            selected <- causal_discover.ATE(X, y, levels, causal_infer = causal_infer.Hajek,
                                            alpha, pscore_mat)
            method_fdp_power <- calc_FDP_power(selected, H0)
            fdp_power <- cbind(fdp_power, method_fdp_power)
            method_names <- c(method_names, "Hajek")
            
            # DR
            selected <- causal_discover.ATE(X, y, levels, causal_infer = causal_infer.DR,
                                            alpha, pscore_mat)
            method_fdp_power <- calc_FDP_power(selected, H0)
            fdp_power <- cbind(fdp_power, method_fdp_power)
            method_names <- c(method_names, "DR")
            
            # Match
            selected <- causal_discover.ATE(X, y, levels, causal_infer = causal_infer.Match,
                                            alpha, pscore_mat)
            method_fdp_power <- calc_FDP_power(selected, H0)
            fdp_power <- cbind(fdp_power, method_fdp_power)
            method_names <- c(method_names, "Match")
            
            # knockoff
            selected <- causal_discover.kn(X, y, levels, alpha)
            method_fdp_power <- calc_FDP_power(selected, H0)
            fdp_power <- cbind(fdp_power, method_fdp_power)
            method_names <- c(method_names, "knockoff")
            
            print(iter)
            

            colnames(fdp_power) <- method_names

            update_progress(expr_name, X_title, sign_strength, sample_size, action = "progress")

            return(fdp_power)
        }

        update_progress(expr_name, X_title, sign_strength, sample_size, action = "end")

        method_names <- colnames(results[[1]])
        results <- abind(results, along=3)

        org_data[[var_i]] <- results

        interp_results <- function(results, method_names, type, alpha){
            results <- unname(results)
            interpretation <- data.frame(methods = factor(1:NROW(results), ordered = T),
                                         method_name = method_names,
                                         mean = rowMeans(results),
                                         std = sqrt(rowMeans((results - rowMeans(results))^2)),
                                         type = type,
                                         alpha = alpha,
                                         fig_x = fig_x_value)
        }

        fdr_power <- rbind(interp_results(results[1, , ], method_names, "FDR", alpha),
                           interp_results(results[2, , ], method_names, "Power", alpha))
        FDR_Power <- rbind(FDR_Power, fdr_power)
    }

    return(list(FDR_Power = FDR_Power, org_data = org_data))
}







