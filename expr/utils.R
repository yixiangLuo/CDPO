# library(mvtnorm)  # rmvnorm
library(here)

library(foreach)
library(doParallel)


corr_noise <- function(n, rho){
    cov_mat <- matrix(0, n, n)
    cov_mat[abs(row(cov_mat) - col(cov_mat)) <= 1] <- -rho * 0.5
    diag(cov_mat) <- 1
    
    R <- chol(cov_mat)
    noise <- t(R) %*% matrix(rnorm(n), nrow = n)
    
    return(noise)
}

makeup_vectors <- function(...){
    vars <- list(...)
    max_length <- max(sapply(vars, length))
    
    env <- parent.frame()
    for(var_i in 1:length(vars)){
        var_name <- names(vars[var_i])
        var_length <- length(vars[[var_i]])
        
        if(var_length < max_length){
            var_list <- list()
            for(i in 1:max_length) var_list <- c(var_list, vars[[var_i]])
            
            env[[var_name]] <- var_list
        } else{
            env[[var_name]] <- as.list(vars[[var_i]])
        }
    }
    
    invisible()
}


calc_FDP_power <- function(rejs, H0, sign_predict = NULL, sign_beta = NULL){
    nrejs <- length(rejs)

    false_discovery <- length(intersect(rejs, which(H0)))
    true_discovery <- nrejs - false_discovery

    FDP <- false_discovery / max(nrejs, 1)
    power <- true_discovery / max(sum(!H0), 1)

    if(!is.null(sign_predict)){
        false_dir <- sum((sign_predict * sign_beta)[rejs] <= 0)
        true_dir <- nrejs - false_dir

        FDP_dir <- false_dir / max(nrejs, 1)
        power_dir <- true_dir / length(sign_beta)
    } else{
        FDP_dir <- NA
        power_dir <- NA
    }

    return(c(FDP, power, FDP_dir, power_dir))
}



try_repeat <- function(expr, default = NULL, n_times = 100){
    success <- F
    for(iter in 1:n_times){
        tryCatch({
            eval(expr, envir = parent.frame())
            success <- T
        }, error = function(msg){}, warning = function(msg){})
        if(success) break
        Sys.sleep(0.01)
    }
    if(!success) eval(default, envir = parent.frame())
    
    return(success)
}

update_count <- function(expr_name, action){
    file_name <- here("data", "temp", paste0("progress-", expr_name, ".Rdata"))
    
    if(action == "start"){
        iters_done <- 0
        start_time <- Sys.time()
        save(iters_done, start_time, file = file_name)
    } else if(action == "progress"){
        success <- try_repeat(load(file = file_name))
        if(success){
            iters_done <- iters_done + 1
            try_repeat(save(iters_done, start_time, file = file_name))
        } else{
            iters_done <- "(can't access progress record)"
        }
    } else if(action == "end"){
        load(file = file_name)
        iters_done <- Sys.time() - start_time  # abuse of var name
        file.remove(file_name)
    }
    
    return(iters_done)
}

print_progress <- function(expr_name, X_title, variable, iters_done, sample_size, action){
    file_name <- here("data", "temp", paste0("progress-", expr_name, ".txt"))
    
    try_num <- ifelse(action == "progress", 100, 1)
    try_repeat(progress_record <- readLines(file_name), 
               default = {progress_record <- NULL},
               n_times = try_num)
    
    if(action == "start"){
        progress_record <- c(progress_record,
                             paste0("-------------- ", X_title, " --------------"),
                             paste0("variable: ", variable, " -- start: ", Sys.time()),
                             paste0("variable: ", variable, " -- 0/", sample_size, " done."))
        writeLines(progress_record, con = file_name)
    } else if(action == "progress"){
        if(length(progress_record) > 1){
            progress_record <- progress_record[1:(length(progress_record)-1)]
            progress_record <- c(progress_record,
                                 paste0("variable: ", variable, " -- ", iters_done, "/", sample_size, " done."))
            try_repeat(writeLines(progress_record, con = file_name))
        }
    } else if(action == "end"){
        runtime <- round(as.double(iters_done), digits = 2)
        time_unit <- units(iters_done)
        progress_record <- c(progress_record,
                             paste0("variable: ", variable, " -- end: ", Sys.time()),
                             paste("--------- runtime:", runtime, time_unit, "---------"),
                             "")
        writeLines(progress_record, con = file_name)
    }
    
    
}

update_progress <- function(expr_name, X_title, alpha, sample_size, action){
    iters_done <- update_count(expr_name, action)
    print_progress(expr_name, X_title, alpha, iters_done, sample_size, action)
}

# borrowed from knockoff
# Evaluate an expression with the given random seed, then restore the old seed
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


# convert a digit to a string
digit_to_char <- function(number){
    return(str_replace(as.character(number), "\\.", "d"))
}

# convert a string to a digit
char_to_digit <- function(char){
    return(str_replace(as.character(char), "d", "\\."))
}

char_to_digit <- function(str){
    gsub(paste0("([0-9]+)d([0-9]+)"),
         paste0("\\1.\\2"),
         str)
}

parse_name <- function(str){
    str <- char_to_digit(str)
    str <- str_replace(str, "_L_", "(")
    str <- str_replace(str, "_R_", ")")
    str <- str_replace(str, "_D_", "-")
    str <- str_replace(str, "_STAR", "*")
}

