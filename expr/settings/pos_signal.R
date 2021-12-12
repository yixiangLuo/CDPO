library(here)

library(knockoff)

source(here("expr", "utils.R"))


experiment <- "pos_signal"

p <- 100
n <- 4*p

sample_size <- 50
n_cores <- 7

pi1 <- 0.2

X_types <- c("Wide", "Moderate", "Long")  # "Wide", "Moderate", "Long"
X_titles <- paste0("X: ", X_types)

n_levels <- c(1, 4, 10)

fig_x_var <- list(name = "Signal Strength", value = c(0.05, 0.1, 0.2, 0.3)) # 0.05, 0.1, 0.2, 0.4


alpha <- 0.2
alphas <- rep(alpha, length(fig_x_var$value))

sig_sign <- 1
noise <- quote(rnorm(1))
# noise <- quote(rt(n, df = 1))













