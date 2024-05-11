# nolint start
setwd("C:/Users/Frank/OneDrive/Documents/Masters/Final Project/Empirical tests/R")

list.of.packages <- c("data.table", "dplyr", "magrittr", "ggplot2", "zoo", "caret", "digest", "Rcpp", "RcppParallel", "RcppArmadillo", "parallel", "reshape2", "CVXR", "quantmod", "matrixcalc", "pwr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load previously installed packages
list.of.packages <- c(list.of.packages, "htmltab")
for (package in list.of.packages) {
    if (package %in% installed.packages()) {
        library(package, character.only = TRUE)
    }
}
sourceCpp("../cpp/LPTN.cpp")

# Load historical stock prices
stock_prices <- fread("./Data/SPY_weekly_log_returns.csv")

date_range = range(stock_prices$Date)
date_from = date_range[1]
date_to = date_range[2]

# Plot S&P 500 historical price
getSymbols("^GSPC", from = date_from, to = date_to, adjust = TRUE)
SnP <- as.data.table(GSPC$GSPC.Adjusted)
names(SnP) <- c("Date", "Price")

min_y = as.integer(0.8*min(SnP$Price))
min_y = min_y - min_y %% 50

ggplot(SnP, aes(x=Date, y=Price)) +
    geom_area(fill="lightgreen") +
    geom_line(color="green", size=1) +
    coord_cartesian(ylim = c(min_y, max(SnP$Price, na.rm = TRUE))) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))


# Price matrix
prices <- as.matrix(stock_prices[, -1, with = FALSE])

# Function to return a random train and test set from the price matrix
sample_matrix <- function(sample_size=25, ndays=2*365, train_p = 0.8){

    start_from <- round(runif(1, min=1, max=nrow(prices)-ndays))
    keep <- seq(from = start_from, to = start_from + ndays - 1, by = 7)
    output <- prices[keep, ]

    keep <- which(!apply(output, 2, function(x){any(is.na(x))}))
    keep <- sample(keep, size = min(sample_size, length(keep)))
    output <- output[, keep]

    split_at <- floor(train_p * nrow(output))
    return(list(train=output[1:split_at, ], test=output[-c(1:split_at), ]))

}

data <- sample_matrix()
X <- data$train
mu <- apply(X, 2, mean)
vcov_matrix <- var(X)

# Function to get the best LPTN likelihood tau parameter
best_tau <- function(X, mu, vcov_mat, tol=0.005, max_iter=100, ncores=1){

    wrapper <- function(p_outlier){fit_LPTN_2(X, mu, vcov_mat, tau=sqrt(qchisq(1-p_outlier, ncol(X))), tol=tol, max_iter=max_iter, ncores=ncores)}
    best <- optimize(wrapper, interval = c(0, 0.995), maximum = TRUE)$maximum
    return(sqrt(qchisq(1-best, ncol(X))))

}

tau <- best_tau(X, mu, vcov_matrix)

# Function to get the LPTN maximum likelihood estimators with the best tau
LPTN_estimators <- function(X, mu, vcov_mat, tol=0.005, max_iter=100, ncores=1){

    tau <- best_tau(X, mu, vcov_mat, tol, max_iter, ncores)
    return(fit_LPTN(X, mu, vcov_mat, tau, tol, max_iter, ncores))

}

robust_estimators <- LPTN_estimators(X, mu, vcov_matrix)

# Mean variance optimisation
MVO <- function(expected_returns, cov_matrix, target_return=0.1/52, how="all"){

    n <- length(expected_returns)
    # Variables
    w <- Variable(n)  # Portfolio weights
    # Objective
    portfolio_variance <- quad_form(w, cov_matrix)
    objective <- Minimize(portfolio_variance)
    # Constraints

    optimal_weights <- rep(NA, length(expected_returns))
    if(how == "all"){
        constraints <- list(sum(w) == 0, # Sum to 0
                            t(w) %*% expected_returns >= target_return) # Expected returns must be greater or equal to target_return
        # Formulate and solve the problem
        problem <- Problem(objective, constraints)
        result <- solve(problem)
        optimal_weights <- result$getValue(w)

        # make sure the weights sum up to 0
        optimal_weights <- optimal_weights - sum(optimal_weights)
        # make sure the sum of abs(weight) sums up to 1
        optimal_weights <- optimal_weights / sum(abs(optimal_weights))

    }

    else if(how == "long"){
        constraints <- list(sum(w) == 1, # Sum to 1
                            w >= 0, # Weights >= 0
                            t(w) %*% expected_returns >= target_return) # Expected returns must be greater or equal to target_return
        # Formulate and solve the problem
        problem <- Problem(objective, constraints)
        result <- solve(problem)
        if(sum(abs(result$getValue(w))) == 0){
            print("Weights are all zero")
            return(optimal_weights)
        }
        # Optimal weights
        optimal_weights <- pmax(result$getValue(w), 0)
        optimal_weights <- optimal_weights/sum(optimal_weights)

    }

    else if(how == "sharpe"){
      precision_mat <- solve(cov_matrix)
      J <- rep(1, length(expected_returns))
      weights <- precision_mat %*% expected_returns
      if(sum(abs(weights)) == 0){
        print("Weights are all zero")
        return(optimal_weights)
      }
      weights <- weights / sum(abs(weights))
      expected_returns <- sum(weights * expected_returns)
      if(expected_returns <= 0){
        print("Weights are all zero")
        return(optimal_weights)
      }
      optimal_weights <- weights

    }

    return(optimal_weights)

}

MVO_results <- MVO(mu, vcov_matrix, mean(mu), how="all")
MVO_results <- MVO(mu, vcov_matrix, mean(mu), how="long")
MVO_results <- MVO(mu, vcov_matrix, mean(mu), how="sharpe")

is_solvable <- function(X) {
  # Attempt to solve the matrix X
  result <- tryCatch({
    solve(X)
    TRUE # If solve() succeeds, return TRUE
  }, error = function(e) {
    FALSE # If an error occurs (e.g., matrix is singular), return FALSE
  })
  
  return(result)
}

# Function to get results
simulate_MVO_portfolio <- function(target_return_scale=0.8, sample_size=25, ndays=2*365, train_p=0.8, how="all"){

    data <- sample_matrix(sample_size, ndays, train_p)
    X <- data$train
    mu <- as.matrix(apply(X, 2, mean))
    vcov_matrix <- var(X)

    if(!is_solvable(vcov_matrix)){
        return(t(matrix(rep(NA, 2*4))))
    }

    if(how != "long"){
        target_return <- target_return_scale*max(abs(mu))
    } else {
        target_return <- target_return_scale*max(mu)
    }
    
    robust_estimators <- LPTN_estimators(X, mu, vcov_matrix)

    weights <- list(empirical=MVO(mu, vcov_matrix, target_return, how), 
                    robust=MVO(robust_estimators$mu, robust_estimators$vcov_matrix, target_return, how))

    returns <- list(empirical=as.vector(data$test %*% weights$empirical), 
                    robust=as.vector(data$test %*% weights$robust),
                    train=as.vector(rep(1/nrow(data$train), nrow(data$train)) %*% data$train),
                    test=as.vector(rep(1/nrow(data$test), nrow(data$test)) %*% data$test))

    output <- t(matrix(rep(NA, 2*length(returns))))
    colnames(output) <- c(paste("mean", names(returns), sep = "_"), paste("sdev", names(returns), sep = "_"))
    output[c(1:4)] <- unlist(lapply(returns, mean))
    output[-c(1:4)] <- unlist(lapply(returns, sd))

    return(output)

}

simulation_result <- simulate_MVO_portfolio(how="all")
simulation_result <- simulate_MVO_portfolio(how="long")
simulation_result <- simulate_MVO_portfolio(how="sharpe")



simulate_results <- function(n_iterations_power=200, sample_size=25, ndays=2*365, train_p=0.8, how="all"){

    print("Computing preliminary sample...", quote=FALSE)
    result_preliminary <- matrix(nrow=n_iterations_power, ncol=8)
    pb <- txtProgressBar(min = 0, max = n_iterations_power, style = 3)
    counter <- 0
    while(counter < n_iterations_power){
        simulation_result <- simulate_MVO_portfolio(how="all")
        if(any(is.na(simulation_result))){
            next
        }
        result_preliminary[counter+1, ] <- simulation_result
        counter <- counter + 1
        setTxtProgressBar(pb, counter)
    }
    colnames(result_preliminary) <- colnames(simulation_result)
    sigma_sq <- apply(result_preliminary[,c(1:2)], 2, var)
    sigma <- sqrt(sum(sigma_sq))

    target_width <- log(1+0.005)/52
    current_width <- qnorm(0.975)*sigma
    n <- ceiling((current_width/target_width)^2)
    print(paste("Sample size: ", n, ".", sep = ""), quote=FALSE)

    result_final <- matrix(nrow=n, ncol=8)
    pb <- txtProgressBar(min = 0, max = n, style = 3)
    counter <- 0
    while(counter < n){
        simulation_result <- simulate_MVO_portfolio(how="all")
        if(any(is.na(simulation_result))){
            next
        }
        result_final[counter+1, ] <- simulation_result
        counter <- counter + 1
        setTxtProgressBar(pb, counter)
        if(counter %% 100 == 0){
            print(apply(result_final[1:(counter-1), ], 2, mean))
        }
    }
    colnames(result_final) <- colnames(simulation_result)
    print(apply(result_final, 2, mean))

    return(result_final)

}

results_all <- simulate_results(how="all")
results_long <- simulate_results(how="long")
results_sharpe <- simulate_results(how="sharpe")

# Save the results
write.csv(results_all, "./Data/MVO_all.csv", row.names = FALSE)
write.csv(results_long, "./Data/MVO_long.csv", row.names = FALSE)
write.csv(results_sharpe, "./Data/MVO_sharpe.csv", row.names = FALSE)

