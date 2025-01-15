library(sn)
library(ggplot2)
library(dplyr)
library(MASS)

n <- 10000
p <- 2
alpha <- 0.1

x1 <- rnorm(n, mean = 0, sd = 1)
x2 <- rnorm(n, mean = 0, sd = 1)
X <- cbind(x1, x2)

beta_0 <- 0.5
beta_1 <- 0.1
beta_2 <- 0.1
Y <- matrix(beta_0 + beta_1*x1 + beta_2*x2) 

results <- c()

sample_regression <- function (n, eps=0.5) {
  Y_eps <- Y[1:n] + rt(n, df=2)
  
  fit <- lm(Y_eps[1:n] ~ x1[1:n] + x2[1:n])
  res <- fit$residuals
  
  beta_0_hat <- unname(fit$coefficients[1])
  
  beta_1_hat <- unname(fit$coefficients[2])
  beta_1_sd = unname(sqrt(diag(vcov(fit)))[2])

  beta_2_hat <- unname(fit$coefficients[3])
  beta_2_sd = unname(sqrt(diag(vcov(fit)))[3])
  
  # Thin
  VCOV <- 1/(n-p-1) * as.numeric(t(res)%*%res) * solve(t(X)%*%X)
  noise <- mvrnorm(1, mu=rep(0, ncol(X)), Sigma=VCOV * eps*(1-eps))
  
  beta_1_train <- beta_1_hat + noise[1] / eps
  beta_1_test <- beta_1_hat - noise[1] / (1 - eps)
  
  beta_2_train <- beta_2_hat + noise[2] / eps
  beta_2_test <- beta_2_hat - noise[2] / (1 - eps)
  
  res_train <- Y_eps - beta_0 - beta_1_train*x1 - beta_2_train*x2
  VCOV_train <- 1/(n-p-1) * as.numeric(t(res_train)%*%res_train) * solve(t(X)%*%X)
  
  res_test <- Y_eps - beta_0 - beta_1_test*x1 - beta_2_test*x2
  VCOV_test <- 1/(n-p-1) * as.numeric(t(res_test)%*%res_test) * solve(t(X)%*%X)
  
  train1_sd <- sqrt(diag(VCOV_train)[1])
  train2_sd <- sqrt(diag(VCOV_train)[2])
  
  test1_sd <- sqrt(diag(VCOV_test)[1])
  test2_sd <- sqrt(diag(VCOV_test)[2])
  
  # train1_sd <- beta_1_sd / sqrt(eps)
  # test1_sd <- beta_1_sd / sqrt(1-eps)
  # train2_sd <- beta_2_sd / sqrt(eps)
  # test2_sd <- beta_2_sd / sqrt(1-eps)
  
  zscore <- qnorm(1 - alpha/2)
  
  # Post selection (uncorrected) inference
  if ( abs(beta_1_hat / beta_1_sd) > abs(beta_2_hat / beta_2_sd) ) {
    contains_beta <- (beta_1 <= beta_1_hat + zscore*beta_1_sd) & (beta_1 >= beta_1_hat - zscore*beta_1_sd)
    raw_bias <- beta_1 - beta_1_hat
  } else {
    contains_beta <- (beta_2 <= beta_2_hat + zscore*beta_2_sd) & (beta_2 >= beta_2_hat - zscore*beta_2_sd) 
    raw_bias <- beta_2 - beta_2_hat
  }
  
  # Post selection (thinning) inference
  if ( abs(beta_1_train / train1_sd) > abs(beta_2_train / train2_sd) ) {
    train_pvalue <- pnorm(abs(beta_1_train) / train1_sd, lower.tail=FALSE)
    test_pvalue <- pnorm(abs(beta_1_test) / test1_sd, lower.tail=FALSE)
    
    contains_train <- (beta_1 <= beta_1_train / eps + zscore*train1_sd) & (beta_1 >= beta_1_train / eps - zscore*train1_sd) 
    contains_test <- (beta_1 <= beta_1_test / (1-eps) + zscore*test1_sd) & (beta_1 >= beta_1_test / (1-eps) - zscore*test1_sd) 
    
    train_bias <- beta_1 - beta_1_train
    test_bias <- beta_1 - beta_1_test
    
  } else {
    train_pvalue <- pnorm(abs(beta_2_train) / train2_sd, lower.tail=FALSE)
    test_pvalue <- pnorm(abs(beta_2_test) / test2_sd, lower.tail=FALSE)
    
    contains_train <- (beta_2 <= beta_2_train / eps + zscore*train2_sd) & (beta_2 >= beta_2_train / (eps) - zscore*train2_sd) 
    contains_test <- (beta_2 <= beta_2_test / (1-eps) + zscore*test2_sd) & (beta_2 >= beta_2_test /(1-eps) - zscore*test2_sd) 
    
    train_bias <- beta_2 - beta_2_train
    test_bias <- beta_2 - beta_2_test
  }
  
  # Check if contains a fixed beta
  beta_1_pvalue <- pnorm(abs(beta_1_hat) / beta_1_sd, lower.tail=FALSE)
  contains_beta_1 <- (beta_1 <= beta_1_hat + zscore*beta_1_sd) & (beta_1 >= beta_1_hat - zscore*beta_1_sd)
  beta_1_bias <- beta_1 - beta_1_hat
  
  # Returns answers
  return(
    # c('Train', train_pvalue, 'Test', test_pvalue, 'Raw', beta_1_pvalue)
    c('Train', contains_train, 'Test', contains_test, 'Beta1', contains_beta_1, "BestBeta", contains_beta)
  )
}

n_reps <- 1000

eps <- 0.5

sample_sizes <- c(20, 50, 100, 1000)
results_list = vector("list", length = length(sample_sizes))

for (i in seq(length(sample_sizes))) {
  n <- sample_sizes[i]
  results <- replicate(n_reps, sample_regression(n, eps))
  
  results_df <- data.frame(
    data=t(array(results, dim = c(2, 3*n_reps)))
  )
  colnames(results_df) <- c("Set", "contains")
  results_df$contains <- as.integer(as.logical(results_df$contains))
  results_df$n <- n
  
  results_list[[i]] <- results_df
}

results_df <- do.call(rbind, results_list)
results_df$n <- as.numeric(results_df$n)

plot_df <- results_df %>%
  group_by(Set, n) %>%
  summarize(
    coverage = mean(contains)
    ) 
# %>%
#   mutate(
#     sd = sqrt(coverage * (1 - coverage) / n),
#     ci_upper = coverage + 1.96*sd,
#     ci_lower = coverage - 1.96*sd
#   )

g <- ggplot(
  data=plot_df,
  aes(x=n, y=coverage, col=Set, group=Set))+
  geom_line(aes(linetype = Set, color = Set, group = Set), size = 1) +
  geom_point(aes(color = Set, group=Set), size = 1.5) +
  # geom_errorbar(aes(ymin=coverage-1.96*se, ymax=coverage+1.96*se), width=.1) +
  geom_hline(yintercept=0.9, linetype="dashed", color = "gray", size=0.5) +
  xlab('Sample Size') +
  ylab('Coverage')
print(g)
# ggsave(paste0('../figures/lm_winner_coverage.png'), width = 4, height = 4, unit = "in")

