library(sn)
library(ggplot2)
library(dplyr)

n <- 10000
p <- 100
X <- matrix(rnorm(n*p, mean=0, sd=1), ncol=p)

betas <- rep(0, p)
betas <- c()
Y <- beta_1 * X[]


Y <- matrix(beta_0 + beta*x1) 

results <- c()

sample_regression <- function (n, eps=0.5) {
  Y_eps <- Y[1:n] + rt(n, df=2)
  
  model_1 <- lm(Y_eps[1:n] ~ x1[1:n] + x2[1:n])
  model_2 <- lm(Y_eps[1:n] ~ x1[1:n] + x3[1:n])
  
  beta_1_hat <- unname(model_1$coefficients[2])
  beta_1_sd = unname(sqrt(diag(vcov(model_1)))[2])
  
  beta_2_hat <- unname(model_2$coefficients[2])
  beta_2_sd = unname(sqrt(diag(vcov(model_2)))[2])
  
  # Thin
  beta_1_train <- rnorm(1, mean=eps*beta_1_hat, sd=beta_1_sd * sqrt(eps*(1-eps)))
  beta_1_test <- beta_1_hat - beta_1_train
  
  beta_2_train <- rnorm(1, mean=eps*beta_2_hat, sd=beta_2_sd * sqrt(eps*(1-eps)))
  beta_2_test <- beta_2_hat - beta_2_train
  
  train1_sd <- beta_1_sd / sqrt(eps)
  test1_sd <- beta_1_sd / sqrt(1-eps)
  train2_sd <- beta_2_sd / sqrt(eps)
  test2_sd <- beta_2_sd / sqrt(1-eps)
  
  if ( sum(model_1$residuals^2) < sum(model_2$residuals^2) ) {
    contains_beta <- (beta <= beta_1_hat + 1.64*beta_1_sd) & (beta >= beta_1_hat - 1.64*beta_1_sd) 
  } else {
    contains_beta <- (beta <= beta_2_hat + 1.64*beta_2_sd) & (beta >= beta_2_hat - 1.64*beta_2_sd) 
  }
  
  # if ( abs(beta_1_train / train1_sd) > abs(beta_2_train / train2_sd) ) {
  #   train_pvalue <- pnorm(abs(beta_1_train) / train1_sd, lower.tail=FALSE)
  #   test_pvalue <- pnorm(abs(beta_1_test) / test1_sd, lower.tail=FALSE)
  #   
  #   contains_train <- (beta_1 <= beta_1_train / eps + 1.64*train1_sd) & (beta_1 >= beta_1_train / eps - 1.64*train1_sd) 
  #   contains_test <- (beta_1 <= beta_1_test / (1-eps) + 1.64*test1_sd) & (beta_1 >= beta_1_test / (1-eps) - 1.64*test1_sd) 
  #   
  # } else {
  #   train_pvalue <- pnorm(abs(beta_2_train) / train2_sd, lower.tail=FALSE)
  #   test_pvalue <- pnorm(abs(beta_2_test) / test2_sd, lower.tail=FALSE)
  #   
  #   contains_train <- (beta_2 <= beta_2_train / eps + 1.64*train2_sd) & (beta_2 >= beta_2_train / (eps) - 1.64*train2_sd) 
  #   contains_test <- (beta_2 <= beta_2_test / (1-eps) + 1.64*test2_sd) & (beta_2 >= beta_2_test /(1-eps) - 1.64*test2_sd) 
  # }
  beta_1_pvalue <- pnorm(abs(beta_1_hat) / beta_1_sd, lower.tail=FALSE)
  contains_beta_1 <- (beta <= beta_1_hat + 1.64*beta_1_sd) & (beta >= beta_1_hat - 1.64*beta_1_sd) 
  return(
    #c('Train', train_pvalue, 'Test', test_pvalue, 'Raw', beta_1_pvalue)
    # c('Train', contains_train, 'Test', contains_test, 'Beta1', contains_beta_1, "BestBeta", contains_beta)
    c('Beta1', contains_beta_1, "BestBeta", contains_beta)
  )
}

n_reps <- 10000

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
ggsave(paste0('../figures/lm_winner_coverage.png'), width = 4, height = 4, unit = "in")

