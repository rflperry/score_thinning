library(sn)
library(ggplot2)
library(dplyr)
library(MASS)
library(glmnet)

N_TOTAL <- 10000
p <- 25
alpha <- 0.05

Z <- matrix(rnorm(p*p, mean = 0, sd = 1), ncol=p)
covar <- t(Z) %*% Z

beta <- c() # c(0, seq(1, 0, by=-0.2))
beta <- c(beta, rep(0, p+1-length(beta)))
X_full <- cbind(
  rep(1, N_TOTAL),
  matrix(rnorm(N_TOTAL*p, mean = 0, sd = 1), ncol=p)
)

Y_noiseless <- X_full %*% beta

sample_regression <- function (n, c=1) {
  
  X <- X_full[1:n,]
  
  # All distributions (except cauchy) have variance 3
  Y <- Y_noiseless[1:n] + rnorm(n, mean=0, sd=1)# rt(n, df=3)
  Y_val <- Y_noiseless[1:n] + rnorm(n, mean=0, sd=1)# rt(n, df=3)
  
  # Conventional data thinning
  sd_Y <- 1
  noise <- rnorm(n, mean=0, sd=1)
  Y_train_thin <- Y + noise * sd_Y * noise * c
  Y_test_thin <- Y - noise * sd_Y * noise / c
  Y_val_thin <- Y_val - noise * sd_Y * noise / c
  
  # Sample splitting
  indices <- sample(1:n, replace=FALSE)
  Y_train_split <- Y[indices[1:(n/2)]]
  Y_test_split <- Y[indices[(n/2 + 1):n]]
  X_train_split <- X[indices[1:(n/2)]]
  X_test_split <- X[indices[(n/2 + 1):n]]

  # Fit models
  fit_split <- lm(Y_train_split ~ 0 + X_train_split)  # intercept is baked into X already
  fit_thin <- lm(Y_train_thin ~ 0 + X) 
  
  # Compute errors
  fit_split_test_error <- mean((Y_test_split - predict(fit_split, as.data.frame(X_test_split)))^2)
  fit_split_val_error <- mean((Y_val - predict(fit_split, as.data.frame(X)))^2)
  fit_thin_test_error <- mean((Y_test_thin - predict(fit_thin, as.data.frame(X)))^2)
  # fit_thin_val_error <- mean((Y_val - predict(fit_thin, as.data.frame(X)))^2)
  fit_thin_val_error <- mean((Y_val_thin - predict(fit_thin, as.data.frame(X)))^2)
  
  # Returns answers
  return(
    c(
      'split_test_error', n, c, fit_split_test_error,
      'split_val_error', n, c, fit_split_val_error,
      'thin_test_error', n, c, fit_thin_test_error,
      'thin_val_error', n, c, fit_thin_val_error
    )
  )
  
  # # Setting the range of lambda values
  # lambda_seq <- 10^seq(2, -2, by = -.25)
  # 
  # # Regular, sample splitting. 2 folds equivalent to our thinning
  # ridge_cv <- cv.glmnet(X, Y, alpha = 0, lambda  = lambda_seq, intercept=FALSE, nfolds = 2)
  # best_lambda_cv <- ridge_cv$lambda.min
  # best_ridge_cv <- glmnet(X, Y, alpha = 0, lambda  = best_lambda_cv, intercept=FALSE)
  # 
  # # Conventional data thinning
  # sd_Y <- 1
  # noise <- rnorm(n, mean=0, sd=1)
  # Y_train <- Y + noise * sd_Y * noise * c
  # Y_test <- Y - noise * sd_Y * noise / c
  # 
  # # Thinning
  # best_ridge_thin <- NA
  # best_mse_thin <- NA
  # for (lambda in lambda_seq) {
  #   ridge_thin <- glmnet(X, Y_train, alpha = 0, lambda  = lambda, intercept=FALSE)
  #   mse <- mean((Y_test - predict(ridge_thin, as.data.frame(X)))^2)
  #   if (is.na(best_mse_thin) | (best_mse_thin > mse)) {
  #     best_mse_thin <- mse
  #     best_ridge_thin <- ridge_thin
  #   }
  # }
  # 
  # # Compute errors
  # cv_val_error <- mean((Y_val - predict(best_ridge_cv, as.data.frame(X)))^2)
  # cv_test_error <- mean((Y_test - predict(best_ridge_cv, as.data.frame(X)))^2)
  # thin_val_error <- mean((Y_val - predict(best_ridge_thin, as.data.frame(X)))^2)
  # thin_test_error <- mean((Y_test - predict(best_ridge_thin, as.data.frame(X)))^2)
}

n_reps <- 100
sample_sizes <- c(1000)
results_list = vector("list", length = length(sample_sizes))

header <- c("Method", "n", "c", "MSE")
for (i in seq(length(sample_sizes))) {
    n <- sample_sizes[i]
    results <- replicate(n_reps, sample_regression(
      n=n, c=1
    ))
    
    results_df <- data.frame(
      data=t(array(results, dim = c(length(header), as.numeric(length(results) / length(header)))))
    )
    results_list[[i]] <- results_df
}

results_df <- do.call(rbind, results_list)
colnames(results_df) <- header
results_df[,2:length(header)] <- sapply(results_df[,2:length(header)], as.numeric)
tail(results_df, n=8)

results_df %>% group_by(Method) %>% summarize(mean_mse = mean(MSE))

plot_df <- results_df %>%
  mutate(
    covered = as.numeric(
      (estimand >= estimate - crit_val * estimate_sd) &
        (estimand <= estimate + crit_val * estimate_sd)
    )
  ) %>%
  subset( Method %in% c("Naive", "Beta1", "Test_DT", "Test_ST") ) %>%
  group_by(Method, n, distribution, var_model) %>%
  summarize(
    coverage = mean(covered),
    mean_bias = mean(estimate - estimand)
  ) 

g <- ggplot(
  data=plot_df,
  aes(x=n, y=coverage, col=Method, group=Method))+
  geom_line(aes(linetype = Method, color = Method, group = Method), size = 1) +
  geom_point(aes(color = Method, group=Method), size = 1.5) +
  # geom_errorbar(aes(ymin=coverage-1.96*se, ymax=coverage+1.96*se), width=.1) +
  geom_hline(yintercept=0.9, linetype="dashed", color = "gray", size=0.5) +
  xlab('Sample Size')  +
  facet_grid(var_model~ distribution)
# ylab('Average')
print(g)
ggsave(paste0('../figures/lm_winner_coverage_null_varXdist_exp.png'), width = 8, height = 8, unit = "in")

