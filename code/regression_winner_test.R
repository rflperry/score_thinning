library(sn)
library(ggplot2)
library(dplyr)
library(MASS)

N_TOTAL <- 10000
p <- 25
alpha <- 0.1

Z <- matrix(rnorm(p*p, mean = 0, sd = 1), ncol=p)
covar <- t(Z) %*% Z
duv <- svd(covar)
ROT <- duv$u

beta <- c() # c(0, seq(1, 0, by=-0.2))
beta <- c(beta, rep(0, p+1-length(beta)))
X_full <- cbind(
  rep(1, N_TOTAL),
  matrix(rnorm(N_TOTAL*p, mean = 0, sd = 1), ncol=p) %*% ROT
  )
Y <- matrix(X_full %*% beta) 

selection <- function(beta, beta_sd) {
  return(
    which.max( abs(beta[2:length(beta)]) / beta_sd[2:length(beta)]) + 1
  )
}

X <- X_full# [,1:(p/2+1)]

sample_regression <- function (n, dist="t", var_model="homogenous", eps=0.5) {
  
  params <- c(dist, var_model, n, eps)
  
  # All distributions (except cauchy) have variance 3
  if (dist == "t") {
    Y_noise <- rt(n, df=3)
  } else if (dist == 'normal'){
    Y_noise <- rnorm(n, mean=0, sd=sqrt(3))
  } else if (dist == 'uniform'){
    Y_noise <- runif(n, min=0, max=sqrt(36))
  } else if (dist == 'binomial'){
    Y_noise <- rbinom(n, 1, prob=0.5) * 3/0.25
  } else if (dist == 'poisson'){
    Y_noise <- rpois(n, lambda=3)
  } else if (dist == 'exponential'){
    Y_noise <- rexp(n, sqrt(1/3))
  } else if (dist == 'cauchy'){
    Y_noise <- rcauchy(n)
  } else {
    stop(paste("Distribution", dist, "is not valid." ))
  }
  
  if (var_model == 'homogenous'){
    Y_eps <- Y[1:n] + Y_noise
  } else if (var_model == 'linear'){
    # Note: variance scales with square of multiplier :)
    Y_eps <- Y[1:n] + Y_noise*sqrt(abs(X[1:n,1]))
  } else if (var_model == 'quadratic'){
    Y_eps <- Y[1:n] + Y_noise*abs(X[1:n,1])
  } else {
    stop(paste("Variance model", var_model, "is not valid." ))
  }
  
  # Conventional data thinning
  # sd_Y <- sd(Y_eps)
  # noise <- rnorm(n, mean=0, sd=sqrt(eps*(1-eps)))
  # Y_eps_train <- Y_eps + noise * sd_Y / eps
  # Y_eps_test <- Y_eps - noise * sd_Y / (1-eps)
  Y_eps_train <- Y_eps * rbeta(n, eps, (1-eps))
  Y_eps_test <- Y_eps - Y_eps_train
  
  fit_train <- lm(Y_eps_train ~ 0 + X[1:n,])  # intercept is baked into X already
  fit_test <- lm(Y_eps_test ~ 0 + X[1:n,])
  
  Si_dt <- selection(
    unname(fit_train$coefficients),
    unname(sqrt(diag(vcov(fit_train))))
    )
  
  confint_dt_train <- unname(confint(fit_train, level=1-alpha)[Si_dt,])
  beta_dt_train <- unname(fit_train$coefficients)[Si_dt]
  beta_dt_train_sd = unname(sqrt(diag(vcov(fit_train))))[Si_dt]

  confint_dt_test <- unname(confint(fit_test, level=1-alpha)[Si_dt,])
  beta_dt_test <- unname(fit_test$coefficients)[Si_dt]
  beta_dt_test_sd = unname(sqrt(diag(vcov(fit_test))))[Si_dt]

  # Other approaches
  fit <- lm(Y_eps ~ 0 + X[1:n,])
  res <- fit$residuals
  
  beta_hat <- unname(fit$coefficients)
  beta_hat_sd = unname(sqrt(diag(vcov(fit))))
  
  # Thin
  VCOV <- 1/(n-p-1) * as.numeric(t(res)%*%res) * solve(t(X[1:n,])%*%X[1:n,])
  noise <- mvrnorm(1, mu=rep(0, ncol(X)), Sigma=VCOV * eps*(1-eps))
  noise[1] <- 0
  
  beta_train <- beta_hat + noise / eps
  beta_test <- beta_hat - noise / (1-eps)
  
  beta_train_sd <- sqrt(diag(VCOV)) / sqrt(eps)
  beta_test_sd <- sqrt(diag(VCOV)) / sqrt(1-eps)
  
  # Post selection (uncorrected) inference
  Si_raw <- selection(beta_train, beta_train_sd)
  raw_bias <- beta[Si_raw] - beta_hat[Si_raw]
  post_select_beta <- beta_hat[Si_raw]
  
  # Post selection (thinning) inference
  Si <- selection(beta_train, beta_train_sd)
  
  train_pvalue <- pnorm(abs(beta_train[Si]) / beta_train_sd[Si], lower.tail=FALSE)
  test_pvalue <- pnorm(abs(beta_test[Si]) / beta_test_sd[Si], lower.tail=FALSE)
  
  train_bias <- beta[Si] - beta_train[Si]
  test_bias <- beta[Si] - beta_test[Si]
  
  train_beta <- beta_train[Si]
  test_beta <- beta_test[Si]

  # Check if estimate a fixed beta
  beta_2_pvalue <- pnorm(abs(beta_hat[2]) / beta_hat_sd[2], lower.tail=FALSE)
  beta_2_bias <- beta[2] - beta_hat[2]
  
  # Returns answers
  return(
    c(
      'Train_ST', params, train_beta, beta_train_sd[Si], beta[Si],
      'Test_ST', params, test_beta, beta_test_sd[Si], beta[Si],
      'Train_DT', params, beta_dt_train, beta_dt_train_sd, beta[Si_dt],
      "Test_DT", params, beta_dt_test, beta_dt_test_sd, beta[Si_dt],
      'Beta1', params, beta_hat[2], beta_hat_sd[2], beta[2],
      "Naive", params, post_select_beta, beta_hat_sd[Si_raw], beta[Si_raw]
    )
  )
}

n_reps <- 1000

eps <- 0.5

sample_sizes <- c(60, 100, 1000)
dists <- c("t", "normal", "poisson", "exponential", "uniform", "binomial")#, "cauchy")
var_models <- c("homogenous")#, "linear", "quadratic")
results_list = vector("list", length = length(sample_sizes) * length(dists))

header <- c("Method", "distribution", "var_model", "n", "eps", "estimate", "estimate_sd", "estimand")

j <- 0
for (i in seq(length(sample_sizes))) {
  for (dist in dists) {
    for (var_model in var_models) {
      n <- sample_sizes[i]
      print(paste(dist, n, var_model))
      results <- replicate(n_reps,sample_regression(
        n=n, dist=dist, eps=eps, var_model = var_model
        ))
      
      results_df <- data.frame(
        data=t(array(results, dim = c(length(header), as.numeric(length(results) / length(header)))))
      )
      j <- j + 1
      results_list[[j]] <- results_df
    }
  }
}

results_df <- do.call(rbind, results_list)
colnames(results_df) <- header
results_df[,4:length(header)] <- sapply(results_df[,4:length(header)], as.numeric)
head(results_df)

crit_val <- qt(1 - alpha/2, df=n-length(beta))

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

