shrinkage_estimator <- function(b, s) {
  z <- b / s
  p <- c(0.32, 0.31, 0.3, 0.07)
  sigma <- c(0.61, 1.42, 2.16, 5.64)
  sigma2 <- sigma^2
  # mixing probabilities
  q <- p * dnorm(z, 0, sqrt(sigma2 + 1))
  q <- q / sum(q)  # Normalize to sum to 1
  # conditional means and variances
  pm <- b * sigma2 / (sigma2 + 1)
  pv <- s^2 * sigma2 / (sigma2 + 1)
  # shrinkage estimator (β̂)
  betahat <- sum(q * pm)
  # second moment and total variance for the mixture
  pm2 <- sum(q * pm^2)
  ps2 <- sum(q * pv)
  # standard deviation (σ) of the mixture distribution
  sigma_hat <- sqrt(ps2 + pm2 - betahat^2)
  # 95% confidence interval around β̂
  ci_lower <- betahat - 1.96 * sigma_hat
  ci_upper <- betahat + 1.96 * sigma_hat
  # Return results as a list
  list(
    betahat = betahat,
    sigma_hat = sigma_hat,
    confidence_interval = c(ci_lower, ci_upper),
    mixture_details = data.frame(q, pm, pv)
  )
}

# Example usage: result <- shrinkage_estimator(b = 0.4, s = 0.3)
# cat("Shrinkage estimator (Beta-hat):", result$betahat, "\n")
# cat("Standard deviation of the estimator (Sigma-hat):", result$sigma_hat, "\n")
# cat("95% Confidence Interval: [", result$confidence_interval[1], ", ", result$confidence_interval[2], "]\
