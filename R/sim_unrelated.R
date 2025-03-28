# simulation study illustrating that the Fisher transformation over-estimates pooled R^2 values
# with linear regression
P <- 10
N <- 100
V <- 0.5 + 0.5 * diag(P)

sim <- function(N, P, rho) {
  V <- rho + (1-rho) * diag(P)
  X <- matrix(rnorm(N*P), N, P) %*% chol(V)
  y <- rnorm(N)
  
  df <- data.frame(y, X)
  
  amp <- mice::ampute(df)$amp
  imp <- mice::mice(amp, method = "norm", print = FALSE)
  
  fit <- with(imp, lm(y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10))
  
  c(obs = summary(lm(y~.,df))$r.squared,
    imp = mice::pool.r.squared(fit)[,1])
}

library(pbapply)

cl <- parallel::makeCluster(4)
parallel::clusterExport(cl, c("sim", "P", "N"))

out <- pbreplicate(1000, sim(N, P, 0.5), cl = cl)

parallel::stopCluster(cl)

rowMeans(out, na.rm = TRUE)
#>        obs        imp 
#> 0.09941488 0.11977473


# simulation study illustrating that the Fisher transformation over-estimates pooled R^2 values
# with logistic regression
P <- 10
N <- 100
rho = 0.5

sim <- function(N, P, rho) {
  V <- rho + (1 - rho) * diag(P)
  X <- matrix(rnorm(N*P), N, P) %*% chol(V)
  y <- rbinom(N, size = 1, prob = 0.5)
  
  df <- data.frame(y, X)
  r2_obs <- glm(y ~ ., data = df, family = "binomial") |> nagelkerke_r2()
  
  amp <- mice::ampute(df)$amp
  imp <- mice::mice(amp, method = c("pmm", rep("norm", P)), print = FALSE)
  
  fit <- with(imp, glm(y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, family = "binomial"))
  r2_imp <- pool_nagelkerke_r2(fit)[,1]
  
  return(c(obs = r2_obs, imp = r2_imp))
}

library(pbapply)

cl <- parallel::makeCluster(4)
parallel::clusterExport(cl, c("sim", "P", "N", "nagelkerke_r2", "pool_nagelkerke_r2", "pool_r2"))

out <- pbreplicate(1000, sim(N, P, 0.5), cl = cl)

parallel::stopCluster(cl)

save(out, file = "results/sim_logistic.RData")

rowMeans(out, na.rm = TRUE)
#>        obs        imp
#>        0.1341351 0.1588817