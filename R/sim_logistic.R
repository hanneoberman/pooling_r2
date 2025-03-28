# simulation study illustrating that the Fisher transformation over-estimates pooled R^2 values
# with logistic regression

P <- 10
N <- 100
rho = 0.5
beta = 0.2

sim <- function(N, P, rho, beta) {
  V <- rho + (1 - rho) * diag(P)
  X <- matrix(rnorm(N*P), N, P) %*% chol(V)
  lin_pred <- X %*% rep(beta, P) + rnorm(N, sd = 0.5)
  p <- 1 / (1 + exp(-lin_pred))
  y <- rbinom(N, size = 1, prob = p)
  
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
parallel::clusterExport(cl, c("sim", "P", "N", "rho", "beta", "nagelkerke_r2", "pool_nagelkerke_r2", "pool_r2"))

out <- pbreplicate(1000, sim(N, P, rho, beta), cl = cl)

parallel::stopCluster(cl)

save(out, file = "results/out.RData")

rowMeans(out, na.rm = TRUE)
