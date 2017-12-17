set.seed(123)

rm(list = ls())

#-------------------------------------------------#
# DATEN GENERIEREN
#-------------------------------------------------#

n <- 1000
beta <- rnorm(n = 5, mean = 0, sd = 10)

X <- matrix(data = rnorm(n * (length(beta) - 1)), nrow = n, ncol = length(beta) - 1)
X <- cbind(1, X)

y <- beta %*% t(X) + rnorm(n)

df <- as.data.frame(cbind(X[, -1], t(y)))
names(df) <- c("x1", "x2", "x3", "x4", "y")

#-------------------------------------------------#
# LINEARES MODELL RECHNEN
#-------------------------------------------------#

ols <- lm(y ~ x1 + x2 + x3 + x4, df)

beta
coef(ols)

#-------------------------------------------------#
# SIMULATION
#-------------------------------------------------#

library(MASS)
library(reshape2)
library(ggplot2)

n_sim <- 1000

beta_hat <- coef(ols)
V_hat <- vcov(ols)

S <- mvrnorm(n_sim, beta_hat, V_hat)

scenario_1 <- c(1, mean(df$x1) - sd(df$x1), mean(df$x2), mean(df$x3), mean(df$x4))
scenario_2 <- c(1, mean(df$x1) + sd(df$x1), mean(df$x2), mean(df$x3), mean(df$x4))

theta_c_1 <- apply(S, 1, function(x) x %*% as.matrix(scenario_1))
theta_c_2 <- apply(S, 1, function(x) x %*% as.matrix(scenario_2))

sigma_est <- sqrt(sum(residuals(ols)^2) / (nrow(df) - length(beta_hat)))

y_hat_1 <- apply(as.matrix(theta_c_1), 1, function(x) rnorm(1, x, sigma_est))
y_hat_2 <- apply(as.matrix(theta_c_2), 1, function(x) rnorm(1, x, sigma_est))


#-------------------------------------------------#
# BAYESIANISCHES MODELL
#-------------------------------------------------#

library(rjags)
library(coda)

model_string <- "model{

  for(i in 1:N){
  
    # Likelihood
    y[i] ~ dnorm(mu[i], sigma)
    mu[i] <-  beta[] %*% X[i, ]

   }

    # Variance priors
    sigma ~ dgamma(0.01, 0.01)
  
    # Beta Priors
    for(j in 1:J){
      beta[j] ~ dnorm(0, 0.0001)
    }
}"

data_list <- list(y = y[1, 1:n],
                  J = length(beta),
                  N = n,
                  X = X
                  )

init_list <- list(beta = rnorm(length(beta)),
                  sigma = runif(1))

jags_model <- jags.model(file = (textConnection(model_string)),
                         data = data_list,
                         inits = init_list,
                         n.chains = 2,
                         n.adapt = 1000
                         )

update(jags_model, 10000)

jags_result <- coda.samples(jags_model, 
                            variable.names = c("beta", "sigma"),
                            n.iter = 2000
                            )

summary(jags_result)

chain_1 <- jags_result[[1]]
chain_2 <- jags_result[[1]]

result <- mcmc.list(as.mcmc(chain_1), as.mcmc(chain_2))

#-------------------------------------------------#
# DIAGNOSTIC
#-------------------------------------------------#

library(ggmcmc)

S <- ggs(jags_result)
ggs_density(S)
ggs_traceplot(S)

#-------------------------------------------------#
# COMPARISION
#-------------------------------------------------#

library(ggplot2)

beta_2_hat <- as.matrix(chain_1[, 3])
plot_beta_2 <- data.frame(beta_2_hat = beta_2_hat[, 1])

ggplot(plot_beta_2, aes(x = beta_2_hat)) + geom_histogram() + 
  geom_vline(xintercept = coef(ols)[3]) + 
  theme_bw() + 
  xlab(expression(hat(beta)[2])) + 
  ylab("Anzahl")

#-------------------------------------------------#
# SIMULATION
#-------------------------------------------------#

beta_hat <- as.matrix(chain_1[, -6])

scenario_1 <- c(1, mean(df$x1) - sd(df$x1), mean(df$x2), mean(df$x3), mean(df$x4))
scenario_2 <- c(1, mean(df$x1) + sd(df$x1), mean(df$x2), mean(df$x3), mean(df$x4))

  
y_hat_1 <- beta_hat %*% scenario_1
y_hat_2 <- beta_hat %*% scenario_2

library(reshape2)

df_plot <- data.frame(y_hat_1, y_hat_2)
df_plot <- melt(df_plot)

ggplot(df_plot, aes(y = value, x = variable)) + geom_boxplot()




