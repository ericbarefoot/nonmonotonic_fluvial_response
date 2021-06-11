## model for levee-building

model {

## Priors

for (i in 1:n_qvLabel) {

  mu.beta0[i] ~ dunif(0.1, 1)
  sigma.beta0[i] ~ dunif(0.1, 5)
  tau.beta0[i] = 1 / (sigma.beta0[i] * sigma.beta0[i])

  mu.beta1[i] ~ dunif(0.1, 1)
  sigma.beta1[i] ~ dunif(0, 5)
  tau.beta1[i] = 1 / (sigma.beta1[i] * sigma.beta1[i])

  for (j in 1:n_rLabel) {
    beta0[i, j] ~ dnorm(mu.beta0[i], tau.beta0[i])T(0.01,)
    beta1[i, j] ~ dnorm(mu.beta1[i], tau.beta1[i])T(0.01,)
  }

}

sigma ~ dunif(0,10)

tau = 1 / (sigma * sigma)

## Model

for (i in 1:n) {
  mu[i] = beta0[qvLabel[i], rLabel[i]] * dt[i] ^ beta1[qvLabel[i], rLabel[i]]
  f[i] ~ dnorm(mu[i], tau)
}

# for (i in 1:n_qvLabel) {
#   mu.interOne[i] = (1 / mu.beta0[i]) ^ (1 / mu.beta1[i])
#   for (j in 1:n_rLabel) {
#     interOne[i, j] = (1 / beta0[i, j]) ^ (1 / beta1[i, j])
#   }
# }

}
