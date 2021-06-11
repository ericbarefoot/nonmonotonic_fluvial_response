#------------------------------------------------------------------------------#
# This model estimates an anova on the
# abundance of levee breaches per leg of experiment
#------------------------------------------------------------------------------#

model {

## Priors

# for breach model

for (i in 1:n_qv) {

  #------------------------------------------------------------------------------#
  # use time and angle as a random effect -- Change Likelihood
  #------------------------------------------------------------------------------#
  # for (aa in 1:n_angle) {
  #   for (tt in 1:n_time) {
  #     breach.mod[i, aa, tt] ~ dnorm(mu.breach[i], tau.breach[i])
  #   }
  # }
  #------------------------------------------------------------------------------#

  #------------------------------------------------------------------------------#
  # use angle as a random effect -- Change Likelihood
  #------------------------------------------------------------------------------#
  # for (aa in 1:n_angle) {
  #   breach.mod[i, aa] ~ dnorm(mu.breach[i], tau.breach[i])
  #   breach.int[i, aa] ~ dnorm(mu.breach[i], tau.breach[i])
  #   # breach.cnt[i, aa] = exp(breach.int[i, aa])
  # }
  #------------------------------------------------------------------------------#

  #------------------------------------------------------------------------------#
  # use time as a random effect -- Change Likelihood
  #------------------------------------------------------------------------------#
  # for (tt in 1:n_time) {
  #   breach.mod[i, tt] ~ dnorm(mu.breach[i], tau.breach[i])
  # }
  #------------------------------------------------------------------------------#

  #------------------------------------------------------------------------------#
  # use no random effect -- Change Likelihood
  #------------------------------------------------------------------------------#
  breach.int[i] ~ dnorm(0, 0.001)
  breach.cnt[i] = exp(breach.int[i])

  # breach.slp[i] ~ dnorm(0, 0.001)
  # breach.slp.dim[i] = exp(breach.slp[i])

  #------------------------------------------------------------------------------#

  # mu.breach[i] ~ dnorm(0, 0.001)
  # sd.breach[i] ~ dunif(0,15)
  # tau.breach[i] = 1 / (sd.breach[i] * sd.breach[i])
  #
  # mu.breach.cnt[i] = exp(mu.breach[i])
}

# sd.breach.resid = 1 / sqrt(tau.breach.resid)
# tau.breach.resid ~ dgamma(0.001, 0.001)

## Model

for (i in 1:n) {

  lambda[i] = exp(breach.int[qv[i]])
  # lambda[i] = exp(breach.int[qv[i]] + breach.slp[qv[i]] * aspect[i])
  count[i] ~ dpois(lambda[i])

}

}
