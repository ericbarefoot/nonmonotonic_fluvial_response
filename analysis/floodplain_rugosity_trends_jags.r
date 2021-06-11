#------------------------------------------------------------------------------#
# This model treats radius as a random effect and duratoin as a fixed category
#------------------------------------------------------------------------------#

model {

## Priors

# for rugosity model

for (i in 1:n_qv) {

  #------------------------------------------------------------------------------#
  # use time and angle as a random effect -- Change Likelihood
  #------------------------------------------------------------------------------#
  # for (aa in 1:n_angle) {
  #   for (tt in 1:n_time) {
  #     rugosity.mod[i, aa, tt] ~ dnorm(mu.rugosity[i], tau.rugosity[i])
  #   }
  # }
  #------------------------------------------------------------------------------#

  #------------------------------------------------------------------------------#
  # use angle as a random effect -- Change Likelihood
  #------------------------------------------------------------------------------#
  # for (aa in 1:n_angle) {
  #   rugosity.mod[i, aa] ~ dnorm(mu.rugosity[i], tau.rugosity[i])
  # }
  #------------------------------------------------------------------------------#

  #------------------------------------------------------------------------------#
  # use time as a random effect -- Change Likelihood
  #------------------------------------------------------------------------------#
  # for (tt in 1:n_time) {
  #   rugosity.mod[i, tt] ~ dnorm(mu.rugosity[i], tau.rugosity[i])
  # }
  #------------------------------------------------------------------------------#

  #------------------------------------------------------------------------------#
  # use no random effect -- Change Likelihood
  #------------------------------------------------------------------------------#
  rugosity.mod[i] ~ dnorm(0, 0.001)
  #------------------------------------------------------------------------------#

  # mu.rugosity[i] ~ dnorm(0, 0.001)
  # sd.rugosity[i] ~ dunif(0,15)
  # tau.rugosity[i] = 1 / (sd.rugosity[i] * sd.rugosity[i])
}

sd.rugosity.resid = 1 / sqrt(tau.rugosity.resid)
tau.rugosity.resid ~ dgamma(0.001, 0.001)

## Model

for (i in 1:n) {

  S[i] = rugosity.mod[qv[i]]
  rugosity[i] ~ dnorm(S[i], tau.rugosity.resid)

}

}
