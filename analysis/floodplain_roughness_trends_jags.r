#------------------------------------------------------------------------------#
# This model treats radius as a random effect and duratoin as a fixed category
#------------------------------------------------------------------------------#

model {

## Priors

# for roughness model

for (i in 1:n_qv) {

  #------------------------------------------------------------------------------#
  # use time and angle as a random effect -- Change Likelihood
  #------------------------------------------------------------------------------#
  # for (aa in 1:n_angle) {
  #   for (tt in 1:n_time) {
  #     roughness.mod[i, aa, tt] ~ dnorm(mu.roughness[i], tau.roughness[i])
  #   }
  # }
  #------------------------------------------------------------------------------#

  #------------------------------------------------------------------------------#
  # use angle as a random effect -- Change Likelihood
  #------------------------------------------------------------------------------#
  # for (aa in 1:n_angle) {
  #   roughness.mod[i, aa] ~ dnorm(mu.roughness[i], tau.roughness[i])
  # }
  #------------------------------------------------------------------------------#

  #------------------------------------------------------------------------------#
  # use time as a random effect -- Change Likelihood
  #------------------------------------------------------------------------------#
  # for (tt in 1:n_time) {
  #   roughness.mod[i, tt] ~ dnorm(mu.roughness[i], tau.roughness[i])
  # }
  #------------------------------------------------------------------------------#

  #------------------------------------------------------------------------------#
  # use no random effect -- Change Likelihood
  #------------------------------------------------------------------------------#
  roughness.mod[i] ~ dnorm(0, 0.001)
  #------------------------------------------------------------------------------#

  # mu.roughness[i] ~ dnorm(0, 0.001)
  # sd.roughness[i] ~ dunif(0,15)
  # tau.roughness[i] = 1 / (sd.roughness[i] * sd.roughness[i])
}

sd.roughness.resid = 1 / sqrt(tau.roughness.resid)
tau.roughness.resid ~ dgamma(0.001, 0.001)

## Model

for (i in 1:n) {

  ruff[i] = roughness.mod[qv[i]]
  short_wave_relief[i] ~ dnorm(ruff[i], tau.roughness.resid)

}

}
