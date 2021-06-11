#------------------------------------------------------------------------------#
# This model treats radius as a random effect and duratoin as a fixed category
#------------------------------------------------------------------------------#

model {

## Priors

# for slope model

for (i in 1:n_qv) {
  #------------------------------------------------------------------------------#
  # use time and angle as a random effect -- Change Likelihood
  #------------------------------------------------------------------------------#
  for (aa in 1:n_angle) {
    for (tt in 1:n_time) {
      slope.mod[i, aa, tt] ~ dnorm(mu.slope[i], tau.slope[i])
    }
  }
  #------------------------------------------------------------------------------#

  #------------------------------------------------------------------------------#
  # use angle as a random effect -- Change Likelihood
  #------------------------------------------------------------------------------#
  # for (aa in 1:n_angle) {
  #   slope.mod[i, aa] ~ dnorm(mu.slope[i], tau.slope[i])
  # }
  #------------------------------------------------------------------------------#

  #------------------------------------------------------------------------------#
  # use time as a random effect -- Change Likelihood
  #------------------------------------------------------------------------------#
  # for (tt in 1:n_time) {
  #   slope.mod[i, tt] ~ dnorm(mu.slope[i], tau.slope[i])
  # }
  #------------------------------------------------------------------------------#

  #------------------------------------------------------------------------------#
  # use no random effect -- Change Likelihood
  #------------------------------------------------------------------------------#
  # slope.mod[i] ~ dnorm(0, 0.001)
  #------------------------------------------------------------------------------#

  mu.slope[i] ~ dnorm(0, 0.001)
  sd.slope[i] ~ dunif(0,15)
  tau.slope[i] = 1 / (sd.slope[i] * sd.slope[i])
}

sd.slope.resid = 1 / sqrt(tau.slope.resid)
tau.slope.resid ~ dgamma(0.001, 0.001)

## Model

for (i in 1:n) {

  S[i] = slope.mod[qv[i], angle[i], time[i]]
  slope[i] ~ dnorm(S[i], tau.slope.resid)

}

}
