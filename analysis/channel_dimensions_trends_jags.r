#------------------------------------------------------------------------------#
# This model treats radius as a random effect and duratoin as a fixed category
#------------------------------------------------------------------------------#

model {

## Priors

# for slope model

for (i in 1:n_qv) {

  #------------------------------------------------------------------------------#
  # use time and radius as a random effect -- Change Likelihood
  #------------------------------------------------------------------------------#
  for (rr in 1:n_radius) {
    for (tt in 1:n_time) {

      depth.mod[i, rr, tt] ~ dnorm(mu.depth[i], tau.depth[i])
      width.mod[i, rr, tt] ~ dnorm(mu.width[i], tau.width[i])

    }
  }
  #------------------------------------------------------------------------------#

  #------------------------------------------------------------------------------#
  # use radius as a random effect -- Change Likelihood
  #------------------------------------------------------------------------------#
  # for (rr in 1:n_radius) {
  #     depth.mod[i, rr] ~ dnorm(mu.depth[i], tau.depth[i])
  #     width.mod[i, rr] ~ dnorm(mu.width[i], tau.width[i])
  # }
  #------------------------------------------------------------------------------#

  #------------------------------------------------------------------------------#
  # use time as a random effect -- Change Likelihood
  #------------------------------------------------------------------------------#
  # for (tt in 1:n_time) {
  #     depth.mod[i, tt] ~ dnorm(mu.depth[i], tau.depth[i])
  #     width.mod[i, tt] ~ dnorm(mu.width[i], tau.width[i])
  # }
  #------------------------------------------------------------------------------#

  #------------------------------------------------------------------------------#
  # use no random effects -- Change Likelihood
  #------------------------------------------------------------------------------#
  # depth.mod[i] ~ dnorm(0, 0.001)
  # width.mod[i] ~ dnorm(0, 0.001)
  #------------------------------------------------------------------------------#

  mu.depth[i] ~ dnorm(0, 0.001)
  sd.depth[i] ~ dunif(0,15)
  tau.depth[i] = 1 / (sd.depth[i] * sd.depth[i])

  mu.width[i] ~ dnorm(0, 0.001)
  sd.width[i] ~ dunif(0,15)
  tau.width[i] = 1 / (sd.width[i] * sd.width[i])

}

sd.depth.resid = 1 / sqrt(tau.depth.resid)
tau.depth.resid ~ dgamma(0.001, 0.001)

sd.width.resid = 1 / sqrt(tau.width.resid)
tau.width.resid ~ dgamma(0.001, 0.001)

## Model

for (i in 1:n) {

  H[i] = depth.mod[qv[i], radius[i], time[i]]
  depth[i] ~ dnorm(H[i], tau.depth.resid)

  B[i] = width.mod[qv[i], radius[i], time[i]]
  width[i] ~ dnorm(B[i], tau.width.resid)

}

}
