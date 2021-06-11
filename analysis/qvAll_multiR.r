library(tidyverse)
library(here)
if (interactive()){
  library(colorout)
}
library(R2jags)
library(tidybayes)
# library(mcmcplots)

completeness = read_csv(
  here('data', 'derived_data', 'strat_completeness.csv'),
  col_types = 'nnnnnn'
) %>% filter(dt < 40)

completeness_model = completeness %>% group_by(qv) %>% sample_frac(0.4) %>%
mutate(
  qvLabel = factor(case_when(
    qv == 1.0 ~ "zero",
    qv == 1.5 ~ "low",
    qv == 3.0 ~ "high"
  ), levels = c('zero', 'low', 'high')),
  rLabel = factor(as.integer(rLabel)),
  # combo = interaction(qvFactor, rLabel, sep = "_")
) %>% ungroup()

start = Sys.time()
print(paste('Fitting JAGS model for stratigraphic completeness. Start time:', start))

comletenessModelData = completeness_model %>% select(dt, f, rLabel, qvLabel) %>% compose_data()

write_rds(comletenessModelData, here('data', 'derived_data', 'sadler_exponent_modelData.rds'))

params_to_save = c(
  'mu.beta0', 'sigma.beta0',
  'mu.beta1', 'sigma.beta1',
  'beta0', 'beta1', 'sigma'
)

multiRSadlerJags = jags.parallel(
  model.file = here('analysis', 'qvAll_multiR_jags.r'),
  data = comletenessModelData, parameters.to.save = params_to_save,
  n.chains = 4, n.iter = 10000, n.burnin = 2500, n.thin = 25
)

end = Sys.time()
print(paste('End time:', end))
dur = difftime(end, start)
print(paste('Total time:', round(dur, 3), attr(dur, 'units')))

write_rds(multiRSadlerJags, here('data', 'derived_data', 'sadler_exponent_posteriors.rds'))
