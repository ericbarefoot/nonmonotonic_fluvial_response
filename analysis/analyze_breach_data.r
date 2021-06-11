# a script to estimate shoreline rugosity measured from topset shorelines.
# Eric Barefoot
# Nov 2020

library(tidyverse)
library(here)
if (interactive()){
  library(colorout)
}
library(R2jags)
library(mcmcplots)
library(tidybayes)


intervalLab <- c(
  `1` = "No Flooding",
  `3` = "High Intensity",
  `1.5` = "Low Intensity"
)

fig_path = here('manuscript', 'figures')

perimData = read_csv(here('data', 'derived_data', 'channel_areas_perimeters.csv')) %>%
mutate(time = floor(time))

breachData10 = read_csv(here('data', 'raw_data', 'breach_counting_10.csv'))
breachData15 = read_csv(here('data', 'raw_data', 'breach_counting_15.csv'))
breachData30 = read_csv(here('data', 'raw_data', 'breach_counting_30.csv'))

breachData = bind_rows(list(breachData10, breachData15, breachData30)) %>%
mutate(t = floor(t))

dat = full_join(perimData, breachData, by = c('time' = 't', 'qv' = 'qv')) %>%
select(-oz)


# dat %>% group_by(qv, area_channel, channel_bank_length, img, time) %>% summarize(n())


each_image = dat %>% mutate(aspect = area_channel / channel_bank_length ) %>%
group_by(aspect, channel_bank_length, qv, img, time) %>%
summarize(n = n()) %>% mutate(nperper = n / aspect, nper = n / channel_bank_length)

write_csv(each_image, here('data', 'derived_data', 'breach_counting_counts.csv'))

avg = each_image %>% group_by(qv) %>%
summarize(
  nper_med = median(nper, na.rm = TRUE),
  nper_q2 = quantile(nper, 0.25, na.rm = TRUE),
  nper_q4 = quantile(nper, 0.75, na.rm = TRUE),
  nperper_med = median(nperper, na.rm = TRUE),
  nperper_q2 = quantile(nperper, 0.25, na.rm = TRUE),
  nperper_q4 = quantile(nperper, 0.75, na.rm = TRUE)
)

# ggplot(data = each_image, aes(x = qv)) +
# geom_boxplot(aes(y = nperper, color = as.factor(qv))) +
# geom_line(aes(x = qv, y = nperper_med), data = avg) +
# # geom_jitter(aes(y = nperper, color = as.factor(qv))) +
# geom_pointrange(aes(x = qv, y = nperper_med, ymax = nperper_q4, ymin = nperper_q2), data = avg) +
# labs(x = 'Flooding Intensity', y = 'Levee breaches per mm of bankline \n (normalized by channel area)', color = 'Flooding Intensity \n State')
#
#
# ggsave(filename = here(fig_path, 'levee_breaches_per_unit_bankline.png'))
#
#
# ggplot(data = each_image, aes(x = qv)) +
# geom_violin(aes(y = nper, color = as.factor(qv))) +
# geom_line(aes(x = qv, y = nper_med), data = avg) +
# geom_jitter(aes(y = nper, color = as.factor(qv))) +
# geom_pointrange(aes(x = qv, y = nper_med, ymax = nper_q4, ymin = nper_q2, color = as.factor(qv)), data = avg)




# kruskal.test(each_image$n, each_image$qv)

#------------------------------------------------------------------------------#
# Now do the rugosity
#------------------------------------------------------------------------------#

breachModel = each_image %>% ungroup %>%
transmute(
  bankLength = channel_bank_length * 0.005,
  count = n,
  qv = factor(qv, levels = c('1', '1.5', '3')),
  # time = factor(time),
  aspect = aspect
) %>% # select(-c(n, nperper, nper, img)) %>%
compose_data()

params_to_save = c(
  'breach.cnt'
  # 'breach.slp.dim'
  # 'mu.breach.cnt',
  # 'breach.int'
)

breachJags = jags(
  model.file = here('analysis', 'levee_breach_abundance_jags.r'),
  data = breachModel, parameters.to.save = params_to_save,
  n.chains = 4, n.iter = 3000, n.burnin = 500, n.thin = 2
)
# breachJags

write_rds(breachJags, here('data', 'derived_data', 'levee_breach_counts_posteriors.rds'))
