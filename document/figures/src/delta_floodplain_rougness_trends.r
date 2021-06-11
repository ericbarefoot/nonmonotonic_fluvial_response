# Eric Barefoot
# Feb 2020

library(tidyverse)
library(here)
if (interactive()) {
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


Rugdata = read_csv(here('analysis', 'chapter_2', 'data', 'shorelineRugosity.csv'))
roughness = read_csv(here('analysis', 'chapter_2', 'data', 'FPrelief.csv'))
FPslp = read_csv(here('analysis', 'chapter_2', 'data', 'FPSlopes.csv'))

# Rugdata = read_csv(here('data/raw_data/shorelineRugosity.csv'))
# roughness = read_csv(here('data/raw_data/FPrelief.csv'))
# FPslp = read_csv(here('data/raw_data/FPSlopes.csv'))

rugSumm = Rugdata %>% pivot_longer(cols = -c(qv, time), names_to = 'type', values_to = 'rug') %>% group_by(type, qv) %>% summarize(
  n = n(),
  mean = mean(rug),
  std = sd(rug),
  stdErr = std / sqrt(n)
) %>%
mutate(type = factor(type,
  levels = c('rugosity'),
  labels = c(
    "Shoreline~Rugosity"
    )
  )
)

ruffSumm = roughness %>% pivot_longer(cols = -c(qv, time), names_to = 'type', values_to = 'relief') %>%
group_by(type, qv) %>% filter(type == 'short_wave_relief') %>%
summarize(
  n = n(),
  mean = mean(relief) * 1000,
  std = sd(relief) * 1000,
  stdErr = std / sqrt(n)
) %>%
mutate(type = factor(type,
  levels = c('short_wave_relief', 'long_wave_relief'),
  labels = c(
    "Floodplain~Relief~(mm)",
    'Long-Wavelength~Relief~(mm)'
    )
  )
)

# slpfilter = FPslp %>% filter(
#   slope != 0,
#   between(angle, 0, 90)
# )

slpSumm = FPslp %>% mutate(slope = abs(slope)) %>%
pivot_longer(cols = -c(qv, time, angle), names_to = 'type', values_to = 'relief') %>% group_by(type, qv) %>% summarize(
  n = n(),
  mean = mean(relief) * 100,
  std = sd(relief) * 100,
  stdErr = std / sqrt(n)
) %>%
mutate(type = factor(type,
  levels = c('slope'),
  labels = c(
    "Slope~(m*m^{-1}%*%10^{-2})"
    )
  )
)

allSumm = bind_rows(list(rugSumm, ruffSumm, slpSumm)) %>% mutate(
  type = ordered(type,
    # levels = c('slope', 'short_wave_relief', 'rugosity'),
    levels = c(
      "Slope~(m*m^{-1}%*%10^{-2})",
      "Floodplain~Relief~(mm)",
      "Shoreline~Rugosity"
    )
  )
)

txtsize = 10

ltrs = c('a', 'b', 'c')
panellabels = tibble(
  type = ordered(c('slope', 'short_wave_relief', 'rugosity'),
      levels = c('slope', 'short_wave_relief', 'rugosity'),
      labels = c(
        "Slope~(m*m^{-1}%*%10^{-2})",
        "Floodplain~Relief~(mm)",
        "Shoreline~Rugosity"
      )
    ),
  panel = ltrs,
  yloc = c(6.5,8.4,0.38)
)

channelGeom_plot = allSumm %>%
ggplot(aes(x = qv, y = mean)) +
geom_line() +
geom_linerange(
  aes(
    ymax = (mean + std),
    ymin = (mean - std)
  )
) +
geom_point(size = 2) +
labs(x = expression(Flood~Variability~(Q[v])), y = NULL) +
facet_wrap(ncol = 3, vars(type), scale = 'free_y', labeller = label_parsed, strip.position = 'left') +
theme_minimal() +
geom_label(aes(x = 2.8, y = yloc, label = panel), label.size = NA, data = panellabels) +
theme(
  strip.background = element_blank(), legend.position = 'bottom', strip.placement = 'outside',
  plot.title = element_text(hjust = 0.5, size = txtsize),
  strip.text.y = element_text(size = txtsize - 1),
  axis.title.y = element_text(size = txtsize - 1),
  axis.text.y = element_text(size = txtsize - 2),
  legend.title = element_text(size = txtsize - 1),
  legend.text = element_text(size = txtsize - 2),
  axis.title.x = element_text(size = txtsize - 1),
  axis.text.x = element_text(size = txtsize - 2),
)

# channelGeom_plot

ggsave('floodplain_roughness.pdf', channelGeom_plot, width = 6, height = 2, path = fig_path)


#------------------------------------------------------------------------------#
# now to model the results
#------------------------------------------------------------------------------#

# Rugdata
# roughness
# FPslp

#------------------------------------------------------------------------------#
# Have to do separately, so slope first
#------------------------------------------------------------------------------#

slopeModel = FPslp %>% mutate(slope = abs(slope)) %>%
mutate(
  qv = factor(qv, levels = c('1', '1.5', '3')),
  time = factor(time),
  # radius = factor(radius),
  angle = factor(round(angle)),
) %>% sample_frac(0.1) %>%
compose_data()

params_to_save = c(
  # 'slope.mod'
  'mu.slope', 'sd.slope.resid'
)

FPSlopeJags = jags(
  model.file = here('analysis', 'chapter_2', 'floodplain_slope_trends_jags.r'),
  data = slopeModel, parameters.to.save = params_to_save,
  n.chains = 4, n.iter = 5000, n.burnin = 1000, n.thin = 2
)

write_rds(FPSlopeJags, here('analysis', 'chapter_2', 'data', 'floodplain_slopes_posteriors.rds'))

# #------------------------------------------------------------------------------#
# # Now do the rugosity
# #------------------------------------------------------------------------------#

rugosityModel = Rugdata %>%
mutate(
  qv = factor(qv, levels = c('1', '1.5', '3')),
  time = factor(time)
) %>% compose_data()

params_to_save = c(
  'rugosity.mod', 'sd.rugosity.resid'
  # 'mu.rugosity'
)

rugosityJags = jags(
  model.file = here('analysis', 'chapter_2', 'floodplain_rugosity_trends_jags.r'),
  data = rugosityModel, parameters.to.save = params_to_save,
  n.chains = 4, n.iter = 3000, n.burnin = 500, n.thin = 2
)

write_rds(rugosityJags, here('analysis', 'chapter_2', 'data', 'shoreline_rugosity_posteriors.rds'))


# #------------------------------------------------------------------------------#
# # Now do the relief
# #------------------------------------------------------------------------------#

roughnessModel = roughness %>%
mutate(
  qv = factor(qv, levels = c('1', '1.5', '3')),
  time = factor(time)
) %>% compose_data()

params_to_save = c(
  'roughness.mod', 'sd.roughness.resid'
  # 'mu.rougness'
)

roughnessJags = jags(
  model.file = here('analysis', 'chapter_2', 'floodplain_roughness_trends_jags.r'),
  data = roughnessModel, parameters.to.save = params_to_save,
  n.chains = 4, n.iter = 3000, n.burnin = 500, n.thin = 2
)

write_rds(roughnessJags, here('analysis', 'chapter_2', 'data', 'floodplain_roughness_posteriors.rds'))
