# Eric Barefoot
# Feb 2020

library(tidyverse)
library(here)
if (interactive()){
  library(colorout)
}

library(conover.test)

library(R2jags)
library(mcmcplots)
library(tidybayes)


intervalLab <- c(
  `1` = "No Flooding",
  `3` = "High Intensity",
  `1.5` = "Low Intensity"
)

fig_path = here('manuscript', 'figures')

chanData = read_csv(here('analysis', 'chapter_2', 'data', 'channeldimensionsxy.csv'))
slopeData = read_csv(here('analysis', 'chapter_2', 'data', 'ChannelSlopes.csv'))

openA = 35

slopeFilter = slopeData %>% filter(
  slope != 0,
  between(angle, 25, 75)
) %>% mutate(slope = abs(slope))

slopeSumm = slopeFilter %>% group_by(qv) %>%
summarize(
  n = n(),
  mean = mean(slope) * 100,
  std = sd(slope) * 100,
  stdErr = std / sqrt(n),
  med = median(slope) * 100,
  q10 = quantile(slope, 0.10) * 100,
  q25 = quantile(slope, 0.25) * 100,
  q75 = quantile(slope, 0.75) * 100,
  q90 = quantile(slope, 0.90) * 100
) %>% select(-n) %>%
pivot_longer(cols = c(mean), names_to = 'stat', values_to = 'value') %>%
mutate(measurement = 'slope')

# chanFilter = chanData %>%
# filter(width > 0.01, depth > 0, radius == 0.5)

chanFilter = chanData %>% group_by(qv, time) %>%
filter(width > 0.01, depth > quantile(depth, 0.85, na.rm = T))
# %>%
# ggplot() + geom_point(aes(width, depth, color = radius)) + facet_wrap(~qv)

chanDepthSumm = chanFilter %>%
group_by(qv) %>%
summarize(
  n = n(),
  mean = mean(depth) * 1000,
  std = sd(depth) * 1000,
  stdErr = std / sqrt(n),
  med = median(depth) * 1000,
  q10 = quantile(depth, 0.10) * 1000,
  q25 = quantile(depth, 0.25) * 1000,
  q75 = quantile(depth, 0.75) * 1000,
  q90 = quantile(depth, 0.90) * 1000
) %>% select(-n) %>%
pivot_longer(cols = c(mean), names_to = 'stat', values_to = 'value') %>%
mutate(measurement = 'depth')

chanWidthSumm = chanFilter %>%
group_by(qv) %>%
summarize(
  n = n(),
  mean = mean(width) * 100,
  std = sd(width) * 100,
  stdErr = std / sqrt(n),
  med = median(width) * 100,
  q10 = quantile(width, 0.10) * 100,
  q25 = quantile(width, 0.25) * 100,
  q75 = quantile(width, 0.75) * 100,
  q90 = quantile(width, 0.90) * 100
) %>% select(-n) %>%
pivot_longer(cols = c(mean), names_to = 'stat', values_to = 'value') %>%
mutate(measurement = 'width')

dimSumm = bind_rows(list(slopeSumm, chanDepthSumm, chanWidthSumm)) %>%
mutate(
  measurement = factor(
    measurement,
      levels = c('slope', 'depth', 'width'),
      labels = c(
        "Slope(m*m^{-1}%*%10^{-2})",
        "Bankfull~Depth~(mm)",
        "Bankfull~Width~(cm)"
      )
    )
  )



txtsize = 10

meas = c('slope', 'depth', 'width')
ltrs = c('a', 'b', 'c')
panellabels = tibble(
  measurement = factor(c('slope', 'depth', 'width'),
      levels = c('slope', 'depth', 'width'),
      labels = c(
        "Slope(m*m^{-1}%*%10^{-2})",
        "Bankfull~Depth~(mm)",
        "Bankfull~Width~(cm)"
      )
    ),
  panel = ltrs,
  yloc = c(3,16,22)
)

channelGeom_plot = dimSumm %>%
ggplot(aes(x = qv, y = value)) +
geom_line() +
geom_linerange(aes(ymax = value + std, ymin = value - std)) +
# geom_linerange(aes(ymax = q75, ymin = q25), size = 2) +
geom_point(size = 2) +
geom_label(aes(x = 2.5, y = yloc, label = panel), label.size = NA, data = panellabels) +
labs(x = expression(Flood~Variability~(Q[v])), y = NULL) +
facet_wrap(ncol = 3, vars(measurement), scale = 'free_y', labeller = label_parsed, strip.position = 'left') +
theme_minimal() +
theme(
  strip.background = element_blank(),
  legend.position = 'bottom', strip.placement = 'outside',
  plot.title = element_text(hjust = 0.5, size = txtsize),
  strip.text.y = element_text(size = txtsize - 1),
  axis.title.y = element_text(size = txtsize - 1),
  axis.text.y = element_text(size = txtsize - 2),
  legend.title = element_text(size = txtsize - 1),
  legend.text = element_text(size = txtsize - 2),
  axis.title.x = element_text(size = txtsize - 1),
  axis.text.x = element_text(size = txtsize - 2),
)

ggsave('channel_geoms.pdf', channelGeom_plot, width = 6, height = 2, path = fig_path)

#------------------------------------------------------------------------------#
# now to model the results
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Have to do separately, so slope first
#------------------------------------------------------------------------------#

slopeModel = slopeFilter %>%
# select(-c(xrb, yrb, xlb, ylb, tID)) %>%
mutate(
  qv = factor(qv, levels = c('1', '1.5', '3')),
  time = factor(time),
  # radius = factor(radius),
  angle = factor(round(angle)),
) %>% sample_frac(0.2) %>%
compose_data()

params_to_save = c(
  # 'slope.mod'
  'mu.slope', 'sd.slope.resid'
)

chanSlopeJags = jags(
  model.file = here('analysis', 'chapter_2',  'channel_slope_trends_jags.r'),
  data = slopeModel, parameters.to.save = params_to_save,
  n.chains = 4, n.iter = 5000, n.burnin = 1000, n.thin = 2
)

write_rds(chanSlopeJags, here('analysis', 'chapter_2', 'data', 'channel_slopes_posteriors.rds'))


#------------------------------------------------------------------------------#
# Now do the channel dimensions
#------------------------------------------------------------------------------#

chanDimsModelData = chanFilter %>%
select(-c(xrb, yrb, xlb, ylb, tID)) %>%
mutate(
  qv = factor(qv, levels = c('1', '1.5', '3')),
  time = factor(time),
  radius = factor(radius),
  # angle = factor(round(angle)),
) %>% compose_data()

params_to_save = c(
  # 'depth.mod',
  # 'width.mod'
  'mu.depth', 'sd.depth.resid',
  'mu.width', 'sd.width.resid'
)

chanDimJags = jags(
  model.file = here('analysis', 'chapter_2',  'channel_dimensions_trends_jags.r'),
  data = chanDimsModelData, parameters.to.save = params_to_save,
  n.chains = 4, n.iter = 3000, n.burnin = 500, n.thin = 2
)

write_rds(chanDimJags, here('analysis', 'chapter_2', 'data', 'channel_dimensions_posteriors.rds'))


# /home/eric/Dropbox/research/transition_time_experiments/analysis/chapter_2/channel_geometry_trends_jags.r

# mcmcplot(biasJags)

# slopeFilter %>% group_by(qv) %>%  mutate(time = time - min(time)) %>% ggplot() +
# geom_point(aes(angle, slope, color = time)) + facet_wrap(~qv)
#
# chanFilter %>% group_by(qv) %>%  mutate(time = time - min(time)) %>%
# ggplot() + geom_point(aes(width, depth, color = time)) + facet_grid(radius~qv)
