
library(tidyverse)
library(modelr)
library(broom)
library(rsample)
library(here)
if (interactive()){
  library(colorout)
}
library(gridExtra)

fig_path = here('manuscript', 'figures')

ov_dat10 = read_csv(here('analysis', 'chapter_2', 'data', 'qv10overlap.csv'))
ov_dat15 = read_csv(here('analysis', 'chapter_2', 'data', 'qv15overlap.csv'))
ov_dat30 = read_csv(here('analysis', 'chapter_2', 'data', 'qv30overlap.csv'))

ov_dat = bind_rows(list(qv10 = ov_dat10, qv15 = ov_dat15, qv30 = ov_dat30), .id = 'interval')

ov_data = ov_dat %>% mutate(dt = (t1 - t0), interval = as.factor(interval)) %>%
filter(ov > 0)

# give a guess for initial parameters
theta_0 = min(ov_data$ov) * 0.5

# Estimate the rest parameters using a linear model

model_0 = lm(log(ov - theta_0) ~ dt, data = ov_data)
alpha_0 = exp(coef(model_0)[1])
beta_0 = coef(model_0)[2]

start_0 = list(alpha = alpha_0, beta = beta_0, theta = theta_0)

andyModel = function(df, start_0) {
  nls(ov ~ alpha * exp(beta * dt) + theta, data = df, start = start_0)
}

interval_fits = ov_data %>%
group_by(interval) %>% nest() %>%
mutate(
  model = map(data, andyModel, start_0),
  coeffs = map(model, tidy, conf.int = TRUE),
  pred = map(model, augment)
)

predictions = interval_fits %>% unnest(., pred) %>%
select(interval, dt, ov, .fitted) %>% mutate(
  qv = case_when(
    interval == 'qv10' ~ 1,
    interval == 'qv15' ~ 1.5,
    interval == 'qv30' ~ 3
  ),
  dataset = interval
)

intervalLab <- c(
  qv10 = "No Flooding",
  qv30 = "High Intensity",
  qv15 = "Low Intensity"
)

coeffs = interval_fits %>% unnest(., coeffs)

overlap_decay_summary = coeffs %>% select(interval, term, estimate, conf.low, conf.high) %>%
pivot_longer(cols = c(estimate, conf.low, conf.high), names_to = 'confbound', values_to = 'value') %>%
pivot_wider(names_from = term, values_from = value) %>%
mutate(M = -beta, am = alpha + theta, pm = theta) %>% select(interval, confbound, M, am, pm) %>%
pivot_longer(cols = c(M, am, pm), names_to = 'term', values_to = 'estimate') %>%
pivot_wider(names_from = confbound, values_from = estimate) %>%
filter(term == 'M') %>%
mutate(
  qv = case_when(
    interval == 'qv10' ~ 1,
    interval == 'qv15' ~ 1.5,
    interval == 'qv30' ~ 3
  ),
  eFoldingT = 1/estimate,
  eFoldingTlwr = 1/conf.low,
  eFoldingTupr = 1/conf.high,
  dataset = 'summary'
) %>%
select(-c(estimate, conf.low, conf.high, term))

bin_mean = ov_data %>%
mutate(bin = cut_width(dt, width = 2, center = 1)) %>%
group_by(qv, bin) %>%
summarize(dt = mean(dt), mov = mean(ov)) %>%
mutate(dataset = factor(
  case_when(
  qv == 1 ~ 'qv10',
  qv == 1.5 ~ 'qv15',
  qv == 3 ~ 'qv30'
),
levels = c('qv10', 'qv15', 'qv30', 'summary')
))

tits = c('No Flooding', 'Low-Intensity Flooding', 'High-Intensity Flooding')

inter = c('qv10', 'qv15', 'qv30')

ov_binned = ov_data %>% mutate(bins = cut_width(dt, width = 2)) %>%
group_by(qv, bins) %>% summarize(
  dt = mean(dt),
  ov_75 = quantile(ov, 0.75),
  ov_25 = quantile(ov, 0.25),
  ov_10 = quantile(ov, 0.10),
  ov_90 = quantile(ov, 0.90),
  ov_01 = quantile(ov, 0.01),
  ov_99 = quantile(ov, 0.99),
  ov = median(ov),
  dataset = case_when(
      qv == 1 ~ 'qv10',
      qv == 1.5 ~ 'qv15',
      qv == 3 ~ 'qv30',
    )
)

txtsize = 10

ltrs = c('a', 'b', 'c', 'd')

decay_plots = list()
j = 1
for (i in inter) {
  pt_data = predictions %>% filter(dataset == i)
  ovBinData = ov_binned %>% filter(dataset == i)
  mvpt_data = bin_mean %>% filter(dataset == i)
  decay_plots[[j]] = ggplot() +
  # geom_point(aes(dt, ov), alpha = 0.1, data = pt_data) +
  geom_ribbon(aes(dt / 50, ymax = ov_99, ymin = ov_01), alpha = 1, fill = '#cccccc', data = ovBinData) +
  geom_ribbon(aes(dt / 50, ymax = ov_90, ymin = ov_10), alpha = 1, fill = '#969696', data = ovBinData) +
  geom_ribbon(aes(dt / 50, ymax = ov_75, ymin = ov_25), alpha = 1, fill = '#636363', data = ovBinData) +
  # geom_line(aes(dt, ov), data = ovBinData, color = 'black', size = 1) +
  # geom_point(aes(dt, mov), color = 'red3', data = mvpt_data) +
  geom_line(aes(dt / 50, .fitted), col = 'red3', size = 1.5, data = pt_data) +
  labs(x = expression(Time~(T[c])), y = 'Channel Overlap (-)', title = tits[j]) +
  # scale_y_continuous(trans = 'log') +
  # scale_x_continuous(trans = 'log') +
  # annotate('rect', x = 1, y = 1, label = ltrs[j]) +
  # annotate('text', x = 1, y = 1, label = ltrs[j]) +
  geom_label(aes(x = 1, y = 0.9), label.size = NA, label = ltrs[j]) +
  xlim(c(0, 50 / 50)) +
  ylim(c(0, 1)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = txtsize),
    strip.text.y = element_text(size = txtsize - 1),
    axis.title.y = element_text(size = txtsize - 1),
    axis.text.y = element_text(size = txtsize - 2),
    legend.title = element_text(size = txtsize - 1),
    legend.text = element_text(size = txtsize - 2),
    axis.title.x = element_text(size = txtsize - 1),
    axis.text.x = element_text(size = txtsize - 2),
  )
  j = j + 1
}

ptrange_plot =
ggplot(overlap_decay_summary) +
geom_line(aes(x = qv, y = eFoldingT / 50)) +
geom_pointrange(aes(x = qv, y = eFoldingT / 50,
  ymax = eFoldingTupr / 50, ymin = eFoldingTlwr / 50)) +
labs(
  x = expression(Flood~Intensity~Q[v]),
  y = expression(italic(e)-Folding~Time~(T[c])),
  # y = 'Channel Decorrelation \n e-Folding Time (hrs)',
  title = '') +
theme_minimal() +
geom_label(aes(x = 3, y = 0.5), label.size = NA, label = ltrs[4]) +
theme(
  plot.title = element_text(hjust = 0.5, size = txtsize),
  strip.text.y = element_text(size = txtsize - 1),
  axis.title.y = element_text(size = txtsize - 1),
  axis.text.y = element_text(size = txtsize - 2),
  legend.title = element_text(size = txtsize - 1),
  legend.text = element_text(size = txtsize - 2),
  axis.title.x = element_text(size = txtsize - 1),
  axis.text.x = element_text(size = txtsize - 2),
)



combo_plot = grid.arrange(
  decay_plots[[1]],
  decay_plots[[2]],
  decay_plots[[3]],
  ptrange_plot,
  nrow = 4
)

ggsave('channel_overlap_decay.pdf', combo_plot, path = fig_path,
  height = 9, width = 3, units = 'in', dpi = 200
)
