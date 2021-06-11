
library(tidyverse)
library(here)
if (interactive()){
  library(colorout)
}
library(gridExtra)
# suppressPackageStartupMessages(
library(R2jags)
library(mcmcplots)
# library(tidybayes)

# plot of completeness estimate

qvLab <- c(
  `1` = "No Flooding",
  `3` = "High Intensity",
  `1.5` = "Low Intensity"
)

fig_path = here('manuscript', 'figures')

# multilevel = read_rds(here('nonmonotonic_chapter', 'data', 'derived_data', 'qvAll_midR_estimates.rds'))
multilevel = read_rds(here('data', 'derived_data', 'sadler_exponent_posteriors.rds'))
# modelData = read_rds(here('nonmonotonic_chapter', 'data', 'derived_data', 'qvAll_midR_modelData.rds'))

completeness = read_csv(
  here('analysis', 'chapter_2', 'data', 'strat_completeness.csv'),
  col_types = 'nnnnnn'
) %>% filter(dt < 100) %>% group_by(qv) %>% sample_frac(0.5)
# ) %>% filter(dt < 100, rLabel == 3) %>% group_by(qv) %>% sample_frac(0.5)

ss = sample((multilevel$BUGSoutput$n.keep * multilevel$BUGSoutput$n.chains), 200)
qvs = c(1, 1.5, 3)
dt = seq(1, 150, 0.5)

sample_data = expand_grid(dt = dt, ss = ss, qv = qvs) %>%
mutate(
  beta0 = case_when(
    qv == 1 ~ multilevel$BUGSoutput$sims.list$mu.beta0[ss, 1],
    qv == 1.5 ~ multilevel$BUGSoutput$sims.list$mu.beta0[ss, 2],
    qv == 3 ~ multilevel$BUGSoutput$sims.list$mu.beta0[ss, 3]
  ),
  beta1 = case_when(
    qv == 1 ~ multilevel$BUGSoutput$sims.list$mu.beta1[ss, 1],
    qv == 1.5 ~ multilevel$BUGSoutput$sims.list$mu.beta1[ss, 2],
    qv == 3 ~ multilevel$BUGSoutput$sims.list$mu.beta1[ss, 3]
  ),
  # intercept = case_when(
  #   qv == 1 ~ multilevel$BUGSoutput$sims.list$interOne[ss, 1],
  #   qv == 1.5 ~ multilevel$BUGSoutput$sims.list$interOne[ss, 2],
  #   qv == 3 ~ multilevel$BUGSoutput$sims.list$interOne[ss, 3]
  # ),
  f.pred = beta0 * dt ^ beta1
)

binned_predictions = sample_data %>% mutate(bins = cut_width(dt, width = 2)) %>%
group_by(qv, bins) %>% summarize(
  dt = mean(dt),
  f_75 = quantile(f.pred, 0.75),
  f_25 = quantile(f.pred, 0.25),
  f_10 = quantile(f.pred, 0.10),
  f_90 = quantile(f.pred, 0.90),
  f.pred = median(f.pred)
)

compl_binned = completeness %>% mutate(bins = cut_width(dt, width = 2)) %>%
group_by(qv, bins) %>% summarize(
  dt = mean(dt),
  f_75 = quantile(f, 0.75),
  f_25 = quantile(f, 0.25),
  f_10 = quantile(f, 0.10),
  f_90 = quantile(f, 0.90),
  f_01 = quantile(f, 0.01),
  f_99 = quantile(f, 0.99),
  f = median(f)
)

# intercepts = sample_data %>% group_by(qv) %>%
# summarize(
#   fVal = 1,
#   i_75 = quantile(intercept, 0.75),
#   i_25 = quantile(intercept, 0.25),
#   i_10 = quantile(intercept, 0.10),
#   i_90 = quantile(intercept, 0.90),
#   i_01 = quantile(intercept, 0.01),
#   i_99 = quantile(intercept, 0.99),
#   intercept = median(intercept),
# )

betaOne = sample_data %>% group_by(qv) %>%
summarize(
  fVal = 1,
  b_75 = quantile(beta1, 0.75),
  b_25 = quantile(beta1, 0.25),
  b_10 = quantile(beta1, 0.10),
  b_90 = quantile(beta1, 0.90),
  b_01 = quantile(beta1, 0.01),
  b_99 = quantile(beta1, 0.99),
  beta1 = median(beta1),
)

# modelData_binned = modelData %>% mutate(bins = cut_width(dt, width = 2)) %>%
# group_by(qv, bins) %>% summarize(
#   dt = mean(dt),
#   f_75 = quantile(f, 0.75),
#   f_25 = quantile(f, 0.25),
#   f_10 = quantile(f, 0.10),
#   f_90 = quantile(f, 0.90),
#   f = median(f)
# )

txtsize = 10

ltrs = c('a', 'b', 'c', 'd')

completeness_decay_plots = list()

qv_toPlot = c(1,1.5,3)

tits = c('No Flooding', 'Low-Intensity Flooding', 'High-Intensity Flooding')

for (i in 1:3) {
  ddd = compl_binned %>% mutate(dt = dt / 50) %>% filter(qv == qv_toPlot[i])
  ppp = binned_predictions %>% mutate(dt = dt / 50) %>% filter(qv == qv_toPlot[i])
  completeness_decay_plots[[i]] = ggplot() +
  geom_ribbon(aes(dt, ymax = f_99, ymin = f_01), alpha = 1, fill = '#cccccc', data = ddd) +
  geom_ribbon(aes(dt, ymax = f_90, ymin = f_10), alpha = 1, fill = '#969696', data = ddd) +
  geom_ribbon(aes(dt, ymax = f_75, ymin = f_25), alpha = 1, fill = '#636363', data = ddd) +
  geom_line(aes(dt, f), data = ddd, color = 'black', size = 1) +
  geom_ribbon(aes(dt, ymax = f_75, ymin = f_25), alpha = 0.5, fill = 'red3', data = ppp) +
  geom_line(aes(dt, f.pred), color = 'red3', size = 1, data = ppp) +
  scale_x_log10() +
  scale_y_log10(limits = c(0.01, 1.1)) +
  # geom_point(aes(x = intercept, y = fVal), color = 'red3', data = intercepts, size = 5, shape = 1) +
  # geom_segment(aes(x = intercept, xend = intercept, y = 0.01, yend = fVal), color = 'red3', linetype = 'dotted', data = intercepts) +
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
  ) +
  geom_label(aes(x = 2.5, y = 0.9), label.size = NA, label = ltrs[i]) +
  # ylim(c(0.01, 1.1)) +
  # facet_grid(cols = vars(qv), labeller = labeller(qv = qvLab))
  labs(x = expression(Time~discretization~(T[c])), y = 'Fraction Complete', title = tits[i])
}

beta1_estimate_plot = ggplot() +
geom_line(aes(x = qv, y = beta1), data = betaOne, size = 1) +
geom_linerange(aes(x = qv, ymax = b_90, ymin = b_10, group = qv), data = betaOne, size = 1) +
geom_linerange(aes(x = qv, ymax = b_75, ymin = b_25, group = qv), data = betaOne, size = 2) +
geom_point(aes(x = qv, y = beta1, group = qv), data = betaOne, size = 3) +
labs(x = expression(Flood~Variability~(Q[v])), y = expression(Sadler~Exponent~(beta)), title = '') +
theme_minimal() +
geom_label(aes(x = 2.9, y = 0.54), label.size = NA, label = ltrs[4]) +
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

# intercept_estimate_plot = ggplot() +
# geom_line(aes(x = qv, y = intercept), data = intercepts, size = 1) +
# geom_linerange(aes(x = qv, ymax = i_90, ymin = i_10, group = qv), data = intercepts, size = 1) +
# geom_linerange(aes(x = qv, ymax = i_75, ymin = i_25, group = qv), data = intercepts, size = 2) +
# geom_point(aes(x = qv, y = intercept, group = qv), data = intercepts, size = 3) +
# labs(x = expression(Flood~Variability~(Q[v])), y = 'Compensation Timescale (hrs)', title = '') +
# theme_minimal()

# completeness_decay_plot
# intercept_estimate_plot

combo_plot = grid.arrange(
  completeness_decay_plots[[1]],
  completeness_decay_plots[[2]],
  completeness_decay_plots[[3]],
  beta1_estimate_plot,
  nrow = 4
)

# combo_plot = grid.arrange(
#   completeness_decay_plots,
#   beta1_estimate_plot,
#   widths = c(3, 1) ,
#   ncol = 2
# )


ggsave('strat_completeness_decay.pdf', combo_plot, path = fig_path,
  height = 9, width = 2.5, units = 'in', dpi = 200
)
