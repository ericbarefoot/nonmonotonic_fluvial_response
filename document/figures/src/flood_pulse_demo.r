library(tidyverse)
library(here)

dd = function(x, a) {
  1/(abs(a) * sqrt(pi)) * exp(-(x/a)^2)
}

Tcycle = 30*60
tres = 1/60
Qcon =  2.721024974 * 0.06309017# gal/min #0.1716 # l/s #
Qbase = Qcon * 0.825 # 82.5% of constant flow equivalent
PulseFactor = 3
Qpeak = PulseFactor * Qbase  # considered generally high for floods. flexible

tt = seq(-Tcycle/2, Tcycle/2, length = Tcycle * tres)
mult = (Qpeak - Qbase) / (Tcycle * (Qcon - Qbase))
a = 1/(mult * sqrt(pi))
curv = dd(tt,a) * (Qcon - Qbase) * Tcycle + Qbase

qv = c(1, 1.5, 3)

Q = c()

for (i in qv) {
    if (i == 1) {
        Q = append(Q, rep(Qcon, length(tt)))
    } else {
        Qpeak = i * Qbase
        mult = (Qpeak - Qbase) / (Tcycle * (Qcon - Qbase))
        a = 1/(mult * sqrt(pi))
        Q = append(Q,  dd(tt, a) * (Qcon - Qbase) * Tcycle + Qbase)
    }
}

flooddata = tibble(
  time = rep(tt, 3),
  # qv = factor(rep(qv, each = length(tt)), levels = c('1', '1.5', '3'), labels = c(
  #     expression(Q[v]==1),
  #     expression(Q[v]==1.5),
  #     expression(Q[v]==3)
  # )),
  qv = factor(rep(qv, each = length(tt)), levels = c('1', '1.5', '3'), labels = c(
      'No Flooding',
      'Low-Intensity Flooding',
      'High-Intensity Flooding'
  )),
  qCon = rep(c(NA, Qcon, Qcon), each = length(tt)),
  q = Q
) %>% pivot_longer(cols = c(qCon, q), names_to = 'curve', values_to = 'Q')

txtsize = 10

floods_plot =
ggplot(flooddata) +
geom_step(aes(x = time * tres + 15, y = Q, group = curve, linetype = curve)) +
labs(x = 'Time (min)', y = expression(Discharge~(L~s^{-1}))) +
theme_minimal() +
theme(
  plot.title = element_text(hjust = 0.5, size = txtsize),
  strip.text.y = element_text(size = txtsize - 1),
  strip.text.x = element_text(size = txtsize - 1),
  axis.title.y = element_text(size = txtsize - 1),
  axis.text.y = element_text(size = txtsize - 2),
  legend.title = element_text(size = txtsize - 1),
  legend.text = element_text(size = txtsize - 2),
  axis.title.x = element_text(size = txtsize - 1),
  axis.text.x = element_text(size = txtsize - 2),
  strip.background = element_blank(),
  legend.position = 'right',
  strip.placement = 'outside'
) +
scale_linetype_manual(
  breaks = c('q', 'qCon'),
  values = c('solid', 'dashed'),
  labels = c("Discharge", "Equivalent\nConstant\nDischarge"),
  guide_legend(title = '')
) +
facet_wrap(ncol = 3, vars(qv))
# facet_wrap(ncol = 3, vars(qv), labeller = label_parsed)

ggsave(here('analysis', 'chapter_2', 'figures', 'output', 'flood_pulse_demo.pdf'),
  width = 6, height = 2
)
#
# ggsave(here('analysis', 'chapter_2', 'figures', 'output', 'flood_pulse_demo.png'),
#   width = 10, height = 4
# )
