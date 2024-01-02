saveRDS(snakemake, '.figure4.R.RDS')
# snakemake <- readRDS('.figure4.R.RDS')

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(zoo)
library(wesanderson)
library(RColorBrewer)

colors <- wes_palette('Zissou1')

###############
# P1 (a)
###############

event_data <- fread(snakemake@input$event_data)
event_data <- event_data %>%
  .[ , Time := Time - unique(EventStart) ] %>%
  .[ , EventStart := 0 ] %>%
  .[ , Type := 'None' ] %>%
  .[ Time >= -0.4 & Time <= -0.2, Type := 'B' ] %>%
  .[ Time >= -0.1 & Time <= 0, Type := 'P' ] %>%
  .[ Time > 0 & Time <= 0.2, Type := 'O' ] %>%
  .[ Time > 0.2 & Time <= 0.4, Type := 'L' ] %>%
  .[ , Region := gsub(x = Region, pattern = '_', replacement = ' ') ] %>%
  .[ , Region := gsub(x = Region, pattern = '23', replacement = '2/3') ] %>%
  .[ , Region := factor(Region, levels = c('S1 L2/3', 'S1 L5', 'M1 L2/3', 'M1 L5')) ]

average_data <- event_data %>%
  .[ , list(TKEO = mean(TKEO), VM = mean(VM), EMG = mean(EMG)), by = list(Region, Time) ]

p1 <- average_data %>%
  .[ , list(VM = rollmean(VM, 250), Time = rollmean(Time, 250)), by = Region ] %>%
  ggplot() +
  annotate('rect', xmin = -0.4, xmax = -0.2, ymin = -Inf, ymax = Inf, fill = colors[ 1 ], size = 0.5, alpha = 0.5) +
  annotate('rect', xmin = -0.1, xmax = 0, ymin = -Inf, ymax = Inf, fill = colors[ 2 ], size = 0.5, alpha = 0.5) +
  annotate('rect', xmin = 0, xmax = 0.1, ymin = -Inf, ymax = Inf, fill = colors[ 3 ], size = 0.5, alpha = 0.5) +
  annotate('rect', xmin = 0.2, xmax = 0.4, ymin = -Inf, ymax = Inf, fill = colors[ 4 ], size = 0.5, alpha = 0.5) +
  geom_line(aes(x = Time, y = VM, color = Region), linewidth = 1) +
  theme_light(base_size = 8) +
  theme(legend.position = 'top') +
  scale_color_manual(name = '', values = brewer.pal(11, 'RdBu')[ c(1, 2, 10, 11) ]) +
  xlab('Time, t (s)') +
  ylab('Membrane potential, V (mV)') +
  annotate('text', x = -0.3, y = -47, label = 'B', size = 2.5) +
  annotate('text', x = -0.05, y = -47, label = 'P', size = 2.5) +
  annotate('text', x = 0.05, y = -47, label = 'O', size = 2.5) +
  annotate('text', x = 0.3, y = -47, label = 'L', size = 2.5)

###############
# P2 (b)
###############

vm <- fread(snakemake@input$vm)

averages <- vm %>%
  .[ , list(Average = mean(Vm), SD = sd(Vm)), by = list(SID, EventID, Type, Region) ] %>%
  .[ , Region := gsub(x = Region, pattern = '_', replacement = ' ') ] %>%
  .[ , Region := gsub(x = Region, pattern = '23', replacement = '2/3') ] %>%
  .[ , Region := factor(Region, levels = c('S1 L2/3', 'S1 L5', 'M1 L2/3', 'M1 L5')) ] %>%
  .[ , Region := factor(Region, levels = c('S1 L2/3', 'S1 L5', 'M1 L2/3', 'M1 L5')) ] %>%
  .[ Type == 'Movement', Type := 'Move' ] %>%
  .[ Type == 'No movement', Type := 'Rest' ] %>%
  .[ , Type := factor(Type, c('Move', 'Rest')) ]

p2 <- averages %>%
  ggplot(aes(x = Region, y = Average)) +
  facet_wrap(~Type) +
  geom_boxplot(aes(color = Region, fill = Region), alpha = 0.25, outlier.alpha = 0.25, outlier.size = 0.5) +
  scale_color_manual(name = '', values = brewer.pal(11, 'RdBu')[ c(1, 2, 10, 11) ]) +
  scale_fill_manual(name = '', values = brewer.pal(11, 'RdBu')[ c(1, 2, 10, 11) ]) +
  theme(legend.position = 'none') +
  stat_compare_means(comparisons = list(
    c('S1 L2/3', 'S1 L5'),
    c('S1 L5', 'M1 L2/3'),
    c('S1 L2/3', 'M1 L2/3'),
    c('S1 L5', 'M1 L5'),
    c('S1 L2/3', 'M1 L5')
  ), method = 't.test', size = 2.5) +
  stat_compare_means(comparisons = list(
    c('M1 L2/3', 'M1 L5')
  ), method = 't.test', size = 2.5) +
  theme_light(base_size = 8) +
  xlab('') +
  theme(legend.position = 'none') +
  ylab('Average membrane potential, V (mV)') +
  theme(strip.background = element_rect(fill = 'white')) +
  theme(strip.text = element_text(colour = 'black'))

###############
# P3 (c)
###############

rest_mean <- averages[ Type == 'Rest', list(RestMean = mean(Average)), by = Region ]

type_average <- event_data %>%
  .[ , list(VM = mean(VM)), by = list(SID, Type, EventId, Region) ] %>%
  .[ Type != 'None' ] %>%
  .[ Type == 'B', Type := 'Baseline' ] %>%
  .[ Type == 'P', Type := 'Pre-movement' ] %>%
  .[ Type == 'O', Type := 'Movement onset' ] %>%
  .[ Type == 'L', Type := 'Late movement' ] %>%
  .[ , Type := factor(Type, levels = c('Baseline', 'Pre-movement', 'Movement onset', 'Late movement')) ] %>%
  .[ , Region := gsub(x = Region, pattern = '_', replacement = ' ') ] %>%
  .[ , Region := gsub(x = Region, pattern = '23', replacement = '2/3') ] %>%
  .[ , Region := factor(Region, levels = c('S1 L2/3', 'S1 L5', 'M1 L2/3', 'M1 L5')) ] %>%
  .[ Region == 'S1 L2/3', VM := VM - rest_mean[ Region == 'S1 L2/3', RestMean ] ] %>%
  .[ Region == 'S1 L5', VM := VM - rest_mean[ Region == 'S1 L5', RestMean ] ] %>%
  .[ Region == 'M1 L2/3', VM := VM - rest_mean[ Region == 'M1 L2/3', RestMean ] ] %>%
  .[ Region == 'M1 L5', VM := VM - rest_mean[ Region == 'M1 L5', RestMean ] ]

p3 <- type_average %>%
  ggplot(aes(x = Region, y = VM, color = Region, fill = Region)) +
  facet_wrap(~Type, ncol = 4) +
  geom_boxplot(alpha = 0.25, outlier.alpha = 0.25, outlier.size = 0.5) +
  stat_compare_means(
    hide.ns = TRUE,
    comparisons = list(
      c('S1 L2/3', 'S1 L5'),
      c('S1 L5', 'M1 L2/3'),
      c('S1 L2/3', 'M1 L2/3'),
      c('S1 L5', 'M1 L5'),
      c('S1 L2/3', 'M1 L5')
    ),
    method = 't.test',
    size = 2.5
  ) +
  stat_compare_means(
    hide.ns = TRUE,
    comparisons = list(c('M1 L2/3', 'M1 L5')),
    method = 't.test',
    size = 2.5
  ) +
  theme_light(base_size = 8) +
  scale_color_manual(name = '', values = brewer.pal(11, 'RdBu')[ c(1, 2, 10, 11) ]) +
  scale_fill_manual(name = '', values = brewer.pal(11, 'RdBu')[ c(1, 2, 10, 11) ]) +
  xlab('') +
  theme(legend.position = 'none') +
  theme(strip.background = element_rect(fill = 'white')) +
  theme(strip.text = element_text(colour = 'black')) +
  ylab('Membrane potential change, V (mV)')

final <- ggarrange(
  ggarrange(p1, p2, labels = letters[ 1:2 ], font.label = list(size = 9)),
  ggarrange(p3, labels = letters[ 3 ], font.label = list(size = 9)),
  nrow = 2
)

ggsave(final, filename = snakemake@output$png, height = 6, width = 8, bg = 'white', dpi = 600)
