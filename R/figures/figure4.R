saveRDS(snakemake, '.figure4.R.RDS')
# snakemake <- readRDS('.figure4.R.RDS')

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(zoo)
library(wesanderson)
library(RColorBrewer)
library(foreach)

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
  ggplot() +
  annotate('rect', xmin = -0.4, xmax = -0.2, ymin = -Inf, ymax = Inf, fill = colors[ 1 ], size = 0.5, alpha = 0.5) +
  annotate('rect', xmin = -0.1, xmax = 0, ymin = -Inf, ymax = Inf, fill = colors[ 2 ], size = 0.5, alpha = 0.5) +
  annotate('rect', xmin = 0, xmax = 0.1, ymin = -Inf, ymax = Inf, fill = colors[ 3 ], size = 0.5, alpha = 0.5) +
  annotate('rect', xmin = 0.2, xmax = 0.4, ymin = -Inf, ymax = Inf, fill = colors[ 4 ], size = 0.5, alpha = 0.5) +
  geom_line(aes(x = Time, y = VM, color = Region), linewidth = 0.5) +
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

type_average <- event_data %>%
  .[ , list(VM = mean(VM)), by = list(SID, Type, Region) ] %>%
  .[ Type != 'None' ] %>%
  .[ Type == 'B', Type := 'Baseline' ] %>%
  .[ Type == 'P', Type := 'Pre-movement' ] %>%
  .[ Type == 'O', Type := 'Movement onset' ] %>%
  .[ Type == 'L', Type := 'Late movement' ] %>%
  .[ , Type := factor(Type, levels = c('Baseline', 'Pre-movement', 'Movement onset', 'Late movement')) ] %>%
  .[ , Region := gsub(x = Region, pattern = '_', replacement = ' ') ] %>%
  .[ , Region := gsub(x = Region, pattern = '23', replacement = '2/3') ] %>%
  .[ , Region := factor(Region, levels = c('S1 L2/3', 'S1 L5', 'M1 L2/3', 'M1 L5')) ] %>%
  .[ , Animal := gsub(x = SID, pattern = '(W[0-9]).*', replacement = '\\1') ] %>%
  .[ , Animal := factor(Animal, levels = c('W1', 'W2', 'W3', 'W4')) ]


change <- foreach(sid = type_average[ , unique(SID) ], .combine = rbind) %do% {
  baseline <- type_average[ SID == sid & Type == 'Baseline', VM ]
  premovement <- type_average[ SID == sid & Type == 'Pre-movement', VM ]
  onset <- type_average[ SID == sid & Type == 'Movement onset', VM ]
  late <- type_average[ SID == sid & Type == 'Late movement', VM ]

  data.table(
    SID = sid,
    Type = c('Pre-movement', 'Movement onset', 'Late movement'),
    Change = c(premovement - baseline, onset - baseline, late - baseline),
    Region = type_average[ SID == sid, unique(Region) ]
  )
} %>%
  .[ , Type := factor(Type, levels = c('Baseline', 'Pre-movement', 'Movement onset', 'Late movement')) ]

p2 <- change %>%
  ggplot(aes(x = Region, y = Change, color = Region)) +
  facet_wrap(~Type, ncol = 4) +
  geom_boxplot(alpha = 0.25, outlier.alpha = 0) +
  geom_jitter(height = 0, width = 0.2, alpha = 0.5) +
  stat_compare_means(
    hide.ns = TRUE,
    comparisons = list(
      c('S1 L2/3', 'S1 L5'),
      c('S1 L5', 'M1 L2/3'),
      c('S1 L2/3', 'M1 L2/3'),
      c('S1 L5', 'M1 L5'),
      c('S1 L2/3', 'M1 L5')
    ),
    method = 'wilcox',
    size = 2.5
  ) +
  stat_compare_means(
    hide.ns = TRUE,
    comparisons = list(c('M1 L2/3', 'M1 L5')),
    method = 'wilcox',
    size = 2.5
  ) +
  theme_light(base_size = 8) +
  scale_color_manual(name = '', values = brewer.pal(11, 'RdBu')[ c(1, 2, 10, 11) ]) +
  scale_fill_manual(name = '', values = brewer.pal(11, 'RdBu')[ c(1, 2, 10, 11) ]) +
  xlab('') +
  theme(legend.position = 'none') +
  theme(strip.background = element_rect(fill = 'white')) +
  theme(strip.text = element_text(colour = 'black')) +
  ylab('Membrane potential change, V (mV)') +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

###############
# P3 (c)
###############

p3 <- type_average %>%
  ggplot(aes(x = Region, y = VM, color = Region)) +
  facet_wrap(~Type, ncol = 4) +
  geom_boxplot(alpha = 0.25, outlier.alpha = 0) +
  geom_jitter(height = 0, width = 0.2, alpha = 0.5) +
  stat_compare_means(
    hide.ns = TRUE,
    comparisons = list(
      c('S1 L2/3', 'S1 L5'),
      c('S1 L5', 'M1 L2/3'),
      c('S1 L2/3', 'M1 L2/3'),
      c('S1 L5', 'M1 L5'),
      c('S1 L2/3', 'M1 L5')
    ),
    method = 'wilcox',
    size = 2.5
  ) +
  stat_compare_means(
    hide.ns = TRUE,
    comparisons = list(c('M1 L2/3', 'M1 L5')),
    method = 'wilcox',
    size = 2.5
  ) +
  theme_light(base_size = 8) +
  scale_color_manual(name = '', values = brewer.pal(11, 'RdBu')[ c(1, 2, 10, 11) ]) +
  scale_fill_manual(name = '', values = brewer.pal(11, 'RdBu')[ c(1, 2, 10, 11) ]) +
  xlab('') +
  theme(legend.position = 'none') +
  theme(strip.background = element_rect(fill = 'white')) +
  theme(strip.text = element_text(colour = 'black')) +
  ylab('Membrane potential, V (mV)') +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

final <- ggarrange(
  ggarrange(p1, p2, labels = letters[ 1:2 ], font.label = list(size = 9), common.legend = TRUE, widths = c(0.4, 0.6)),
  ggarrange(p3, labels = letters[ 3 ], font.label = list(size = 9)),
  nrow = 2
)

ggsave(final, filename = snakemake@output$png, height = 6, width = 8, bg = 'white', dpi = 600)
