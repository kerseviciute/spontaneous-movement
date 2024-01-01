library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(wesanderson)
library(cowplot)

colors <- wes_palette('Zissou1')

example_data <- fread('output/spontaneous-movement/figures/event_data.csv')

example_data <- example_data %>%
  .[ , Time := Time - unique(EventStart) ] %>%
  .[ , EventStart := 0 ] %>%
  .[ , Type := 'None' ] %>%
  .[ Time >= -0.4 & Time <= -0.2, Type := 'Baseline' ] %>%
  .[ Time >= -0.1 & Time <= 0, Type := 'Pre-movement' ] %>%
  .[ Time > 0 & Time <= 0.2, Type := 'Movement onset' ] %>%
  .[ Time > 0.2 & Time <= 0.4, Type := 'Late movement' ]

sample_data <- example_data[ , list(SID, Region) ] %>%
  unique()

average_data <- example_data %>%
  .[ , list(TKEO = mean(TKEO), VM = mean(VM), EMG = mean(EMG)), by = list(Region, Time) ]

plot <- list()

for (region in c('S1_L23', 'S1_L5', 'M1_L23', 'M1_L5')) {
  region_average <- average_data[ Region == region ]

  region_label <- gsub(x = region, pattern = '_', replacement = ' ')
  region_label <- gsub(x = region_label, pattern = '23', replacement = '2/3')

  max_vm <- region_average[ , max(VM) ] + 0.25

  p1 <- region_average %>%
    ggplot() +
    annotate('rect', xmin = -0.4, xmax = -0.2, ymin = -Inf, ymax = Inf, fill = colors[ 1 ], size = 0.5, alpha = 0.5) +
    annotate('rect', xmin = -0.1, xmax = 0, ymin = -Inf, ymax = Inf, fill = colors[ 2 ], size = 0.5, alpha = 0.5) +
    annotate('rect', xmin = 0, xmax = 0.1, ymin = -Inf, ymax = Inf, fill = colors[ 3 ], size = 0.5, alpha = 0.5) +
    annotate('rect', xmin = 0.2, xmax = 0.4, ymin = -Inf, ymax = Inf, fill = colors[ 4 ], size = 0.5, alpha = 0.5) +
    geom_line(aes(x = Time, y = VM), color = colors[ 5 ], linewidth = 0.25) +
    theme_light(base_size = 8) +
    xlab('Time, t (s)') +
    ylab('Membrane\npotential, V (mV)') +
    xlim(-0.4, 0.4) +
    ggtitle(region_label) +
    theme(plot.title = element_text(hjust = 0.5, size = 9, face = 'bold')) +
    scale_x_continuous(breaks = c(-0.4, -0.2, 0, 0.2, 0.4)) +
    annotate('text', x = -0.3, y = max_vm, label = 'B', size = 2.5) +
    annotate('text', x = -0.05, y = max_vm, label = 'P', size = 2.5) +
    annotate('text', x = 0.05, y = max_vm, label = 'O', size = 2.5) +
    annotate('text', x = 0.3, y = max_vm, label = 'L', size = 2.5)

  type_average <- example_data[ Region == region ] %>%
    .[ , list(VM = mean(VM)), by = list(SID, Type, EventId) ] %>%
    .[ Type != 'None' ] %>%
    .[ Type == 'Baseline', Type := 'B' ] %>%
    .[ Type == 'Pre-movement', Type := 'P' ] %>%
    .[ Type == 'Movement onset', Type := 'O' ] %>%
    .[ Type == 'Late movement', Type := 'L' ] %>%
    .[ , Type := factor(Type, levels = c('B', 'P', 'O', 'L')) ]

  p2 <- type_average %>%
    ggplot(aes(x = Type, y = VM, color = Type)) +
    geom_boxplot(outlier.alpha = 0) +
    geom_jitter(height = 0, width = 0.2, alpha = 0.2) +
    theme_light(base_size = 8) +
    xlab('') +
    scale_colour_manual(values = c('B' = colors[ 1 ], 'P' = colors[ 2 ], 'O' = colors[ 3 ], 'L' = colors[ 4 ])) +
    stat_compare_means(
      comparisons = list(c('B', 'P'), c('P', 'O')),
      method = 't.test', size = 2.5) +
        stat_compare_means(
      comparisons = list(c('O', 'L')),
      method = 't.test', size = 2.5) +
    theme(legend.position = 'none') +
    ylab('Average membrane potential, V (mV)')

  ap_time <- fread('figures/all_ap_time.csv') %>%
    .[ , SID := gsub(x = Trial, pattern = '.*(W[0-9]+_C[0-9]+)', replacement = '\\1') ] %>%
    merge(sample_data) %>%
    .[ , Time := Time - 0.4 ] %>%
    .[ Region == region ] %>%
    .[ , Trial := as.factor(Trial) %>% as.numeric() ]

  p3 <- ap_time %>%
    ggplot() +
    annotate('rect', xmin = -0.4, xmax = -0.2, ymin = -Inf, ymax = Inf, fill = colors[ 1 ], size = 0.5, alpha = 0.5) +
    annotate('rect', xmin = -0.1, xmax = 0, ymin = -Inf, ymax = Inf, fill = colors[ 2 ], size = 0.5, alpha = 0.5) +
    annotate('rect', xmin = 0, xmax = 0.1, ymin = -Inf, ymax = Inf, fill = colors[ 3 ], size = 0.5, alpha = 0.5) +
    annotate('rect', xmin = 0.2, xmax = 0.4, ymin = -Inf, ymax = Inf, fill = colors[ 4 ], size = 0.5, alpha = 0.5) +
    geom_histogram(aes(x = Time), fill = 'white', color = 'black', binwidth = 0.025) +
    theme_light(base_size = 8) +
    ylab('AP count') +
    xlab('Time, t (s)') +
    theme(strip.background = element_rect(fill = 'white')) +
    theme(strip.text = element_text(colour = 'black')) +
    scale_x_continuous(breaks = c(-0.4, -0.2, 0, 0.2, 0.4)) +
    expand_limits(y = c(0, 70)) +
    annotate('text', x = -0.3, y = 69, label = 'B', size = 2.5) +
    annotate('text', x = -0.05, y = 69, label = 'P', size = 2.5) +
    annotate('text', x = 0.05, y = 69, label = 'O', size = 2.5) +
    annotate('text', x = 0.3, y = 69, label = 'L', size = 2.5)

  max_trial <- ap_time[ , unique(Trial) ] %>% length() + 1

  p4 <- ap_time %>%
    ggplot() +
    annotate('rect', xmin = -0.4, xmax = -0.2, ymin = -Inf, ymax = Inf, fill = colors[ 1 ], size = 0.5, alpha = 0.5) +
    annotate('rect', xmin = -0.1, xmax = 0, ymin = -Inf, ymax = Inf, fill = colors[ 2 ], size = 0.5, alpha = 0.5) +
    annotate('rect', xmin = 0, xmax = 0.1, ymin = -Inf, ymax = Inf, fill = colors[ 3 ], size = 0.5, alpha = 0.5) +
    annotate('rect', xmin = 0.2, xmax = 0.4, ymin = -Inf, ymax = Inf, fill = colors[ 4 ], size = 0.5, alpha = 0.5) +
    geom_point(aes(x = Time, y = Trial), color = 'black', size = 0.5, alpha = 0.75) +
    theme_light(base_size = 8) +
    xlab('Time, t (s)') +
    theme(axis.text.y = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(strip.background = element_blank()) +
    theme(strip.text = element_blank()) +
    ylab('Movement event') +
    annotate('text', x = -0.3, y = max_trial + max_trial * 0.08, label = 'B', size = 2.5) +
    annotate('text', x = -0.05, y = max_trial + max_trial * 0.08, label = 'P', size = 2.5) +
    annotate('text', x = 0.05, y = max_trial + max_trial * 0.08, label = 'O', size = 2.5) +
    annotate('text', x = 0.3, y = max_trial + max_trial * 0.08, label = 'L', size = 2.5) +
    theme(panel.grid = element_blank()) +
    expand_limits(y = c(0, max_trial + max_trial * 0.1))

  p5 <- ggdraw() +
    draw_plot(p1, 0, 0.75, 1, 0.25) +
    draw_plot(p2, 0.068, 0.4, 1 - 0.068, 0.35) +
    draw_plot(p3, 0.08, 0.2, 1 - 0.08, 0.2) +
    draw_plot(p4, 0.138, 0, 1 - 0.138, 0.2)

  # Only for the first in list
  if (region == 'S1_L23') {
    p5 <- p5 +
      draw_plot_label(letters[ 1:4 ], x = c(0, 0, 0, 0), y = c(1, 0.75, 0.4, 0.2), size = 9)
  }

  plot[[region]] <- p5
}

final <- ggarrange(plotlist = plot, ncol = 4)

ggsave(final, filename = 'figure3.png', height = 7.4, width = 8, bg = 'white', dpi = 600)
